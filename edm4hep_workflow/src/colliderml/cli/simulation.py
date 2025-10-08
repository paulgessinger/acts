from pathlib import Path
import time
import json
from typing import Annotated
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import math
import os
import sys
import logging
import tempfile
import contextlib

from colliderml.config import Config, SimulationConfig
from colliderml.util import HepMC3Meta, hadd, hepmc_normalize, parse_hepmc3_file
import typer

from colliderml import constants
from colliderml.cli import args
import colliderml.logging


class DDSimException(RuntimeError):
    def __init__(self, message: str, exit_code: int, output_file: str | None = None):
        super().__init__(message)
        self.message = message
        self.exit_code = exit_code
        self.output_file = output_file


def do_simulation(
    *,
    input_files: list[Path],
    output_file: Path,
    seed: int,
    events: int,
    config: SimulationConfig,
) -> int:
    colliderml.logging.configure_logging(logging.INFO)
    logger = colliderml.logging.get_logger(__name__)

    logger.info(f"Simulating and writing to {output_file} from worker")

    from DDSim.DD4hepSimulation import DD4hepSimulation
    from acts.examples.odd import getOpenDataDetectorDirectory
    from g4units import GeV

    odd_dir = getOpenDataDetectorDirectory()
    logger.info(f"Using ODD from {odd_dir}")
    odd_xml = odd_dir / "xml" / "OpenDataDetector.xml"

    ddsim = DD4hepSimulation()

    # Configure DD4hep
    if isinstance(ddsim.compactFile, list):
        ddsim.compactFile = [str(odd_xml)]
    else:
        ddsim.compactFile = str(odd_xml)

    # @TODO: Deal with custom particle handler
    # https://github.com/OpenDataDetector/ColliderML/blob/75ad4313a7f2b1bf86ea140393d3fd9c348a0fba/scripts/simulation/ddsim_run.py#L175
    # ddsim.part.userParticleHandler = "Geant4TCUserParticleHandler"
    # Geant4FullTruthParticleHandler that's the one Daniel wrote
    ddsim.part.userParticleHandler = "Geant4FullTruthParticleHandler"
    ddsim.part.keepAllParticles = False

    # truncate calo particles ?
    ddsim.enableDetailedShowerMode = False

    logger.info(f"Simulating all events from {len(input_files)} files to {output_file}")

    ddsim.part.minimalKineticEnergy = config.minimal_kinetic_energy

    # ddsim.inputFiles = [str(input_file.resolve())]
    ddsim.inputFiles = [str(p.resolve()) for p in input_files]
    ddsim.numberOfEvents = events
    # ddsim.skipNEvents = skip

    ddsim.outputFile = str(output_file.resolve())

    ddsim.random.seed = seed

    ec = ddsim.run()

    return ec


def job_wrapper(*, output_file: Path, **kwargs):
    job_out = output_file.parent / f"{output_file.name}.log"
    if job_out.exists():
        job_out.unlink()

    with job_out.open("w") as log_file:
        os.dup2(log_file.fileno(), sys.stdout.fileno())
        os.dup2(log_file.fileno(), sys.stderr.fileno())
        return do_simulation(output_file=output_file, **kwargs), output_file


def main(
    input_files: Annotated[list[Path], typer.Argument(metavar="INPUT")],
    output: Annotated[Path, typer.Option("--output", "-o")],
    seed: int = constants.SEED_DEFAULT,
    events: Annotated[
        int | None,
        typer.Option(
            "--events",
            "-n",
            help="Number of events to process. (default = all input events)",
        ),
    ] = None,
    procs: Annotated[int, typer.Option("--processes", "-p")] = 1,
    resplit_files: bool = True,
    scratch_dir: Path | None = None,
    force: Annotated[bool, typer.Option("--force", "-f")] = False,
    config_path: Annotated[Path | None, typer.Option("--config")] = None,
):
    logger = colliderml.logging.get_logger(__name__)

    if config_path is None:
        logger.warning("Using default configuration")
        config = Config()
    else:
        logger.info("Loading configuration from %s", config_path)
        config = Config.load(config_path)

    logger.info("Starting simulation with base seed %s", seed)
    logger.info("Simulation config: %s", config.simulation)

    output.parent.mkdir(parents=True, exist_ok=True)

    if len(input_files) == 1:
        logger.info("Single input file provided, limiting events is supported")
    else:
        logger.info("Multiple input files provided, limiting events is not supported")

        if events is not None and not resplit_files:
            logger.error(
                "Number of events cannot be limited with multiple input files without resplitting"
            )
            raise typer.Exit(1)

    if output.exists():
        if not force:
            logger.error(f"Output file {output} exists, use --force to overwrite")
            raise typer.Exit(1)

        output.unlink()

    logger.info("Split processing into %d processes", procs)

    # special case single process case so we can get direct output
    if procs == 1:
        total_events = sum(HepMC3Meta.for_file(f).num_events for f in input_files)

        do_simulation(
            input_files=input_files,
            output_file=output,
            seed=seed,
            events=events or total_events,
            config=config.simulation,
        )
    else:
        total_events = sum(HepMC3Meta.for_file(f).num_events for f in input_files)

        mp_context = multiprocessing.get_context("spawn")

        with ProcessPoolExecutor(
            max_workers=procs, mp_context=mp_context
        ) as executor, contextlib.ExitStack() as stack:
            if scratch_dir is None:
                tempdir = Path(stack.enter_context(tempfile.TemporaryDirectory()))
            else:
                scratch_dir.mkdir(parents=True, exist_ok=True)
                tempdir = scratch_dir

            logger.debug(f"Using temporary directory {tempdir}")

            if len(input_files) != procs and resplit_files:
                logger.info("Resplitting input files to match number of processes")

                effective_events = events or total_events
                events_per_proc = math.ceil(effective_events / procs)

                input_files = hepmc_normalize(
                    files=input_files,
                    # write uncompressed so that we're quick
                    output=tempdir / "split_input.hepmc3",
                    max_events=events,
                    events_per_file=events_per_proc,
                )

                logger.debug(f"Resplit into {len(input_files)} files")

            futures = []
            for i, file in enumerate(input_files):
                file = file.resolve()
                meta = HepMC3Meta.for_file(file)
                parsed = parse_hepmc3_file(file)
                job_output = tempdir / f"{parsed.prefix}.edm4hep.root"
                job_seed = seed + i

                logger.debug(
                    "%s -> %s (%d events, seed: %d)",
                    file,
                    job_output,
                    meta.num_events,
                    job_seed,
                )

                futures.append(
                    executor.submit(
                        job_wrapper,
                        input_files=[file],
                        output_file=job_output,
                        seed=job_seed,
                        events=meta.num_events,
                        config=config.simulation,
                    )
                )

            for f in as_completed(futures):
                try:
                    ec, out_file = f.result()
                    if ec == 0:
                        logger.info(f"Job completed successfully, output at {out_file}")
                    else:
                        logger.error(f"Job failed with exit code {ec}, see {out_file}")
                except Exception as e:
                    logger.exception("Job failed with exception: %s", e)

            output_files = []
            for f in futures:
                _, out_file = f.result()
                output_files.append(out_file)

            logger.info("Merging output files to %s", output)
            hadd(output_files, output)

            file_size = output.stat().st_size / (1024**2)
            logger.info(f"Wrote output file {output} ({file_size:.1f} MB)")
