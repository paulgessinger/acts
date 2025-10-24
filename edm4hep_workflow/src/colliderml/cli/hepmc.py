import os
from colliderml import constants
from colliderml.cli import args
import typer
from pathlib import Path
from typing import Annotated
import enum

import colliderml.logging
from colliderml.config import PileupConfig, PileupStrategy
from colliderml.util import human_readable_size

app = typer.Typer()


@app.command()
def merge(
    hard_scatter: Annotated[
        Path,
        typer.Option("--hard-scatter", "--hs", dir_okay=False, exists=True),
    ],
    output: Path,
    pileup: Annotated[
        Path | None, typer.Option("--pileup", "--pu", dir_okay=False, exists=True)
    ] = None,
    npileup: int = 0,
    nhard_scatter: int = 1,
    seed: int = constants.SEED_DEFAULT,
    jobs: args.JOBS = 1,
    config_path: Annotated[Path | None, typer.Option("--config")] = None,
    max_events: int | None = None,
):
    import acts
    import acts.examples
    import acts.examples.hepmc3

    logger = colliderml.logging.get_logger(__name__)
    logger.info("HepMC3 merging")

    logger.info("Using configuration: %s", config_path)
    config = PileupConfig.load(config_path)
    logger.debug("Config: %s", config)

    compression = acts.examples.hepmc3.compressionFromFilename(output)
    format = acts.examples.hepmc3.formatFromFilename(output)

    vertex_mean = config.vertex_mean.to_acts()
    vertex_stddev = config.vertex_stddev.to_acts()

    logger.info("Hard scatter input: %s", hard_scatter)
    logger.info("Pileup input: %s", pileup)
    logger.info("Merging strategy: %s", config.strategy)
    logger.info("Number of hard scatter events: %d", nhard_scatter)
    logger.info("Number of pileup events: %d", npileup)
    logger.info("Output compression: %s", compression.name)
    logger.info("Output format: %s", format.name)
    logger.info("Output file: %s", output)
    logger.info("Random seed: %d", seed)
    logger.info("Number of jobs: %d", jobs)
    logger.info("Vertex mean: %s -> %s", config.vertex_mean, vertex_mean)
    logger.info("Vertex stddev: %s -> %s", config.vertex_stddev, vertex_stddev)
    if jobs > 1:
        logger.info("With >1 jobs events may be written out of order")

    s = acts.examples.Sequencer(
        numThreads=jobs, logLevel=acts.logging.INFO, events=max_events
    )

    rng = acts.examples.RandomNumbers(seed=seed)

    HepMC3Reader = acts.examples.hepmc3.HepMC3Reader
    inputs = [HepMC3Reader.Input.Fixed(hard_scatter, nhard_scatter)]

    if pileup is not None:
        if npileup <= 0:
            logger.error("Number of pileup events must be > 0 if pileup file is given")
            raise typer.Exit(1)

        if config.strategy == PileupStrategy.fixed:
            inputs.append(HepMC3Reader.Input.Fixed(pileup, npileup))
        else:
            inputs.append(HepMC3Reader.Input.Poisson(pileup, npileup))

    s.addReader(
        HepMC3Reader(
            inputs=inputs,
            level=acts.logging.INFO,
            outputEvent="hepmc3_event",
            checkEventNumber=False,  # This is not generally guaranteed for arbitrary inputs
            randomNumbers=rng,
            vertexGenerator=acts.examples.GaussianVertexGenerator(
                stddev=vertex_stddev,
                mean=vertex_mean,
            ),
        )
    )

    print("OUTPUT", output)

    s.addWriter(
        acts.examples.hepmc3.HepMC3Writer(
            inputEvent="hepmc3_event",
            outputPath=output,
            level=acts.logging.INFO,
            writeEventsInOrder=False,
        )
    )

    s.run()

    # Get output file size and report
    if not output.exists():
        logger.error("Output file %s does not exist", output)
        raise typer.Exit(1)

    output_size = output.stat().st_size
    logger.info("Wrote output to %s (%s)", output, human_readable_size(output_size))

    # Also check for sidecar metadata file
    sidecar_path = Path(str(output) + ".json")
    if sidecar_path.exists():
        logger.info("Sidecar metadata file: %s", sidecar_path)
