from typing import Annotated
from pathlib import Path
import os
import json

from colliderml.util import HepMC3Meta
import typer

from colliderml import logging
from colliderml.constants import SEED_DEFAULT
from colliderml.cli import args

from colliderml.config import Pythia8SampleConfig as SampleConfig


def main(
    sample_file: Annotated[Path, typer.Argument(..., exists=True, dir_okay=False)],
    output: Path,
    events: args.EVENTS,
    jobs: args.JOBS = os.cpu_count() or 1,
    seed: int = SEED_DEFAULT,
):

    config = SampleConfig.load(sample_file)
    logger = logging.get_logger(__name__)
    logger.info("Loaded sample configuration from %s", sample_file)
    logger.info("Sample label: %s", config.label)
    logger.info("Settings: %s", config.settings)
    logger.info("Beam: %s -> 💥 <- %s", config.pdg_beam0, config.pdg_beam1)
    logger.info("Jobs: %d", jobs)

    import acts
    import acts.examples
    import acts.examples.hepmc3
    from acts import UnitConstants as u

    s = acts.examples.Sequencer(numThreads=jobs, events=events)

    rng = acts.examples.RandomNumbers(seed=seed)

    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=1),
                vertex=acts.examples.FixedVertexGenerator(),
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=acts.logging.INFO,
                    pdgBeam0=config.pdg_beam0,
                    pdgBeam1=config.pdg_beam1,
                    cmsEnergy=config.cms_energy,
                    settings=config.settings,
                    # printLongEventListing=printLongEventListing,
                    # printShortEventListing=printShortEventListing,
                ),
            )
        ],
        outputEvent="hepmc3_event",
        randomNumbers=rng,
    )
    s.addReader(evGen)

    # should we do particle selection here? probably not?
    # addGenParticleSelection(
    #     s,
    #     ParticleSelectorConfig(
    #         rho=(0.0, 24 * u.mm),
    #         absZ=(0.0, 1.0 * u.m),
    #         eta=(-3.0, 3.0),
    #         pt=(150 * u.MeV, None),
    #     ),
    # )

    compression = acts.examples.hepmc3.compressionFromFilename(output)
    meta_file = output.with_suffix(output.suffix + ".json")

    if compression is not acts.examples.hepmc3.Compression.none:
        output = output.parent / output.stem

    logger.info("Writing output to %s (compression=%s)", output, compression)

    s.addWriter(
        acts.examples.hepmc3.HepMC3Writer(
            acts.logging.INFO,
            inputEvent="hepmc3_event",
            outputPath=output,
            compression=compression,
        )
    )

    s.run()

    logger.info("Writing meta file to %s", meta_file)

    meta = HepMC3Meta(num_events=events)

    with meta_file.open("w") as mf:
        json.dump(meta.model_dump(), mf, indent=4)
