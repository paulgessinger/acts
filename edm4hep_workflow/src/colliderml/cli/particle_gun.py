import json
from typing import Annotated
from pathlib import Path
import os
import enum

import typer

from colliderml import logging
from colliderml.constants import SEED_DEFAULT
from colliderml.cli import args
from colliderml.util import HepMC3Meta


class ParticleType(enum.StrEnum):
    mu = "mu"
    pi = "pi"
    el = "el"
    g = "g"
    k = "k"

    @property
    def pdg(self) -> "acts.PdgParticle":
        import acts

        return {
            "mu": acts.PdgParticle.eMuon,
            "pi": acts.PdgParticle.ePionPlus,
            "el": acts.PdgParticle.eElectron,
            "g": acts.PdgParticle.eGamma,
            "k": acts.PdgParticle.eKaonPlus,
        }[self]


def main(
    output: Path,
    particle_type: Annotated[ParticleType, typer.Option("--type")],
    events: args.EVENTS,
    pt: Annotated[
        str, typer.Option(help="Transverse momentum in GeV or range min,max")
    ],
    jobs: args.JOBS = os.cpu_count() or 1,
    seed: int = SEED_DEFAULT,
):

    import acts
    import acts.examples
    import acts.examples.hepmc3
    from acts import UnitConstants as u

    logger = logging.get_logger(__name__)

    s = acts.examples.Sequencer(numThreads=jobs, events=events)

    rng = acts.examples.RandomNumbers(seed=seed)

    pdg = particle_type.pdg

    ptvals = pt.split(",")
    if len(ptvals) == 1:
        ptmin = float(ptvals[0])
        ptmax = ptmin
    elif len(ptvals) == 2:
        ptmin = float(ptvals[0])
        ptmax = float(ptvals[1])
    else:
        logger.error("Invalid pt specification: %s", pt)
        raise typer.Exit(1)

    logger.info("Transverse momentum range: %.2f - %.2f GeV", ptmin, ptmax)

    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=1),
                vertex=acts.examples.FixedVertexGenerator(),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(ptmin * u.GeV, ptmax * u.GeV),
                    pTransverse=True,
                    pdg=pdg,
                    eta=(-3, 3),
                    etaUniform=True,
                    phi=(0, 360 * u.degree),
                    randomizeCharge=True,
                    numParticles=1,
                ),
            )
        ],
        outputEvent="hepmc3_event",
        randomNumbers=rng,
    )
    s.addReader(evGen)

    # s.addAlgorithm(
    #     acts.examples.hepmc3.HepMC3InputConverter(
    #         level=acts.logging.INFO,
    #         inputEvent=evGen.config.outputEvent,
    #         outputParticles="particles_generated",
    #         outputVertices="vertices_truth",
    #     )
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
