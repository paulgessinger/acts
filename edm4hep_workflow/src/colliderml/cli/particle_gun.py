import json
from typing import Annotated
from pathlib import Path
import os
import enum

import typer

from colliderml import logging
from colliderml.constants import SEED_DEFAULT
from colliderml.cli import args


class ParticleType(enum.StrEnum):
    mu = "mu"
    pi = "pi"
    el = "el"

    @property
    def pdg(self) -> "acts.PdgParticle":
        import acts

        return {
            "mu": acts.PdgParticle.eMuon,
            "pi": acts.PdgParticle.ePionPlus,
            "el": acts.PdgParticle.eElectron,
        }[self]


def main(
    output: Path,
    particle_type: Annotated[ParticleType, typer.Option("--type")],
    events: args.EVENTS,
    pt: Annotated[float, typer.Option(help="Transverse momentum in GeV")],
    jobs: args.JOBS = os.cpu_count() or 1,
    seed: int = SEED_DEFAULT,
):

    import acts
    import acts.examples
    import acts.examples.hepmc3

    logger = logging.get_logger(__name__)

    from acts import UnitConstants as u

    s = acts.examples.Sequencer(numThreads=jobs, events=events)

    rng = acts.examples.RandomNumbers(seed=seed)

    pdg = particle_type.pdg

    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=1),
                vertex=acts.examples.GaussianVertexGenerator(
                    # stddev=acts.Vector4(50 * u.um, 50 * u.um, 150 * u.mm, 20 * u.ns),
                    stddev=acts.Vector4(0 * u.um, 0 * u.um, 0 * u.mm, 0 * u.ns),
                    mean=acts.Vector4(0, 0, 0, 0),
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(pt * u.GeV, pt * u.GeV),
                    pTransverse=True,
                    pdg=pdg,
                    eta=(-3, 3),
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

    s.addAlgorithm(
        acts.examples.hepmc3.HepMC3InputConverter(
            level=acts.logging.INFO,
            inputEvent=evGen.config.outputEvent,
            outputParticles="particles_generated",
            outputVertices="vertices_truth",
        )
    )

    compression = None
    meta_file = output.with_suffix(output.suffix + ".json")

    if output.suffix != ".hepmc3":
        stem = output.stem
        suffix = output.suffix

        print(suffix, stem)

        available = acts.examples.hepmc3.availableCompressionModes()

        for comp in available:
            ext = acts.examples.hepmc3.compressionExtension(comp)
            if suffix == ext:
                compression = comp
                break
        if compression is None:
            raise RuntimeError(f"Invalid output format {output}")
        output = output.parent / stem
    else:
        compression = acts.examples.hepmc3.Compression.none

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

    with meta_file.open("w") as mf:
        json_content = {"num_events": events}
        json.dump(json_content, mf, indent=4)
