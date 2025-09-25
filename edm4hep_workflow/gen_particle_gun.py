#!/usr/bin/env python3

import acts
import acts.examples
import acts.examples.hepmc3

from acts import UnitConstants as u

from edm4hep_workflow.args import make_parser

p = make_parser()
p.add_argument("--type", required=True, choices=["mu", "pi", "el"])
p.add_argument("--pt", type=float, required=True)
p.add_argument("--events", type=int, required=True)

args = p.parse_args()

s = acts.examples.Sequencer(numThreads=args.jobs, events=args.events)

rng = acts.examples.RandomNumbers(seed=42)

pdg = {
    "mu": acts.PdgParticle.eMuon,
    "pi": acts.PdgParticle.ePionPlus,
    "el": acts.PdgParticle.eElectron,
}[args.type]

evGen = acts.examples.EventGenerator(
    level=acts.logging.INFO,
    generators=[
        acts.examples.EventGenerator.Generator(
            multiplicity=acts.examples.FixedMultiplicityGenerator(n=1),
            vertex=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(50 * u.um, 50 * u.um, 150 * u.mm, 20 * u.ns),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            particles=acts.examples.ParametricParticleGenerator(
                p=(args.pt * u.GeV, args.pt * u.GeV),
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
if args.output.suffix != ".hepmc3":
    stem = args.output.stem
    suffix = args.output.suffix

    available = acts.examples.hepmc3.availableCompressionModes()

    for comp in available:
        ext = acts.examples.hepmc3.compressionExtension(comp)
        if suffix == ext:
            compression = comp
            break
    if compression is None:
        raise RuntimeError(f"Invalid output format {args.output}")
    args.output = args.output.parent / stem
else:
    compression = acts.examples.hepmc3.Compression.none


s.addWriter(
    acts.examples.hepmc3.HepMC3Writer(
        acts.logging.INFO,
        inputEvent="hepmc3_event",
        outputPath=args.output,
        compression=compression,
    )
)

s.run()
