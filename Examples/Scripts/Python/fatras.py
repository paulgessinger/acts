#!/usr/bin/env python3

import acts
import acts.examples

u = acts.UnitConstants

# Preliminaries
rnd = acts.examples.RandomNumbers(seed=42)

gdc = acts.examples.GenericDetector.Config()
detector = acts.examples.GenericDetector()
trackingGeometry, contextDecorators = detector.finalize(gdc, None)

field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

# Input
vtxGen = acts.examples.GaussianVertexGenerator()
vtxGen.stddev = acts.Vector4(0, 0, 0, 0)

ptclGen = acts.examples.ParametricParticleGenerator(
    p=(1 * u.GeV, 10 * u.GeV), eta=(-2, 2)
)

g = acts.examples.EventGenerator.Generator()
g.multiplicity = acts.examples.FixedMultiplicityGenerator()
g.vertex = vtxGen
g.particles = ptclGen

evGen = acts.examples.EventGenerator(
    level=acts.logging.VERBOSE,
    generators=[g],
    outputParticles="particles_input",
    randomNumbers=rnd,
)

# Simulation
alg = acts.examples.FatrasAlgorithm(
    level=acts.logging.VERBOSE,
    inputParticles="particles_input",
    outputParticlesInitial="particles_initial",
    outputParticlesFinal="particles_final",
    outputSimHits="simhits",
    randomNumbers=rnd,
    trackingGeometry=trackingGeometry,
    magneticField=field,
)

# Output
csv = acts.examples.CsvParticleWriter(
    outputDir="csv", inputParticles="particles_final", outputStem="particles_final"
)

# Sequencer
s = acts.examples.Sequencer(events=10, numThreads=1)

s.addReader(evGen)
s.addAlgorithm(alg)
s.addWriter(csv)

s.run()


print("DONE")
