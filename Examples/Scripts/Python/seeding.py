#!/usr/bin/env python3

import acts
import acts.examples

u = acts.UnitConstants

# Preliminaries
rnd = acts.examples.RandomNumbers()

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
    level=acts.logging.INFO,
    generators=[g],
    outputParticles="particles_input",
    randomNumbers=rnd,
)

# Simulation
alg = acts.examples.FatrasAlgorithm(
    level=acts.logging.INFO,
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
root = acts.examples.RootParticleWriter(
    inputParticles="particles_final", filePath="fatras_particles.root"
)

# Sequencer
s = acts.examples.Sequencer(events=100, numThreads=-1, logLevel=acts.logging.INFO)

s.addReader(evGen)
s.addAlgorithm(alg)
# s.addWriter(csv)
# s.addWriter(root)

s.run()


print("DONE")
