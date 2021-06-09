#!/usr/bin/env python3

import acts
import acts.examples

u = acts.UnitConstants

# Preliminaries
rndCfg = acts.examples.RandomNumbers.Config()
rndCfg.seed = 42
rnd = acts.examples.RandomNumbers(rndCfg)

gdc = acts.examples.GenericDetector.Config()
detector = acts.examples.GenericDetector()
trackingGeometry, contextDecorators = detector.finalize(gdc, None)

field = acts.ConstantBField(acts.Vector3(0, 0, 2*u.T))

# Input
vtxGen = acts.examples.GaussianVertexGenerator()
vtxGen.stddev = acts.Vector4(0, 0, 0, 0)

ptclGenCfg = acts.examples.ParametricParticleGenerator.Config()
ptclGen = acts.examples.ParametricParticleGenerator(ptclGenCfg)

evGenCfg = acts.examples.EventGenerator.Config()
g = acts.examples.EventGenerator.Generator()
g.multiplicity = acts.examples.FixedMultiplicityGenerator()
g.vertex = vtxGen
g.particles = ptclGen
evGenCfg.generators = [g]
evGenCfg.outputParticles = "particles_input"
evGenCfg.randomNumbers = rnd
evGen = acts.examples.EventGenerator(evGenCfg, acts.logging.VERBOSE)

# Simulation
algCfg = acts.examples.FatrasAlgorithm.Config()
algCfg.inputParticles = "particles_input"
algCfg.outputParticlesInitial = "particles_initial"
algCfg.outputParticlesFinal = "particles_final"
algCfg.outputSimHits = "simhits"
algCfg.randomNumbers = rnd
algCfg.trackingGeometry = trackingGeometry
algCfg.magneticField = field
alg = acts.examples.FatrasAlgorithm(algCfg, acts.logging.VERBOSE)

# Output
csvCfg = acts.examples.CsvParticleWriter.Config()
csvCfg.outputDir = "csv"
csvCfg.inputParticles = "particles_final"
csvCfg.outputStem = "particles_final"
csv = acts.examples.CsvParticleWriter(csvCfg)

# Sequencer
sCfg = acts.examples.Sequencer.Config()
sCfg.events = 10
sCfg.numThreads = 1
s = acts.examples.Sequencer(sCfg)

s.addReader(evGen)
s.addAlgorithm(alg)
s.addWriter(csv)

s.run()


print("DONE")