#!/usr/bin/env python3

import acts
import acts.examples

u = acts.UnitConstants

# Preliminaries
rnd = acts.examples.RandomNumbers()

nhard = 1
npileup = 200
beam0 = acts.PdgParticle.eProton
beam1 = acts.PdgParticle.eProton
cmsEnergy = 14 * u.TeV

vertexGenerator = acts.examples.GaussianVertexGenerator(
    stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
)

generators = []
if nhard > 0:
    generators.append(
        acts.examples.EventGenerator.Generator(
            multiplicity=acts.examples.FixedMultiplicityGenerator(n=nhard),
            vertex=vertexGenerator,
            particles=acts.examples.Pythia8Generator(
                level=acts.logging.INFO,
                pdgBeam0=beam0,
                pdgBeam1=beam1,
                cmsEnergy=cmsEnergy,
                settings=["HardQCD:all = on"],
            ),
        )
    )
if npileup > 0:
    generators.append(
        acts.examples.EventGenerator.Generator(
            multiplicity=acts.examples.FixedMultiplicityGenerator(n=npileup),
            vertex=vertexGenerator,
            particles=acts.examples.Pythia8Generator(
                level=acts.logging.INFO,
                pdgBeam0=beam0,
                pdgBeam1=beam1,
                cmsEnergy=cmsEnergy,
                settings=["SoftQCD:all = on"],
            ),
        )
    )


# Input
evGen = acts.examples.EventGenerator(
    level=acts.logging.INFO,
    generators=generators,
    outputParticles="particles_input",
    randomNumbers=rnd,
)

# Sequencer
s = acts.examples.Sequencer(events=1000, numThreads=-1, logLevel=acts.logging.INFO)

s.addReader(evGen)
# s.addAlgorithm(
#     acts.examples.ParticlesPrinter(
#         level=acts.logging.INFO, inputParticles="particles_input"
#     )
# )
s.addWriter(
    acts.examples.RootParticleWriter(
        level=acts.logging.INFO,
        inputParticles=evGen.config.outputParticles,
        filePath="output/pythia8_particles.root",
    )
)


s.run()


print("DONE")
