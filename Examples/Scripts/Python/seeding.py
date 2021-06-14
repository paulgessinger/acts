#!/usr/bin/env python3


import acts
import acts.examples

u = acts.UnitConstants

# Preliminaries
rnd = acts.examples.RandomNumbers()

detector, trackingGeometry, _ = acts.examples.GenericDetector.create()
# detector, trackingGeometry, _ = acts.examples.DD4hepDetector.create(
#     xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
# )

field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

# Input
evGen = acts.examples.EventGenerator(
    level=acts.logging.INFO,
    generators=[
        acts.examples.EventGenerator.Generator(
            multiplicity=acts.examples.FixedMultiplicityGenerator(n=2),
            vertex=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
            ),
            particles=acts.examples.ParametricParticleGenerator(
                p=(1 * u.GeV, 10 * u.GeV),
                eta=(-2, 2),
                phi=(0, 90 * u.degree),
                randomizeCharge=True,
                numParticles=4,
            ),
        )
    ],
    outputParticles="particles_input",
    randomNumbers=rnd,
)

# Read input from input collection (e.g. Pythia8 output)

# evGen = acts.examples.RootParticleReader(
#     level=acts.logging.INFO,
#     particleCollection="particles_input",
#     inputDir="output",
#     inputFile="pythia8_particles.root",
# )


# Simulation
simAlg = acts.examples.FatrasAlgorithm(
    level=acts.logging.INFO,
    inputParticles=evGen.config.outputParticles,
    # inputParticles=evGen.config.particleCollection,
    outputParticlesInitial="particles_initial",
    outputParticlesFinal="particles_final",
    outputSimHits="simhits",
    randomNumbers=rnd,
    trackingGeometry=trackingGeometry,
    magneticField=field,
    generateHitsOnSensitive=True,
)

# Digitization
digiCfg = acts.examples.DigitizationConfig(
    acts.examples.readDigiConfigFromJson(
        "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
    )
)
digiCfg.trackingGeometry = trackingGeometry
digiCfg.randomNumbers = rnd
digiCfg.inputSimHits = simAlg.config.outputSimHits
digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)

spAlg = acts.examples.SpacePointMaker(
    level=acts.logging.INFO,
    inputSourceLinks=digiCfg.outputSourceLinks,
    inputMeasurements=digiCfg.outputMeasurements,
    outputSpacePoints="spacepoints",
    trackingGeometry=trackingGeometry,
    geometrySelection=acts.examples.readJsonGeometryList(
        "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
    ),
)

seedingAlg = acts.examples.SeedingAlgorithm(
    level=acts.logging.INFO,
    inputSpacePoints=[spAlg.config.outputSpacePoints],
    outputSeeds="seeds",
    outputProtoTracks="prototracks",
)

parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
    level=acts.logging.INFO,
    inputProtoTracks=seedingAlg.config.outputProtoTracks,
    inputSpacePoints=[spAlg.config.outputSpacePoints],
    inputSourceLinks=digiCfg.outputSourceLinks,
    outputTrackParameters="estimatedparameters",
    outputProtoTracks="prototracks_estimated",
    trackingGeometry=trackingGeometry,
    magneticField=field,
)

# Output
tfPerf = acts.examples.TrackFinderPerformanceWriter(
    level=acts.logging.INFO,
    inputProtoTracks=seedingAlg.config.outputProtoTracks,
    inputParticles=simAlg.config.outputParticlesFinal,
    inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
    outputDir="output",
    outputFilename="performance_seeding_trees.root",
)

seedPerf = acts.examples.SeedingPerformanceWriter(
    level=acts.logging.INFO,
    inputProtoTracks=seedingAlg.config.outputProtoTracks,
    inputParticles=simAlg.config.outputParticlesFinal,
    inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
    outputDir="output",
    outputFilename="performance_seeding_trees.root",
)

tpWriter = acts.examples.RootTrackParameterWriter(
    level=acts.logging.INFO,
    inputTrackParameters=parEstimateAlg.config.outputTrackParameters,
    inputProtoTracks=parEstimateAlg.config.outputProtoTracks,
    inputParticles=simAlg.config.outputParticlesFinal,
    inputSimHits=simAlg.config.outputSimHits,
    inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
    inputMeasurementSimHitsMap=digiCfg.outputMeasurementSimHitsMap,
    outputDir="output",
    outputFilename="estimatedparams.root",
    outputTreename="estimatedparams",
)

# Sequencer
s = acts.examples.Sequencer(events=100, numThreads=-1, logLevel=acts.logging.INFO)

s.addReader(evGen)
s.addAlgorithm(simAlg)
s.addAlgorithm(digiAlg)
s.addAlgorithm(spAlg)
s.addAlgorithm(seedingAlg)
s.addAlgorithm(parEstimateAlg)

s.addWriter(
    acts.examples.RootParticleWriter(
        level=acts.logging.INFO,
        inputParticles=evGen.config.outputParticles,
        # inputParticles=evGen.config.particleCollection,
        filePath="output/evgen_particles.root",
    )
)

s.addWriter(
    acts.examples.RootParticleWriter(
        level=acts.logging.INFO,
        inputParticles=simAlg.config.outputParticlesFinal,
        filePath="output/fatras_particles_final.root",
    )
)

s.addWriter(
    acts.examples.RootParticleWriter(
        level=acts.logging.INFO,
        inputParticles=simAlg.config.outputParticlesInitial,
        filePath="output/fatras_particles_initial.root",
    )
)
s.addWriter(tfPerf)
s.addWriter(seedPerf)
s.addWriter(tpWriter)

s.run()


print("DONE")
