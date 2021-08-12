#!/usr/bin/env python3
import os

import acts
import acts.examples


u = acts.UnitConstants


def runSeeding(trackingGeometry, field, outputDir):

    logLevel = acts.logging.WARNING

    # Input
    rnd = acts.examples.RandomNumbers(seed=42)
    evGen = acts.examples.EventGenerator(
        level=logLevel,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=2),
                vertex=acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(1 * u.GeV, 10 * u.GeV),
                    eta=(-2, 2),
                    phi=(0, 360 * u.degree),
                    randomizeCharge=True,
                    numParticles=500,
                ),
            )
        ],
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    # Simulation
    simAlg = acts.examples.FatrasAlgorithm(
        level=logLevel,
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
        ),
        trackingGeometry=trackingGeometry,
        randomNumbers=rnd,
        inputSimHits=simAlg.config.outputSimHits,
    )
    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, logLevel)

    selAlg = acts.examples.TruthSeedSelector(
        level=logLevel,
        ptMin=1 * u.GeV,
        eta=(-2.5, 2.5),
        nHitsMin=9,
        inputParticles=simAlg.config.outputParticlesFinal,
        inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
        outputParticles="particles_selected",
    )

    inputParticles = selAlg.config.outputParticles

    spAlg = acts.examples.SpacePointMaker(
        level=logLevel,
        inputSourceLinks=digiCfg.outputSourceLinks,
        inputMeasurements=digiCfg.outputMeasurements,
        outputSpacePoints="spacepoints",
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.readJsonGeometryList(
            "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
    )

    seedingAlg = acts.examples.SeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=[spAlg.config.outputSpacePoints],
        outputSeeds="seeds",
        outputProtoTracks="prototracks",
    )

    parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
        level=logLevel,
        inputProtoTracks=seedingAlg.config.outputProtoTracks,
        inputSpacePoints=[spAlg.config.outputSpacePoints],
        inputSourceLinks=digiCfg.outputSourceLinks,
        outputTrackParameters="estimatedparameters",
        outputProtoTracks="prototracks_estimated",
        trackingGeometry=trackingGeometry,
        magneticField=field,
    )

    s = acts.examples.Sequencer(
        events=100,
        numThreads=-1,
        logLevel=logLevel,
    )

    s.addReader(evGen)
    s.addAlgorithm(simAlg)
    s.addAlgorithm(digiAlg)
    s.addAlgorithm(selAlg)
    s.addAlgorithm(spAlg)
    s.addAlgorithm(seedingAlg)
    s.addAlgorithm(parEstimateAlg)

    seedPerf = acts.examples.SeedingPerformanceCollector(
        level=logLevel,
        inputProtoTracks=seedingAlg.config.outputProtoTracks,
        inputParticles=inputParticles,
        inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
    )

    s.addWriter(seedPerf)

    s.run()

    return (
        seedPerf.nTotalSeeds,
        seedPerf.nTotalMatchedSeeds,
        seedPerf.nTotalParticles,
        seedPerf.nTotalMatchedParticles,
        seedPerf.nTotalDuplicatedParticles,
    )


if "__main__" == __name__:
    from common import getOpenDataDetector

    # detector, trackingGeometry, _ = getOpenDataDetector()
    detector, trackingGeometry, _ = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    (
        nTotalSeeds,
        nTotalMatchedSeeds,
        nTotalParticles,
        nTotalMatchedParticles,
        nTotalDuplicatedParticles,
    ) = runSeeding(trackingGeometry, field, outputDir=os.getcwd())

    efficiency = nTotalMatchedParticles / nTotalParticles
    fakeRate = (nTotalSeeds - nTotalMatchedSeeds) / nTotalSeeds
    duplicationRate = nTotalDuplicatedParticles / nTotalMatchedParticles
    aveNDuplicatedSeeds = (
        nTotalMatchedSeeds - nTotalMatchedParticles
    ) / nTotalMatchedParticles

    print("nTotalSeeds               = ", nTotalSeeds)
    print("nTotalMatchedSeeds        = ", nTotalMatchedSeeds)
    print("nTotalParticles           = ", nTotalParticles)
    print("nTotalMatchedParticles    = ", nTotalMatchedParticles)
    print("nTotalDuplicatedParticles = ", nTotalDuplicatedParticles)

    print("Efficiency (nMatchedParticles / nAllParticles) = ", efficiency)
    print("Fake rate (nUnMatchedSeeds / nAllSeeds) = ", fakeRate)
    print(
        "Duplication rate (nDuplicatedMatchedParticles / nMatchedParticles) = ",
        duplicationRate,
    )
    print(
        "Average number of duplicated seeds ((nMatchedSeeds - nMatchedParticles) "
        "/ nMatchedParticles) = ",
        aveNDuplicatedSeeds,
    )
