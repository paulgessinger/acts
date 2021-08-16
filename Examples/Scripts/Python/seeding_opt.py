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

    gridConfig = acts.SpacePointGridConfig(
        rMax=100 * u.mm,
        deltaRMax=60 * u.mm,
        zMin=-2000 * u.mm,
        zMax=2000 * u.mm,
        cotThetaMax=7.40627,
        minPt=500 * u.MeV,
        bFieldInZ=1.99724 * u.T,
    )

    seedFilterConfig = acts.SeedFilterConfig(
        deltaRMin=1 * u.mm,
        maxSeedsPerSpM=1,
    )

    seedFinderConfig = acts.SeedfinderConfig(
        rMax=gridConfig.rMax,
        deltaRMin=seedFilterConfig.deltaRMin,
        deltaRMax=gridConfig.deltaRMax,
        collisionRegionMin=-250 * u.mm,
        collisionRegionMax=250 * u.mm,
        zMin=gridConfig.zMin,
        zMax=gridConfig.zMax,
        maxSeedsPerSpM=seedFilterConfig.maxSeedsPerSpM,
        cotThetaMax=gridConfig.cotThetaMax,
        bFieldInZ=gridConfig.bFieldInZ,
        minPt=gridConfig.minPt,
        beamPos=acts.Vector2(0 * u.mm, 0 * u.mm),
        sigmaScattering=50,
        radLengthPerSeed=0.1,
        impactMax=3 * u.mm,
    )

    seedingAlg = acts.examples.SeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=[spAlg.config.outputSpacePoints],
        outputSeeds="seeds",
        outputProtoTracks="prototracks",
        gridConfig=gridConfig,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
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
        events=10,
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


def run_trial(trk_geo, field):
    (
        nTotalSeeds,
        nTotalMatchedSeeds,
        nTotalParticles,
        nTotalMatchedParticles,
        nTotalDuplicatedParticles,
    ) = runSeeding(trackingGeometry, field, outputDir=os.getcwd())

    efficiency = 0
    fakeRate = float("inf")
    duplicationRate = float("inf")
    aveNDuplicatedSeeds = float("inf")

    try:
        efficiency = nTotalMatchedParticles / nTotalParticles
    except:
        pass
    try:
        fakeRate = (nTotalSeeds - nTotalMatchedSeeds) / nTotalSeeds
    except:
        pass
    try:
        duplicationRate = nTotalDuplicatedParticles / nTotalMatchedParticles
    except:
        pass
    try:
        aveNDuplicatedSeeds = (
            nTotalMatchedSeeds - nTotalMatchedParticles
        ) / nTotalMatchedParticles
    except:
        pass

    return (
        nTotalSeeds,
        nTotalMatchedSeeds,
        nTotalParticles,
        nTotalMatchedParticles,
        nTotalDuplicatedParticles,
        efficiency,
        fakeRate,
        duplicationRate,
        aveNDuplicatedSeeds,
    )


from deap import algorithms, base, creator, tools
import numpy
import array
import random


if "__main__" == __name__:
    detector, trackingGeometry, _ = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    toolbox = base.Toolbox()

    def evalOneMax(individual):
        return (sum(individual),)

    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", array.array, typecode="b", fitness=creator.FitnessMax)

    # Attribute generator
    toolbox.register("attr_bool", random.randint, 0, 1)
    # Structure initializers
    toolbox.register(
        "individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, 100
    )
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    toolbox.register("evaluate", evalOneMax)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=3)

    pop = toolbox.population(n=300)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)

    pop, log = algorithms.eaSimple(
        pop,
        toolbox,
        cxpb=0.5,
        mutpb=0.2,
        ngen=40,
        stats=stats,
        halloffame=hof,
        verbose=True,
    )

    # (
    #     nTotalSeeds,
    #     nTotalMatchedSeeds,
    #     nTotalParticles,
    #     nTotalMatchedParticles,
    #     nTotalDuplicatedParticles,
    #     efficiency,
    #     fakeRate,
    #     duplicationRate,
    #     aveNDuplicatedSeeds,
    # ) = run_trial(trackingGeometry, field)

    # print("nTotalSeeds               = ", nTotalSeeds)
    # print("nTotalMatchedSeeds        = ", nTotalMatchedSeeds)
    # print("nTotalParticles           = ", nTotalParticles)
    # print("nTotalMatchedParticles    = ", nTotalMatchedParticles)
    # print("nTotalDuplicatedParticles = ", nTotalDuplicatedParticles)

    # print("Efficiency (nMatchedParticles / nAllParticles) = ", efficiency)
    # print("Fake rate (nUnMatchedSeeds / nAllSeeds) = ", fakeRate)
    # print(
    #     "Duplication rate (nDuplicatedMatchedParticles / nMatchedParticles) = ",
    #     duplicationRate,
    # )
    # print(
    #     "Average number of duplicated seeds ((nMatchedSeeds - nMatchedParticles) "
    #     "/ nMatchedParticles) = ",
    #     aveNDuplicatedSeeds,
    # )
