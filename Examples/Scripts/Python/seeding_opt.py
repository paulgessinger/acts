#!/usr/bin/env python3
import os
import dataclasses
import functools

import acts
import acts.examples


u = acts.UnitConstants


def runSeeding(trackingGeometry, field, gridConfig, seedFilterConfig, seedFinderConfig):

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
    simAlg = acts.examples.FatrasSimulation(
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


class Attribute:
    def __init__(self, _type, min, max):
        assert isinstance(min, _type) and isinstance(max, _type)
        assert type(min) == type(max), f"{type(min)}, {type(max)}"
        self.min = min
        self.max = max


limits = {"maxSeedsPerSpM": Attribute(int, 1, 10)}


@dataclasses.dataclass
class SeedingConfig:
    maxSeedsPerSpM: int = 1
    rMax: float = 200 * u.mm
    deltaRMin: float = 1 * u.mm
    deltaRMax: float = 60 * u.mm
    sigmaScattering: float = 10
    impactMax: float = 3 * u.mm
    minPt: float = 500 * u.MeV


def run_trial(trk_geo, field, config: SeedingConfig):

    gridConfig = acts.SpacePointGridConfig(
        bFieldInZ=1.99724 * u.T,
        minPt=config.minPt,
        rMax=config.rMax,
        zMax=2000 * u.mm,
        zMin=-2000 * u.mm,
        deltaRMax=config.deltaRMax,
        cotThetaMax=7.40627,  # 2.7 eta
    )

    seedFilterConfig = acts.SeedFilterConfig(
        maxSeedsPerSpM=config.maxSeedsPerSpM, deltaRMin=config.deltaRMin
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
        sigmaScattering=config.sigmaScattering,
        radLengthPerSeed=0.1,
        minPt=gridConfig.minPt,
        bFieldInZ=gridConfig.bFieldInZ,
        beamPos=acts.Vector2(0 * u.mm, 0 * u.mm),
        impactMax=config.impactMax,
    )

    (
        nTotalSeeds,
        nTotalMatchedSeeds,
        nTotalParticles,
        nTotalMatchedParticles,
        nTotalDuplicatedParticles,
    ) = runSeeding(
        trackingGeometry,
        field,
        gridConfig=gridConfig,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
    )

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


if "__main__" == __name__:
    detector, trackingGeometry, _ = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    from orion.client import build_experiment

    storage = {
        "database": {
            "type": "ephemeraldb",
        },
    }

    space = {
        "maxSeedsPerSpM": "uniform(1, 10)",
        "rMax": "uniform(100, 200)",
        "deltaRMax": "uniform(10, 100)",
        "deltaRMin": "uniform(1, 10)",
    }

    experiment = build_experiment(
        "random-exp",
        space=space,
        storage=storage,
        # algorithms={"tpe": {"n_initial_points": 5}},
    )

    def evaluate(maxSeedsPerSpM, rMax, deltaRMax, deltaRMin):
        print("evaluate")

        if deltaRMax < deltaRMin:
            deltaRMin, deltaRMax = deltaRMax, deltaRMin

        config = SeedingConfig(
            maxSeedsPerSpM=round(maxSeedsPerSpM),
            rMax=rMax * u.mm,
            deltaRMax=deltaRMax * u.mm,
            deltaRMin=deltaRMin * u.mm,
        )
        (
            nTotalSeeds,
            nTotalMatchedSeeds,
            nTotalParticles,
            nTotalMatchedParticles,
            nTotalDuplicatedParticles,
            efficiency,
            fakeRate,
            duplicationRate,
            aveNDuplicatedSeeds,
        ) = run_trial(trackingGeometry, field, config=config)

        K = 1000
        effScore = efficiency - 10 * (fakeRate * duplicationRate) / K
        print(efficiency, fakeRate, duplicationRate, effScore)

        objective = 1 - effScore

        return [
            {"name": "score", "type": "objective", "value": objective},
        ]


    print("begin workon")
    experiment.workon(evaluate, max_trials=50)
    print("workon done")

    experiment.plot.regret().show()

    df = experiment.to_pandas()

    best = df.iloc[df.objective.idxmin()]
    print(best)
