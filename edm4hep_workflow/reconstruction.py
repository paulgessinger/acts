#!/usr/bin/env python3
from pathlib import Path
import re
import typer
import acts
from typing import Annotated

from config import Config


def main(
    input: Path,
    output: Path,
    jobs: int = -1,
    seed: int = 998877,
    logLevel: str = "INFO",
    config_path: Annotated[Path | None, typer.Option("--config")] = None,
):
    from acts import UnitConstants as u
    import acts.examples
    from acts.examples import Sequencer
    from acts.examples.edm4hep import (
        EDM4hepSimInputConverter,
        PodioReader,
        PodioWriter,
        PodioMeasurementInputConverter,
    )
    from acts.examples.odd import getOpenDataDetectorDirectory, getOpenDataDetector
    from acts.examples.simulation import (
        addSimParticleSelection,
        addDigitization,
        addDigiParticleSelection,
        ParticleSelectorConfig,
    )

    from acts.examples.reconstruction import (
        addSeeding,
        addCKFTracks,
        TrackSelectorConfig,
        CkfConfig,
        SeedingAlgorithm,
        SeedFinderConfigArg,
    )

    print(input, "->", output)

    config = Config.load(config_path)

    logLevel = getattr(acts.logging, logLevel.upper(), acts.logging.INFO)

    rnd = acts.examples.RandomNumbers(seed=seed)

    s = Sequencer(numThreads=jobs)
    s.config.logLevel = acts.logging.DEBUG

    # Get detector and field
    geoDir = getOpenDataDetectorDirectory()

    # Set performance output directory based on flag
    # perf_output = output_dir if getattr(config, "performance_metrics", False) else None

    # Load material map
    oddMaterialMap = geoDir / "data/odd-material-maps.root"

    oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

    oddSeedingSel = geoDir / "config/odd-seeding-config.json"

    # Get detector
    detector = getOpenDataDetector(odd_dir=geoDir, materialDecorator=oddMaterialDeco)
    trackingGeometry = detector.trackingGeometry()
    field = detector.field

    # Configure EDM4hep reader and converter
    # Step 1: PodioReader to read the EDM4hep file
    podioReader = PodioReader(
        level=logLevel,
        inputPath=str(input),
        outputFrame="events",
        category="events",
    )
    s.addReader(podioReader)

    # Step 2: EDM4hepSimInputConverter algorithm to convert EDM4hep data to ACTS format
    edm4hepConverter = EDM4hepSimInputConverter(
        level=logLevel,
        inputFrame="events",
        inputSimHits=[
            "PixelBarrelReadout",
            "PixelEndcapReadout",
            "ShortStripBarrelReadout",
            "ShortStripEndcapReadout",
            "LongStripBarrelReadout",
            "LongStripEndcapReadout",
        ],
        outputParticlesGenerator="particles_input",
        outputParticlesSimulation="particles_simulated",
        outputSimHits="simhits",
        outputSimHitAssociation="simhit_associations",
        outputSimVertices="simvertices",
        sortSimHitsInTime=True,
        dd4hepDetector=detector,
        trackingGeometry=trackingGeometry,
        **config.sim_hit_reading.model_dump(),
    )
    s.addAlgorithm(edm4hepConverter)
    s.addWhiteboardAlias("particles", "particles_input")

    # Add sim particle selection (filters particles from simulation)
    addSimParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-4.0, 4.0),
            pt=(150 * u.MeV, None),
            removeNeutral=True,
        ),
        logLevel=logLevel,
    )

    measConv = PodioMeasurementInputConverter(
        level=logLevel,
        inputFrame=podioReader.config.outputFrame,
        inputMeasurements="ActsMeasurements",
        outputMeasurements="measurements",
        outputMeasurementParticlesMap="measurement_particles_map",
        outputMeasurementSimHitsMap="measurement_simhits_map",
        outputParticleMeasurementsMap="particle_measurements_map",
        outputSimHitMeasurementsMap="simhit_measurements_map",
        inputSimHits=edm4hepConverter.config.outputSimHits,
        inputSimHitAssociation=edm4hepConverter.config.outputSimHitAssociation,
    )
    s.addAlgorithm(measConv)

    # Add digi particle selection (filters particles with sufficient measurements)
    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            # we are only interested in the hard scatter vertex
            # primaryVertexId=(1, 2),
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            # using something close to 1 to include for sure
            pt=(0.999 * u.GeV, None),
            measurements=(6, None),
            removeNeutral=True,
            removeSecondaries=False,
            # nMeasurementsGroupMin=measurementCounter,
        ),
        logLevel=logLevel,
    )

    # @TODO: Add track finding etc here

    # Add seeding
    addSeeding(
        s,
        trackingGeometry,
        field,
        seedingAlgorithm=SeedingAlgorithm.Default,
        particleHypothesis=acts.ParticleHypothesis.pion,
        seedFinderConfigArg=SeedFinderConfigArg(
            r=(33 * u.mm, 200 * u.mm),
            # kills efficiency at |eta|~2
            deltaR=(1 * u.mm, 300 * u.mm),
            collisionRegion=(-250 * u.mm, 250 * u.mm),
            z=(-2000 * u.mm, 2000 * u.mm),
            maxSeedsPerSpM=40,
            sigmaScattering=5,
            radLengthPerSeed=0.1,
            minPt=0.5 * u.GeV,
            impactMax=3 * u.mm,
            zBinEdges=[-1600, -1000, -600, 0, 600, 1000, 1600],
        ),
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0.1 / u.GeV,
            1 * u.ns,
        ],
        initialSigmaQoverPt=0.1 * u.e / u.GeV,
        initialSigmaPtRel=0.1,
        initialVarInflation=[1e0, 1e0, 1e0, 1e0, 1e0, 1e0],
        geoSelectionConfigFile=oddSeedingSel,
        outputDirRoot=None,
        logLevel=logLevel,
    )

    # Add CKF tracking
    addCKFTracks(
        s,
        trackingGeometry,
        field,
        trackSelectorConfig=TrackSelectorConfig(
            pt=(0.7 * u.GeV, None),
            absEta=(None, 3.5),
            nMeasurementsMin=6,
            maxHolesAndOutliers=3,
        ),
        ckfConfig=CkfConfig(
            chi2CutOffMeasurement=15.0,
            chi2CutOffOutlier=25.0,
            numMeasurementsCutOff=1,
            seedDeduplication=True,
            stayOnSeed=True,
        ),
        twoWay=True,
        outputDirRoot=output.parent,
        outputDirCsv=None,
        writeCovMat=None,
        writeTrackStates=None,
        writeTrackSummary=True,
        writePerformance=True,
        logLevel=logLevel,
    )

    podioWriter = PodioWriter(
        level=logLevel,
        outputPath=str(output),
        inputFrame="events",
        category=podioReader.config.category,
        collections=[],
    )
    s.addWriter(podioWriter)

    s.run()

    performance_files = ["performance_finding_ckf", "performance_fitting_ckf"]
    name_base = re.sub(r"\.edm4hep\.root$", "", str(output.name))

    for file in performance_files:
        performance_path = output.parent / f"{file}.root"
        if performance_path.exists():
            performance_target = performance_path.parent / f"{file}_{name_base}.root"
            print(performance_path, "->", performance_target)
            performance_path.rename(performance_target)


typer.run(main)
