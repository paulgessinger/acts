#!/usr/bin/env python3
from pathlib import Path
import typer
import acts


def main(
    input: Path,
    output: Path,
    jobs: int = -1,
    seed: int = 998877,
    logLevel: str = "INFO",
):
    from acts import UnitConstants as u
    import acts.examples
    from acts.examples import Sequencer
    from acts.examples.edm4hep import EDM4hepSimInputConverter
    from acts.examples.podio import PodioReader, PodioWriter
    from acts.examples.odd import getOpenDataDetectorDirectory, getOpenDataDetector
    from acts.examples.simulation import (
        addSimParticleSelection,
        addDigitization,
        addDigiParticleSelection,
        ParticleSelectorConfig,
    )

    print(input, "->", output)

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

    oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"

    oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

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
        outputSimVertices="simvertices",
        dd4hepDetector=detector,
        trackingGeometry=trackingGeometry,
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

    # Add digitization if enabled
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=oddDigiConfig,
        # outputDirRoot=perf_output if getattr(config, "output_root", True) else None,
        outputDirCsv=None,
        rnd=rnd,
        logLevel=logLevel,
    )

    # Add digi particle selection (filters particles with sufficient measurements)
    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(1.0 * u.GeV, None),
            eta=(-3.0, 3.0),
            measurements=(9, None),
            removeNeutral=True,
        ),
        logLevel=logLevel,
    )

    podioWriter = PodioWriter(
        level=logLevel,
        outputPath=str(output),
        inputFrame="events",
        category=podioReader.config.category,
    )
    s.addWriter(podioWriter)

    s.run()


typer.run(main)
