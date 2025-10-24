from pathlib import Path
import typer
import os
from typing import Annotated

from colliderml.cli import args
from colliderml.constants import SEED_DEFAULT
from colliderml.config import Config
import colliderml.logging


def main(
    input: Path,
    output: Path,
    jobs: args.JOBS = os.cpu_count() or 1,
    seed: int = SEED_DEFAULT,
    logLevel: str = "INFO",
    skip: int = 0,
    events: Annotated[
        int | None,
        typer.Option(help="Number of events to process. (default = all input events)"),
    ] = None,
    # @TODO: Rethink this configuration here
    config_path: Annotated[Path | None, typer.Option("--config")] = None,
):
    logger = colliderml.logging.get_logger(__name__)
    logger.info("Digitization workflow")
    logger.info("Using configuration: %s", config_path)
    logger.debug("Seed: %d", seed)
    logger.debug("Number of jobs: %d", jobs)
    logger.debug("Input file: %s", input)
    logger.debug("Output file: %s", output)
    logger.debug("Skip events: %d", skip)
    logger.debug("Number of events: %s", events)

    import acts
    from acts import UnitConstants as u
    import acts.examples
    from acts.examples import Sequencer
    from acts.examples.edm4hep import (
        EDM4hepSimInputConverter,
        PodioReader,
        PodioWriter,
        PodioMeasurementOutputConverter,
    )
    from acts.examples.odd import getOpenDataDetectorDirectory, getOpenDataDetector
    from acts.examples.simulation import (
        addSimParticleSelection,
        addDigitization,
        addDigiParticleSelection,
        ParticleSelectorConfig,
    )

    config = Config.load(config_path)
    logger.debug("Config: %s", config)

    logger.debug("%s -> %s", input, output)

    logLevel = getattr(acts.logging, logLevel.upper(), acts.logging.INFO)

    rnd = acts.examples.RandomNumbers(seed=seed)

    s = Sequencer(
        numThreads=jobs,
        skip=skip,
        events=events,
        trackFpes=False,
    )
    s.config.logLevel = acts.logging.DEBUG

    # Get detector and field
    geoDir = getOpenDataDetectorDirectory()

    # Set performance output directory based on flag
    # perf_output = output_dir if getattr(config, "performance_metrics", False) else None

    # Load material map
    oddMaterialMap = geoDir / "data/odd-material-maps.root"

    oddDigiConfig = geoDir / config.digitization.config_file
    if not oddDigiConfig.exists():
        logger.error("Digitization config file %s does not exist", oddDigiConfig)
        raise typer.Exit(1)

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
        outputSimHitAssociation="simhit_associations",
        outputSimVertices="simvertices",
        dd4hepDetector=detector,
        trackingGeometry=trackingGeometry,
        **config.sim_hit_reading.model_dump(),
    )
    s.addAlgorithm(edm4hepConverter)
    s.addWhiteboardAlias("particles", "particles_input")

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

    measConv = PodioMeasurementOutputConverter(
        level=logLevel,
        inputMeasurements="measurements",
        outputMeasurements="ActsMeasurements",
        inputMeasurementSimHitsMap="measurement_simhits_map",
        inputSimHitAssociation=edm4hepConverter.config.outputSimHitAssociation,
    )
    s.addAlgorithm(measConv)

    podioWriter = PodioWriter(
        level=logLevel,
        outputPath=str(output),
        inputFrame="events",
        category=podioReader.config.category,
        collections=measConv.collections,
    )
    s.addWriter(podioWriter)

    s.run()
