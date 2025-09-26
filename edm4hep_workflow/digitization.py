#!/usr/bin/env python3
from pathlib import Path
import typer
import acts
from typing import Annotated

from edm4hep_workflow.config import Config


def main(
    input: Path,
    output: Path,
    jobs: int = -1,
    seed: int = 998877,
    logLevel: str = "INFO",
    skip: int = 0,
    events: int | None = None,
    config_path: Annotated[Path | None, typer.Option("--config")] = None,
):
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

    print(input, "->", output)

    logLevel = getattr(acts.logging, logLevel.upper(), acts.logging.INFO)

    rnd = acts.examples.RandomNumbers(seed=seed)

    s = Sequencer(numThreads=jobs, skip=skip, events=events)
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
        outputSimHitAssociation="simhit_associations",
        outputSimVertices="simvertices",
        dd4hepDetector=detector,
        trackingGeometry=trackingGeometry,
        **config.sim_hit_reading.model_dump(),
    )
    s.addAlgorithm(edm4hepConverter)
    s.addWhiteboardAlias("particles", "particles_input")

    # Do we digitize all particles? Otherwise, we don't need this yet
    # addSimParticleSelection(
    #     s,
    #     ParticleSelectorConfig(
    #         rho=(0.0, 24 * u.mm),
    #         absZ=(0.0, 1.0 * u.m),
    #         eta=(-4.0, 4.0),
    #         pt=(150 * u.MeV, None),
    #         removeNeutral=True,
    #     ),
    #     logLevel=logLevel,
    # )

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

    # Removed since reconstruction is done later
    # addDigiParticleSelection(
    #     s,
    #     ParticleSelectorConfig(
    #         # we are only interested in the hard scatter vertex
    #         # primaryVertexId=(1, 2),
    #         rho=(0.0, 24 * u.mm),
    #         absZ=(0.0, 1.0 * u.m),
    #         eta=(-3.0, 3.0),
    #         # using something close to 1 to include for sure
    #         pt=(0.999 * u.GeV, None),
    #         measurements=(6, None),
    #         removeNeutral=True,
    #         removeSecondaries=False,
    #         # nMeasurementsGroupMin=measurementCounter,
    #     ),
    #     logLevel=logLevel,
    # )

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


typer.run(main)
