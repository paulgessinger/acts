#!/usr/bin/env python3
from pathlib import Path

import acts
from acts.examples import (
    readDigiConfigFromJson,
    DigitizationConfigurator,
    writeDigiConfigToJson,
    GenericDetectorFactory,
    DigiConfigContainer,
)


u = acts.UnitConstants


def runDigitizationConfig(
    trackingGeometry,
    input: Path,
    output: Path,
):
    inputConfig = readDigiConfigFromJson(str(input))

    digiConfigurator = DigitizationConfigurator()
    digiConfigurator.compactify = True
    digiConfigurator.inputDigiComponents = inputConfig

    trackingGeometry.visitSurfaces(digiConfigurator)

    outputConfig = DigiConfigContainer(digiConfigurator.outputDigiComponents)

    writeDigiConfigToJson(outputConfig, str(output))


if "__main__" == __name__:
    detector = GenericDetectorFactory().buildDetector()
    trackingGeometry = detector.gen1Geometry()

    runDigitizationConfig(
        trackingGeometry=trackingGeometry,
        input=Path(__file__).parent
        / "../../Algorithms/Digitization/share/default-smearing-config-generic.json",
        output=Path.cwd() / "digi-config-out.json",
    )
