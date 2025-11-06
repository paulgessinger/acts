import os
import pytest
from pathlib import Path

from helpers import (
    edm4hepEnabled,
    podioEnabled,
    AssertCollectionExistsAlg,
)

import acts
from acts import UnitConstants as u

from acts.examples import (
    Sequencer,
)


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_podio_measurement_output_converter(tmp_path, fatras):
    """Smoke test for PodioMeasurementOutputConverter"""
    from acts.examples.edm4hep import (
        PodioMeasurementOutputConverter,
        PodioWriter,
    )

    s = Sequencer(numThreads=1, events=10)
    _, simAlg, digiAlg = fatras(s)

    out = tmp_path / "measurements_podio.root"

    # Test converter instantiation
    converter = PodioMeasurementOutputConverter(
        level=acts.logging.INFO,
        inputMeasurements=digiAlg.config.outputMeasurements,
        outputMeasurements="ActsMeasurements",
    )
    s.addAlgorithm(converter)

    s.addWriter(
        PodioWriter(
            level=acts.logging.INFO,
            outputPath=str(out),
            category="events",
            collections=converter.collections,
        )
    )

    s.run()

    assert os.path.isfile(out)
    assert os.stat(out).st_size > 10


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_podio_measurement_input_converter(tmp_path, fatras):
    """Smoke test for PodioMeasurementInputConverter"""
    from acts.examples.edm4hep import (
        PodioMeasurementOutputConverter,
        PodioMeasurementInputConverter,
        PodioWriter,
        PodioReader,
    )

    s = Sequencer(numThreads=1, events=10)
    _, simAlg, digiAlg = fatras(s)

    out = tmp_path / "measurements_podio.root"

    # Write data first
    converter_out = PodioMeasurementOutputConverter(
        level=acts.logging.INFO,
        inputMeasurements=digiAlg.config.outputMeasurements,
        outputMeasurements="ActsMeasurements",
    )
    s.addAlgorithm(converter_out)

    s.addWriter(
        PodioWriter(
            level=acts.logging.INFO,
            outputPath=str(out),
            category="events",
            collections=converter_out.collections,
        )
    )

    s.run()

    # Read back
    s2 = Sequencer(numThreads=1)

    s2.addReader(
        PodioReader(
            level=acts.logging.INFO,
            inputPath=str(out),
            outputFrame="events",
            category="events",
        )
    )

    converter_in = PodioMeasurementInputConverter(
        level=acts.logging.INFO,
        inputFrame="events",
        inputMeasurements="ActsMeasurements",
        outputMeasurements="measurements_read",
    )
    s2.addAlgorithm(converter_in)

    alg = AssertCollectionExistsAlg("measurements_read", "check_alg", acts.logging.INFO)
    s2.addAlgorithm(alg)

    s2.run()

    assert alg.events_seen == 10


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_podio_measurement_converter_config():
    """Smoke test to verify converter configuration properties"""
    from acts.examples.edm4hep import (
        PodioMeasurementOutputConverter,
        PodioMeasurementInputConverter,
    )

    # Test that we can create converters with minimal config
    # This tests that Python bindings are properly exposed

    # Output converter
    out_config = PodioMeasurementOutputConverter.Config()
    out_config.inputMeasurements = "test_input"
    out_config.outputMeasurements = "test_output"

    converter_out = PodioMeasurementOutputConverter(
        level=acts.logging.INFO,
        **vars(out_config)
    )

    assert converter_out.config.inputMeasurements == "test_input"
    assert converter_out.config.outputMeasurements == "test_output"

    # Input converter
    in_config = PodioMeasurementInputConverter.Config()
    in_config.inputFrame = "test_frame"
    in_config.inputMeasurements = "test_input"
    in_config.outputMeasurements = "test_output"

    converter_in = PodioMeasurementInputConverter(
        level=acts.logging.INFO,
        **vars(in_config)
    )

    assert converter_in.config.inputFrame == "test_frame"
    assert converter_in.config.inputMeasurements == "test_input"
    assert converter_in.config.outputMeasurements == "test_output"
