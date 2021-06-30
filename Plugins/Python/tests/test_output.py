from typing import Type
import pytest

import acts
from acts.examples import (
    ObjPropagationStepsWriter,
    RootPropagationStepsWriter,
    RootParticleWriter,
    TrackFinderPerformanceWriter,
)


# def test_obj_propagation_step_writer():
#     with pytest.raises(TypeError):
#         ObjPropagationStepsWriter()
#     assert ObjPropagationStepsWriter(acts.logging.INFO, collection="some_collection")


# @pytest.mark.parametrize(
#     "writer,extra",
#     [
#         (RootPropagationStepsWriter, {"collection": "collection"}),
#         (RootParticleWriter, {"inputParticles": "collection"}),
#         (
#             TrackFinderPerformanceWriter,
#             {
#                 "inputProtoTracks": "prototracks",
#                 "inputMeasurementParticlesMap": "measmap",
#                 "inputParticles": ""
#             },
#         ),
#     ],
# )
# def test_root_writer_interface(tmp_path, writer, extra):
#     with pytest.raises(TypeError):
#         writer()

#     f = tmp_path / "test.root"
#     assert not f.exists()
#     assert writer(acts.logging.INFO, filePath=str(f), **extra)
#     assert f.exists()
