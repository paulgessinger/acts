from typing import Type
import pytest

import acts
from acts.examples import (
    RootParticleWriter,
    RootParticleReader,
    BareAlgorithm,
    Sequencer,
)


def test_root_particle_reader(tmp_path, ptcl_gun):
    # need to write out some particles first
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.ERROR)
    evGen = ptcl_gun(s)

    file = tmp_path / "particles.root"
    s.addWriter(
        RootParticleWriter(
            acts.logging.ERROR,
            inputParticles=evGen.config.outputParticles,
            filePath=str(file),
        )
    )

    s.run()

    del s  # to properly close the root file

    # reset sequencer for reading

    s2 = Sequencer(numThreads=1, logLevel=acts.logging.ERROR)

    s2.addReader(
        RootParticleReader(
            acts.logging.ERROR,
            particleCollection="input_particles",
            filePath=str(file),
        )
    )

    n = 0

    class CheckAlg(BareAlgorithm):
        def execute(self, ctx):
            assert ctx.eventStore.exists("input_particles")
            nonlocal n
            n += 1
            return acts.examples.ProcessCode.SUCCESS

    s2.addAlgorithm(CheckAlg("check_alg", acts.logging.ERROR))

    s2.run()

    assert n == 10


# def test_csv_particle_writer(tmp_path, conf_const, ptcl_gun):
#     s = Sequencer(numThreads=1, events=10)
#     evGen = ptcl_gun(s)

#     out = tmp_path / "csv"

#     out.mkdir()

#     s.addWriter(
#         conf_const(
#             CsvParticleWriter,
#             acts.logging.INFO,
#             inputParticles=evGen.config.outputParticles,
#             outputStem="particle",
#             outputDir=str(out),
#         )
#     )

#     s.run()

#     assert len([f for f in out.iterdir() if f.is_file()]) == s.config.events
#     assert all(f.stat().st_size > 0 for f in out.iterdir())


# def test_root_meas_writer(tmp_path, fatras, trk_geo):
#     s = Sequencer(numThreads=1, events=10)
#     evGen, simAlg, digiAlg = fatras(s)

#     out = tmp_path / "meas.root"

#     assert not out.exists()

#     config = RootMeasurementWriter.Config(
#         inputMeasurements=digiAlg.config.outputMeasurements,
#         inputClusters=""
#         if digiAlg.config.isSimpleSmearer
#         else digiAlg.config.outputClusters,
#         inputSimHits=simAlg.config.outputSimHits,
#         inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
#         filePath=str(out),
#         trackingGeometry=trk_geo,
#     )
#     config.addBoundIndicesFromDigiConfig(digiAlg.config)
#     s.addWriter(RootMeasurementWriter(level=acts.logging.INFO, config=config))
#     s.run()

#     assert out.exists()
#     assert out.stat().st_size > 0


# def test_csv_meas_writer(tmp_path, fatras, trk_geo):
#     s = Sequencer(numThreads=1, events=10)
#     evGen, simAlg, digiAlg = fatras(s)

#     out = tmp_path / "csv"
#     out.mkdir()

#     config = CsvMeasurementWriter.Config(
#         inputMeasurements=digiAlg.config.outputMeasurements,
#         inputClusters=""
#         if digiAlg.config.isSimpleSmearer
#         else digiAlg.config.outputClusters,
#         inputSimHits=simAlg.config.outputSimHits,
#         inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
#         outputDir=str(out),
#     )
#     s.addWriter(CsvMeasurementWriter(level=acts.logging.INFO, config=config))
#     s.run()

#     assert len([f for f in out.iterdir() if f.is_file()]) == s.config.events * 3
#     assert all(f.stat().st_size > 0 for f in out.iterdir())


# # @pytest.mark.parametrize(
# #     "writer,extra",
# #     [
# #         (RootPropagationStepsWriter, {"collection": "collection"}),
# #         (RootParticleWriter, {"inputParticles": "collection"}),
# #         (
# #             TrackFinderPerformanceWriter,
# #             {
# #                 "inputProtoTracks": "prototracks",
# #                 "inputMeasurementParticlesMap": "measmap",
# #                 "inputParticles": ""
# #             },
# #         ),
# #     ],
# # )
# # def test_root_writer_interface(tmp_path, writer, extra):
# #     with pytest.raises(TypeError):
# #         writer()

# #     f = tmp_path / "test.root"
# #     assert not f.exists()
# #     assert writer(acts.logging.INFO, filePath=str(f), **extra)
# #     assert f.exists()
