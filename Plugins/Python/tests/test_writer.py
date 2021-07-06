from typing import Type
import inspect

import pytest


import acts
from acts.examples import (
    ObjPropagationStepsWriter,
    TrackFinderPerformanceWriter,
    SeedingPerformanceWriter,
    RootPropagationStepsWriter,
    RootParticleWriter,
    RootTrackParameterWriter,
    RootMaterialTrackWriter,
    RootMaterialWriter,
    RootPlanarClusterWriter,
    RootSimHitWriter,
    RootTrajectoryStatesWriter,
    RootTrajectorySummaryWriter,
    RootVertexPerformanceWriter,
    RootMeasurementWriter,
    CsvParticleWriter,
    CsvPlanarClusterWriter,
    CsvSimHitWriter,
    CsvMultiTrajectoryWriter,
    CsvTrackingGeometryWriter,
    CsvMeasurementWriter,
    TrackParamsEstimationAlgorithm,
    Sequencer,
)


@pytest.mark.obj
def test_obj_propagation_step_writer(tmp_path, trk_geo, conf_const, basic_prop_seq):
    with pytest.raises(TypeError):
        ObjPropagationStepsWriter()

    obj = tmp_path / "obj"
    obj.mkdir()

    s, alg = basic_prop_seq(trk_geo)
    w = conf_const(
        ObjPropagationStepsWriter,
        acts.logging.INFO,
        collection=alg.config.propagationStepCollection,
        outputDir=str(obj),
    )

    s.addWriter(w)

    s.run()

    assert len([f for f in obj.iterdir() if f.is_file()]) == s.config.events
    for f in obj.iterdir():
        assert f.stat().st_size > 0


@pytest.mark.csv
def test_csv_particle_writer(tmp_path, conf_const, ptcl_gun):
    s = Sequencer(numThreads=1, events=10)
    evGen = ptcl_gun(s)

    out = tmp_path / "csv"

    out.mkdir()

    s.addWriter(
        conf_const(
            CsvParticleWriter,
            acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            outputStem="particle",
            outputDir=str(out),
        )
    )

    s.run()

    assert len([f for f in out.iterdir() if f.is_file()]) == s.config.events
    assert all(f.stat().st_size > 0 for f in out.iterdir())


@pytest.mark.root
def test_root_prop_step_writer(tmp_path, trk_geo, conf_const, basic_prop_seq):
    with pytest.raises(TypeError):
        RootPropagationStepsWriter()

    file = tmp_path / "prop_steps.root"
    assert not file.exists()

    s, alg = basic_prop_seq(trk_geo)
    w = conf_const(
        RootPropagationStepsWriter,
        acts.logging.INFO,
        collection=alg.config.propagationStepCollection,
        filePath=str(file),
    )

    s.addWriter(w)

    s.run()

    assert file.exists()
    assert file.stat().st_size > 0


@pytest.mark.root
def test_root_particle_writer(tmp_path, conf_const, ptcl_gun):
    s = Sequencer(numThreads=1, events=10)
    evGen = ptcl_gun(s)

    file = tmp_path / "particles.root"

    assert not file.exists()

    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            filePath=str(file),
        )
    )

    s.run()

    assert file.exists()
    assert file.stat().st_size > 0


@pytest.mark.root
def test_root_meas_writer(tmp_path, fatras, trk_geo):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "meas.root"

    assert not out.exists()

    config = RootMeasurementWriter.Config(
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputClusters=""
        if digiAlg.config.isSimpleSmearer
        else digiAlg.config.outputClusters,
        inputSimHits=simAlg.config.outputSimHits,
        inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
        filePath=str(out),
        trackingGeometry=trk_geo,
    )
    config.addBoundIndicesFromDigiConfig(digiAlg.config)
    s.addWriter(RootMeasurementWriter(level=acts.logging.INFO, config=config))
    s.run()

    assert out.exists()
    assert out.stat().st_size > 0


@pytest.mark.csv
def test_csv_meas_writer(tmp_path, fatras, trk_geo):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "csv"
    out.mkdir()

    config = CsvMeasurementWriter.Config(
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputClusters=""
        if digiAlg.config.isSimpleSmearer
        else digiAlg.config.outputClusters,
        inputSimHits=simAlg.config.outputSimHits,
        inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
        outputDir=str(out),
    )
    s.addWriter(CsvMeasurementWriter(level=acts.logging.INFO, config=config))
    s.run()

    assert len([f for f in out.iterdir() if f.is_file()]) == s.config.events * 3
    assert all(f.stat().st_size > 0 for f in out.iterdir())


@pytest.mark.parametrize(
    "writer",
    [
        RootPropagationStepsWriter,
        RootParticleWriter,
        TrackFinderPerformanceWriter,
        SeedingPerformanceWriter,
        RootTrackParameterWriter,
        RootMaterialTrackWriter,
        RootMeasurementWriter,
        RootMaterialWriter,
        RootPlanarClusterWriter,
        RootSimHitWriter,
        RootTrajectoryStatesWriter,
        RootTrajectorySummaryWriter,
        RootVertexPerformanceWriter,
        SeedingPerformanceWriter,
    ],
)
@pytest.mark.root
def test_root_writer_interface(writer, conf_const, tmp_path, trk_geo):
    assert hasattr(writer, "Config")

    config = writer.Config

    assert hasattr(config, "filePath")
    assert hasattr(config, "fileMode")

    f = tmp_path / "target.root"
    assert not f.exists()

    kw = {"level": acts.logging.INFO, "filePath": str(f)}

    for k, _ in inspect.getmembers(config):
        if k.startswith("input"):
            kw[k] = "collection"
        if k == "trackingGeometry":
            kw[k] = trk_geo

    assert conf_const(writer, **kw)

    assert f.exists()


@pytest.mark.parametrize(
    "writer",
    [
        CsvParticleWriter,
        CsvMeasurementWriter,
        CsvPlanarClusterWriter,
        CsvSimHitWriter,
        CsvMultiTrajectoryWriter,
        CsvTrackingGeometryWriter,
    ],
)
@pytest.mark.csv
def test_csv_writer_interface(writer, conf_const, tmp_path, trk_geo):
    assert hasattr(writer, "Config")

    config = writer.Config

    assert hasattr(config, "outputDir")

    kw = {"level": acts.logging.INFO, "outputDir": str(tmp_path)}

    for k, _ in inspect.getmembers(config):
        if k.startswith("input"):
            kw[k] = "collection"
        if k == "trackingGeometry":
            kw[k] = trk_geo
        if k == "outputStem":
            kw[k] = "stem"

    assert conf_const(writer, **kw)
