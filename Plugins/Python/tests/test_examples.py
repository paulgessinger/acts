from Examples.Scripts.Python.geantino_recording import runGeantinoRecording
from pathlib import Path

import pytest

import acts
from acts.examples import Sequencer

u = acts.UnitConstants


@pytest.fixture
def field():
    return acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))


@pytest.fixture
def seq():
    return Sequencer(events=10, numThreads=1)


def assert_csv_output(csv_path, stem):
    # print(list(csv_path.iterdir()))
    assert len([f for f in csv_path.iterdir() if f.name.endswith(stem + ".csv")]) > 0
    assert all([f.stat().st_size > 0 for f in csv_path.iterdir()])


def test_fatras(trk_geo, tmp_path, field, seq):
    from fatras import runFatras

    csv = tmp_path / "csv"
    csv.mkdir()

    root_files = [
        "fatras_particles_final.root",
        "fatras_particles_initial.root",
        "hits.root",
    ]

    assert len(list(csv.iterdir())) == 0
    for rf in root_files:
        assert not (tmp_path / rf).exists()

    runFatras(trk_geo, field, str(tmp_path), s=seq).run()

    assert_csv_output(csv, "particles_final")
    assert_csv_output(csv, "particles_initial")
    assert_csv_output(csv, "hits")
    for rf in root_files:
        assert (tmp_path / rf).exists()
        assert (tmp_path / rf).stat().st_size > 0


def test_seeding(tmp_path, trk_geo, field, seq):
    from seeding import runSeeding

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    csv = tmp_path / "csv"
    csv.mkdir()

    root_files = [
        "estimatedparams.root",
        "performance_seeding_trees.root",
        "performance_seeding_hists.root",
        "evgen_particles.root",
        "fatras_particles_final.root",
        "fatras_particles_initial.root",
    ]

    for fn in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(csv.iterdir())) == 0

    runSeeding(trk_geo, field, outputDir=str(tmp_path), s=seq).run()

    for fn in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 0

    assert_csv_output(csv, "evgen_particles")
    assert_csv_output(csv, "evgen_particles")
    assert_csv_output(csv, "fatras_particles_final")
    assert_csv_output(csv, "fatras_particles_initial")


@pytest.mark.slow
def test_pythia8(tmp_path, seq):
    from pythia8 import runPythia8

    (tmp_path / "csv").mkdir()

    assert not (tmp_path / "pythia8_particles.root").exists()
    assert len(list((tmp_path / "csv").iterdir())) == 0

    runPythia8(str(tmp_path), s=seq).run()

    assert (tmp_path / "pythia8_particles.root").exists()
    assert (tmp_path / "pythia8_particles.root").stat().st_size > 0

    assert len(list((tmp_path / "csv").iterdir())) > 0
    assert_csv_output(tmp_path / "csv", "particles")


def test_propagation(tmp_path, trk_geo, field, seq):
    from propagation import runPropagation

    obj = tmp_path / "obj"
    obj.mkdir()

    root_files = ["propagation_steps.root"]

    for fn in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(obj.iterdir())) == 0

    runPropagation(trk_geo, field, str(tmp_path), s=seq).run()

    for fn in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 0

    assert len(list(obj.iterdir())) > 0


@pytest.mark.slow
def test_geantino_recording(tmp_path, seq):
    from geantino_recording import runGeantinoRecording

    root_files = ["geant-material-tracks.root"]

    dd4hepSvc = acts.examples.dd4hep.DD4hepGeometryService(
        xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    )
    dd4hepG4Construction = acts.examples.geant4.dd4hep.DD4hepDetectorConstruction(
        dd4hepSvc
    )

    for fn in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    runGeantinoRecording(dd4hepG4Construction, str(tmp_path)).run()

    for fn in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 0
