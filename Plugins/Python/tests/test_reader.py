from typing import Type
import inspect

import pytest

import acts
from acts.examples import (
    RootParticleWriter,
    RootParticleReader,
    RootMaterialTrackReader,
    RootTrajectorySummaryReader,
    CsvParticleWriter,
    CsvParticleReader,
    BareAlgorithm,
    Sequencer,
)


class AssertCollectionExistsAlg(BareAlgorithm):
    events_seen = 0

    def __init__(self, collection, *args, **kwargs):
        self.collection = collection
        BareAlgorithm.__init__(self, *args, **kwargs)

    def execute(self, ctx):
        assert ctx.eventStore.exists(self.collection)
        self.events_seen += 1
        return acts.examples.ProcessCode.SUCCESS


def test_root_particle_reader(tmp_path, conf_const, ptcl_gun):
    # need to write out some particles first
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.ERROR)
    evGen = ptcl_gun(s)

    file = tmp_path / "particles.root"
    s.addWriter(
        conf_const(
            RootParticleWriter,
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
        conf_const(
            RootParticleReader,
            acts.logging.ERROR,
            particleCollection="input_particles",
            filePath=str(file),
        )
    )

    alg = AssertCollectionExistsAlg("input_particles", "check_alg", acts.logging.ERROR)
    s2.addAlgorithm(alg)

    s2.run()

    assert alg.events_seen == 10


def test_csv_particle_reader(tmp_path, conf_const, ptcl_gun):
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.ERROR)
    evGen = ptcl_gun(s)

    out = tmp_path / "csv"

    out.mkdir()

    s.addWriter(
        conf_const(
            CsvParticleWriter,
            acts.logging.ERROR,
            inputParticles=evGen.config.outputParticles,
            outputStem="particle",
            outputDir=str(out),
        )
    )

    s.run()

    # reset the seeder
    s = Sequencer(numThreads=1, logLevel=acts.logging.ERROR)

    s.addReader(
        conf_const(
            CsvParticleReader,
            acts.logging.ERROR,
            inputDir=str(out),
            inputStem="particle",
            outputParticles="input_particles",
        )
    )

    alg = AssertCollectionExistsAlg("input_particles", "check_alg", acts.logging.ERROR)

    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


@pytest.mark.parametrize(
    "reader",
    [RootParticleReader, RootTrajectorySummaryReader],
)
@pytest.mark.root
def test_root_reader_interface(reader, conf_const, tmp_path):
    assert hasattr(reader, "Config")

    config = reader.Config

    assert hasattr(config, "filePath")

    kw = {"level": acts.logging.INFO, "filePath": str(tmp_path / "file.root")}

    assert conf_const(reader, **kw)
