#!/usr/bin/env python3
import pathlib, acts, acts.examples
import acts.examples.dd4hep
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector

# acts.examples.dump_args_calls(locals())  # show python binding calls

u = acts.UnitConstants
outputDir = pathlib.Path.cwd() / "odd_output"
outputDir.mkdir(exist_ok=True)

oddDir = getOpenDataDetectorDirectory()

oddMaterialMap = oddDir / "data/odd-material-maps.root"
oddDigiConfig = oddDir / "config/odd-digi-smearing-config.json"
oddSeedingSel = oddDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector, trackingGeometry, decorators = getOpenDataDetector(
    getOpenDataDetectorDirectory(), mdecorator=oddMaterialDeco
)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addFatras,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    addCKFTracks,
    CKFPerformanceConfig,
    addVertexFitting,
    VertexFinder,
)

s = acts.examples.Sequencer(events=1, numThreads=-1, logLevel=acts.logging.INFO)

vtxGen = acts.examples.GaussianVertexGenerator(
    stddev=acts.Vector4(10 * u.um, 10 * u.um, 50 * u.mm, 0),
    mean=acts.Vector4(0, 0, 0, 0),
)

addParticleGun(
    s,
    MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
    EtaConfig(-3.0, 3.0, uniform=True),
    ParticleConfig(4, acts.PdgParticle.eMuon, randomizeCharge=True),
    vtxGen=vtxGen,
    multiplicity=50,
    rnd=rnd,
)
addFatras(
    s,
    trackingGeometry,
    field,
    outputDirRoot=outputDir,
    rnd=rnd,
)
addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=oddDigiConfig,
    outputDirRoot=outputDir,
    rnd=rnd,
)
addSeeding(
    s,
    trackingGeometry,
    field,
    geoSelectionConfigFile=oddSeedingSel,
    outputDirRoot=outputDir,
)
addCKFTracks(
    s,
    trackingGeometry,
    field,
    CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
    outputDirRoot=outputDir,
)
s.addAlgorithm(
    acts.examples.TrackSelector(
        level=acts.logging.INFO,
        inputTrackParameters="fittedTrackParameters",
        outputTrackParameters="trackparameters",
        outputTrackIndices="outputTrackIndices",
        removeNeutral=True,
        absEtaMax=2.5,
        loc0Max=4.0 * u.mm,  # rho max
        ptMin=500 * u.MeV,
    )
)
addVertexFitting(
    s,
    field,
    vertexFinder=VertexFinder.AMVF,
    outputDirRoot=outputDir,
    logLevel=acts.logging.DEBUG,
)

s.run()
