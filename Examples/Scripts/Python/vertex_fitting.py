#!/usr/bin/env python3
from pathlib import Path
from typing import Optional

from acts.examples import (
    Sequencer,
    ParticleSelector,
    ParticleSmearing,
    TruthVertexFinder,
    VertexFitterAlgorithm,
    IterativeVertexFinderAlgorithm,
    RootParticleReader,
)

import acts

from acts import UnitConstants as u

from common import addPythia8, getOpenDataDetector


def runVertexFitting(
    trackingGeometry,
    decorators,
    field,
    outputDir: Path,
    inputParticlePath: Optional[Path] = None,
    truthVertexFinder: bool = True,
    s=None,
):
    s = s or Sequencer(events=100, numThreads=-1)

    logger = acts.logging.getLogger("VertexFittingExample")

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    if inputParticlePath is None:
        logger.info("Generating particles using Pythia8")
        evGen = addPythia8(s, rnd)
        inputParticles = evGen.config.outputParticles
    else:
        logger.info("Reading particles from %s", inputParticlePath)
        assert inputParticlePath.exists()
        inputParticles = "particles_read"
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath),
                particleCollection=inputParticles,
            )
        )

    ptclSelector = ParticleSelector(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
        outputParticles="particles_selected",
        removeNeutral=True,
        absEtaMax=2.5,
        rhoMax=4.0,
        ptMin=500 * u.MeV,
    )
    s.addAlgorithm(ptclSelector)

    ptclSmearing = ParticleSmearing(
        level=acts.logging.INFO,
        inputParticles=ptclSelector.config.outputParticles,
        outputTrackParameters="trackparameters",
        randomNumbers=rnd,
    )
    s.addAlgorithm(ptclSmearing)

    if truthVertexFinder:
        logger.info("Using truth vertex finder")
        findVertices = TruthVertexFinder(
            level=acts.logging.INFO,
            inputParticles=ptclSelector.config.outputParticles,
            outputProtoVertices="protovertices",
            excludeSecondaries=True,
        )
    else:
        logger.info("Using iterative vertex finder")
        findVertices = IterativeVertexFinderAlgorithm(
            level=acts.logging.INFO,
            bField=field,
            inputTrackParameters=ptclSmearing.config.outputTrackParameters,
            outputProtoVertices="protovertices",
        )
    s.addAlgorithm(findVertices)

    fitVertices = VertexFitterAlgorithm(
        level=acts.logging.INFO,
        bField=field,
        inputTrackParameters=ptclSmearing.config.outputTrackParameters,
        inputProtoVertices=findVertices.config.outputProtoVertices,
        outputVertices="fittedvertices",
    )
    s.addAlgorithm(fitVertices)

    return s


if "__main__" == __name__:
    detector, trackingGeometry, decorators = getOpenDataDetector()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("pythia8_particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None

    runVertexFitting(
        trackingGeometry,
        decorators,
        field,
        truthVertexFinder=False,
        inputParticlePath=inputParticlePath,
        outputDir=Path.cwd(),
    ).run()
