#!/usr/bin/env python3
import os

import acts
import acts.examples

import acts.examples.dd4hep
import acts.examples.geant4
import acts.examples.geant4.dd4hep

u = acts.UnitConstants


def runGeantinoRecording(g4_construction, outputDir, s=None):

    g4AlgCfg = acts.examples.geant4.GeantinoRecording.Config()
    g4AlgCfg.detectorConstruction = g4_construction
    g4AlgCfg.tracksPerEvent = 100

    g4Alg = acts.examples.geant4.GeantinoRecording(
        level=acts.logging.INFO, config=g4AlgCfg
    )

    s = s or acts.examples.Sequencer(events=100, numThreads=1)
    s.addAlgorithm(g4Alg)

    s.addWriter(
        acts.examples.RootMaterialTrackWriter(
            prePostStep=True,
            recalculateTotals=True,
            collection=g4Alg.config.outputMaterialTracks,
            filePath=outputDir + "/" + g4Alg.config.outputMaterialTracks + ".root",
            level=acts.logging.INFO,
        )
    )

    return s


if "__main__" == __name__:
    dd4hepSvc = acts.examples.dd4hep.DD4hepGeometryService(
        xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    )
    dd4hepG4Construction = acts.examples.geant4.dd4hep.DD4hepDetectorConstruction(
        dd4hepSvc
    )
    runGeantinoRecording(
        g4_construction=dd4hepG4Construction, outputDir=os.getcwd()
    ).run()
