#!/usr/bin/env python3
import os
from pathlib import Path

from acts.examples import (
    TGeoDetector,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
    Interval,
)

import acts

from acts import MaterialMapJsonConverter, UnitConstants as u


def runITk(
    trackingGeometry,
    decorators,
    outputDir,
    events=1,
    outputObj=True,
    outputCsv=False,
):

    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0

        context = AlgorithmContext(ialg, ievt, eventStore)

        for cdr in decorators:
            r = cdr.decorate(context)
            if r != ProcessCode.SUCCESS:
                raise RuntimeError("Failed to decorate event context")

        if outputCsv:
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=os.path.join(outputDir, "csv"),
                writePerEvent=True,
            )
            writer.write(context)

        if outputObj:
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO, outputDir=os.path.join(outputDir, "obj")
            )
            writer.write(context, trackingGeometry)


if "__main__" == __name__:
    # detector, trackingGeometry, decorators = AlignedDetector.create()
    # detector, trackingGeometry, decorators = GenericDetector.create()
    # detector, trackingGeometry, decorators = getOpenDataDetector()

    # runGeometry(trackingGeometry, decorators, outputDir=os.getcwd())

    geo_example_dir = Path.cwd() / "../acts-detector-example"
    assert geo_example_dir.exists(), "Detector example input directory missing"

    # detector, trackingGeometry, decorators = TGeoDetector.create(
    #     fileName=str(geo_example_dir / "atlas/itk-hgtd/ATLAS-ITk-HGTD.tgeo.root")
    # )

    Volume = TGeoDetector.Config.Volume
    LayerTriplet = TGeoDetector.Config.LayerTriplet

    itkConfig = TGeoDetector.Config(
        fileName=str(geo_example_dir / "atlas/itk-hgtd/ATLAS-ITk-HGTD.tgeo.root"),
        buildBeamPipe=True,
        unitScalor=1.0,  # explicit units
        beamPipeRadius=29.0 * u.mm,
        beamPipeHalflengthZ=3000.0 * u.mm,
        beamPipeLayerThickness=0.8 * u.mm,
        volumes=[
            Volume(
                name="Pixel::Pixel",
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                sensitiveNames=LayerTriplet(["Pixel::siLog"]),
                sensitiveAxes=LayerTriplet("YZX"),
                rRange=LayerTriplet((0 * u.mm, 135 * u.mm)),
                zRange=LayerTriplet(
                    negative=(-3000 * u.mm, -250 * u.mm),
                    central=(-250 * u.mm, 250 * u.mm),
                    positive=(250 * u.mm, 3000 * u.mm),
                ),
                splitTolR=LayerTriplet(negative=-1.0, central=5 * u.mm, positive=-1.0),
                splitTolZ=LayerTriplet(
                    negative=10 * u.mm, central=-1.0, positive=10 * u.mm
                ),
            )
        ],
    )
