#!/usr/bin/env python3

from pathlib import Path
import typer

import edm4hep
import ROOT

# if ROOT.gSystem.Load("libActsPodioEdmDict") != 0:
#     raise RuntimeError("Failed to load libActsPodioEdmDict")
# from ROOT import ActsPodioEdm
import podio
import cppyy


def main(input: Path, category: str = "events"):
    print(input.resolve())
    assert input.exists(), "Input file does not exist"
    reader = podio.root_io.Reader([input])

    assert category in reader.categories, f"Category {category} not in file"

    for frame in reader.get(category):
        print(frame.getAvailableCollections())
        # meas = frame.get("ActsMeasurements")


typer.run(main)
