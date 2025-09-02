#!/usr/bin/env python3

from pathlib import Path
import time

from DDSim.DD4hepSimulation import DD4hepSimulation
from acts.examples.odd import getOpenDataDetectorDirectory
from g4units import GeV

from args import make_parser

p = make_parser()
p.add_argument("--input", type=Path, required=True)
p.add_argument("--minimalKineticEnergy", type=float, default=1.0)
p.add_argument("--seed", type=int, default=int(time.time()))

args = p.parse_args()

args.output.parent.mkdir(parents=True, exist_ok=True)


odd_dir = getOpenDataDetectorDirectory()
odd_xml = odd_dir / "xml" / "OpenDataDetector.xml"

ddsim = DD4hepSimulation()

# Configure DD4hep
if isinstance(ddsim.compactFile, list):
    ddsim.compactFile = [str(odd_xml)]
else:
    ddsim.compactFile = str(odd_xml)

# @TODO: Deal with custom particle handler

# https://github.com/OpenDataDetector/ColliderML/blob/75ad4313a7f2b1bf86ea140393d3fd9c348a0fba/scripts/simulation/ddsim_run.py#L175

count_file = args.input.with_suffix(args.input.suffix + ".count")
if not count_file.exists():
    raise RuntimeError(f"Count file {count_file} does not exist")

n_events = int(count_file.read_text())

print(n_events, "events in", count_file)

ddsim.part.minimalKineticEnergy = args.minimalKineticEnergy * GeV

ddsim.inputFiles = [str(args.input.resolve())]
ddsim.numberOfEvents = n_events

ddsim.outputFile = str(args.output.resolve())

ddsim.random.seed = args.seed

ddsim.run()
