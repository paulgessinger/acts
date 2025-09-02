#!/usr/bin/env python3

import pyhepmc
from pathlib import Path
import sys

input_file, output_file = sys.argv[1:]

print(input_file, "->", output_file)
with pyhepmc.open(input_file, "r") as f:

    def _count(i):
        for _ in i:
            yield 1

    Path(output_file).write_text(str(sum(_count(f))))
