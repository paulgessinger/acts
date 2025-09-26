#!/usr/bin/env python3

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "typer",
#   "matplotlib",
#   "numpy",
# ]
# ///

from pathlib import Path
import typer
import re
import matplotlib.pyplot as plt
import numpy as np


def main(infile: Path):
    positions = []

    with infile.open() as fh:
        for line in fh:
            # m = re.match(r".*endpoint: (.*), (.*), (.*)", line)
            # │09:54:43    EDM4hepSimIn   VERBOSE     - at -0.00830637   0.0045211    -56.5942
            m = re.match(r".*at +([-\d.]+) +([-\d.]+) +([-\d.]+) +([-\d.]+)", line)
            if m is None:
                continue
            print(line)
            print(m.groups())
            # assert m is not None, f"Could not parse line: {line}"
            x, y, z, t = map(float, m.groups())
            positions.append((x, y, z, t))

    positions = np.array(positions)

    # positions = positions[:100]
    # print(positions)

    xs, ys, zs, ts = positions.T
    rs = np.sqrt(xs**2 + ys**2)

    fig, (ax_xy, ax_rz) = plt.subplots(1, 2, figsize=(12, 8))

    cmap = plt.cm.viridis

    idxs = np.argsort(ts) / len(ts)

    # print(xs)
    # print(idxs)

    colors = cmap(idxs)

    ax_xy.scatter(xs, ys, c=colors, s=10)
    ax_xy.set(xlabel="x [mm]", ylabel="y [mm]", title="XY plane")

    ax_rz.scatter(rs, zs, c=colors, s=10)
    ax_rz.set(xlabel="z [mm]", ylabel="r [mm]", title="RZ plane")

    plt.show()


typer.run(main)
