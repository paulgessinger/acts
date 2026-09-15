#!/usr/bin/env python3
"""Verify that the core benchmark preserves every exported covariance element."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile


def verify(binary, tracks):
    metadata = tracks.with_name(tracks.name.replace("-tracks.csv", "-metadata.txt"))
    with tempfile.TemporaryDirectory() as directory:
        snapshot = Path(directory) / "parsed.txt"
        process = subprocess.run(
            [
                str(binary),
                "--input",
                str(tracks),
                "--metadata",
                str(metadata),
                "--seeder-only",
                "--warmup",
                "0",
                "--repetitions",
                "1",
                "--output-input",
                str(snapshot),
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        result = json.loads(process.stdout)
        parsed = [
            [float(v) for v in line.split()]
            for line in snapshot.read_text().splitlines()
        ]
    with tracks.open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len(parsed) == result["tracks"]
    names = ["d0", "z0", "phi", "theta", "qop"]
    correlations = 0
    for row, values in zip(rows, parsed):
        expected = [float(row[name]) for name in names] + [0.0]
        for i in range(6):
            for j in range(6):
                if i == 5 or j == 5:
                    value = float(i == j)
                else:
                    key = "var_" + names[i] if i == j else "cov_" + names[i] + names[j]
                    value = float(row[key])
                    correlations += i != j and value != 0.0
                expected.append(value)
        assert (
            values == expected
        ), f'Parameter/covariance mismatch for track {row["trackId"]}'
    return dict(
        file=str(tracks),
        tracks=len(rows),
        nonzero_off_diagonal=correlations,
        seed_z_width=result["seed_z_width"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    files = sorted(args.directory.glob("event*-tracks.csv"))
    if not files:
        parser.error("No exported events found")
    for tracks in files:
        print(json.dumps(verify(args.binary, tracks)))


if __name__ == "__main__":
    main()
