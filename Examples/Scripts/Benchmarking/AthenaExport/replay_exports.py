#!/usr/bin/env python3
"""Replay exported events against two ACTS builds and an optional reference."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--baseline-binary", type=Path)
    parser.add_argument("--baseline-library", type=Path, required=True)
    parser.add_argument("--optimized-library", type=Path, required=True)
    parser.add_argument("--reference-library", type=Path)
    parser.add_argument("--reference-binary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    if args.reference_binary and not args.reference_library:
        parser.error("--reference-binary requires --reference-library")
    libraries = dict(baseline=args.baseline_library, optimized=args.optimized_library)
    binaries = dict(baseline=args.baseline_binary or args.binary, optimized=args.binary)
    if args.reference_library:
        libraries["reference"] = args.reference_library
        binaries["reference"] = args.reference_binary or args.binary
    files = sorted(args.directory.glob("event*-tracks.csv"))
    if args.limit is not None:
        files = files[: args.limit]
    if not files:
        parser.error("No exported events found")
    records = []
    with tempfile.TemporaryDirectory() as directory:
        for event, tracks in enumerate(files):
            metadata = tracks.with_name(
                tracks.name.replace("-tracks.csv", "-metadata.txt")
            )
            for mode in ["seeder", "amvf"]:
                results = {}
                order = (
                    ["baseline", "optimized"]
                    if event % 2 == 0
                    else ["optimized", "baseline"]
                )
                if args.reference_library:
                    variants = list(libraries)
                    offset = event % len(variants)
                    order = variants[offset:] + variants[:offset]
                    if (event // len(variants)) % 2:
                        order.reverse()
                for variant in order:
                    dump = Path(directory) / f"{variant}.txt"
                    command = [
                        "taskset",
                        "-c",
                        str(args.cpu),
                        str(binaries[variant]),
                        "--input",
                        str(tracks),
                        "--metadata",
                        str(metadata),
                        "--athena-era",
                        "run4",
                        "--warmup",
                        str(args.warmup),
                        "--repetitions",
                        str(args.repetitions),
                    ]
                    command += (
                        ["--seeder-only"]
                        if mode == "seeder"
                        else ["--output-vertices", str(dump)]
                    )
                    library = libraries[variant]
                    environment = dict(
                        os.environ,
                        LD_LIBRARY_PATH=str(library)
                        + ":"
                        + os.environ.get("LD_LIBRARY_PATH", ""),
                    )
                    result = subprocess.run(
                        command,
                        env=environment,
                        text=True,
                        capture_output=True,
                        check=True,
                    )
                    timing = json.loads(result.stdout)
                    results[variant] = timing
                    records.append(
                        dict(event=event, mode=mode, variant=variant, **timing)
                    )
                if mode == "seeder":
                    assert all(
                        result["seed_z_width"] == results["baseline"]["seed_z_width"]
                        for result in results.values()
                    ), tracks
                else:
                    expected = (Path(directory) / "baseline.txt").read_bytes()
                    assert all(
                        (Path(directory) / f"{variant}.txt").read_bytes() == expected
                        for variant in results
                    ), tracks
                print(
                    f"{tracks.name} {mode}: "
                    + "; ".join(
                        f"{variant} {results[variant]['median_ms']:.3f} ms"
                        for variant in libraries
                    )
                    + "; identical output",
                    flush=True,
                )
                args.output.write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
