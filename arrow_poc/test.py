#!/usr/bin/env python3
"""
Arrow C Data Interface zero-copy proof-of-concept.

Build the module first (from arrow_poc/):
    cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Then run (pointing PYTHONPATH at the build dir where the .so lives):
    PYTHONPATH=build python test.py
"""

import pyarrow as pa
import polars as pl
import arrow_bridge  # the pybind11 extension built by CMake


def main() -> None:
    print(f"pyarrow version: {pa.__version__}\n")

    # C++ builds the RecordBatch and hands us raw C Data Interface pointers.
    # array_ptr, schema_ptr = arrow_bridge.make_record_batch()
    # print(f"Raw C pointers from C++:")
    # print(f"  ArrowArray*  = {array_ptr:#018x}")
    # print(f"  ArrowSchema* = {schema_ptr:#018x}\n")

    batch = arrow_bridge.make_record_batch(num_tracks=1000000)

    print(f"Schema:\n  {batch.schema}\n")
    print(f"RecordBatch ({batch.num_rows} rows × {batch.num_columns} columns):")

    print("Pandas:")
    df = batch.to_pandas()
    print(df)
    print(df[df["pt"] > 1])

    print("Polars:")
    df = pl.from_arrow(batch, rechunk=False)
    print(df)


if __name__ == "__main__":
    main()
