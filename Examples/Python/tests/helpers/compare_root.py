#!/usr/bin/env python3
import sys
import argparse
from pathlib import Path

import uproot
import awkward as ak


def compare_root_files(file_a: Path, file_b: Path) -> bool:

    result = True

    rf_a = uproot.open(file_a)
    rf_b = uproot.open(file_b)

    assert set(rf_a.keys()) == set(rf_b.keys()), "Files have different contents"

    for tree_name in sorted(rf_a.keys()):
        tree_a = rf_a[tree_name]
        tree_b = rf_b[tree_name]

        assert set(tree_a.keys()) == set(
            tree_a.keys()
        ), f"Tree {tree_name} has different keys:\n   A: {list(tree_a.keys())} != B: {list(tree_b.keys())}"

        keys = list(tree_a.keys())

        branches_a = tree_a.arrays(library="ak")
        branches_b = tree_b.arrays(library="ak")

        for row_a, row_b in zip(
            zip(*[branches_a[b] for b in keys]),
            zip(*[branches_b[b] for b in keys]),
        ):

            for key, (obj_a, obj_b) in zip(keys, zip(row_a, row_b)):

                if type(obj_a) != type(obj_b):
                    print("Type mismatch:", key)
                    print(" A =>", type(obj_a))
                    print(" B =>", type(obj_b))
                    result = False
                    continue

                same = True
                #  print(key)
                if isinstance(obj_a, ak.highlevel.Array):
                    ba = obj_a == obj_b
                    #  print("ba", ba.to_list())
                    same = ak.all(ba)
                else:
                    same = obj_a == obj_b

                if not same:
                    print("Value mismatch:", key)
                    print(" A =>", obj_a.to_list())
                    print(" B =>", obj_b.to_list())

                    result = False
                    print("Lengths: A|B:", len(obj_a), len(obj_b))

                    for idx, (v_a, v_b) in enumerate(zip(obj_a, obj_b)):
                        if v_a != v_b:
                            print(" -> mismatch at idx", idx, ": A|B: ", v_a, "!=", v_b)
        if result:
            print("Files are identical")
    return result


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Compare two root files, not this is not invariant under event reordering"
    )

    p.add_argument("file_a", type=Path, help="First file to compare")

    p.add_argument("file_b", type=Path, help="Second file to compare")

    args = p.parse_args()

    res = compare_root_files(
        file_a=args.file_a,
        file_b=args.file_b,
    )

    if not res:
        sys.exit(1)
