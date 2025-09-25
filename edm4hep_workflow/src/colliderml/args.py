import argparse
from pathlib import Path


def make_parser():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--jobs", "-j", type=int, default=-1)
    return p
