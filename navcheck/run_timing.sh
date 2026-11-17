#!/usr/bin/env bash
# Usage: run_timing.sh <binary A> <binary B> <rounds>
S=$(cd "$(dirname "$0")/.." && pwd)
G=$S/geo/itk-geometry.json
T="python3 -I $S/harness/timing.py $1 $2 $G $3"
for c in prompt backward displaced stress; do $T category=$c n=3000 batches=14; done
$T category=prompt n=3000 batches=14 cov=1
$T category=prompt n=3000 batches=14 stepper=sympy cov=1
