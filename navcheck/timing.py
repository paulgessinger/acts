"""Interleaved A/B timing of two navcheck binaries.

Usage: timing.py <binary A> <binary B> <geometry> <rounds> [key=value ...]

Each round runs both binaries once, in random order. A run loads the geometry,
does one warm-up batch and then prints the time per track of each batch. The
comparison is the ratio B/A of per-round medians, with a bootstrap interval.
"""

import random, statistics, subprocess, sys

a, b, geometry, rounds = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
extra = [x for x in sys.argv[5:] if not x.startswith("B:")]
extra_b = [x[2:] for x in sys.argv[5:] if x.startswith("B:")]
random.seed(1234)


def run(binary, more=()):
    out = subprocess.run(
        [binary, f"geometry={geometry}", "mode=bench", *extra, *more],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    rows = [l.split() for l in out.splitlines() if l.strip()]
    times = [float(r[0]) for r in rows]
    # The fastest batch is the least disturbed one
    return min(times), int(rows[0][1])


ratios, ta, tb, steps = [], [], [], {}
for _ in range(rounds):
    order = [("a", a), ("b", b)]
    random.shuffle(order)
    res = {}
    for name, binary in order:
        res[name], steps[name] = run(binary, extra_b if name == "b" else ())
    ta.append(res["a"])
    tb.append(res["b"])
    ratios.append(res["b"] / res["a"])

boot = []
for _ in range(4000):
    sample = random.choices(ratios, k=len(ratios))
    boot.append(statistics.median(sample))
boot.sort()
lo, hi = boot[int(0.025 * len(boot))], boot[int(0.975 * len(boot))]
med = statistics.median(ratios)
print(
    f"{' '.join(extra):52s} A={statistics.median(ta):7.3f}us B={statistics.median(tb):7.3f}us "
    f"steps A={steps['a']} B={steps['b']} "
    f"median ratio {100 * (med - 1):+.2f}% [{100 * (lo - 1):+.2f}, {100 * (hi - 1):+.2f}] "
    f"best-run ratio {100 * (min(tb) / min(ta) - 1):+.2f}%"
)
