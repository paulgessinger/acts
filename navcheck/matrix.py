"""Run the navigation check over a matrix of configurations in parallel."""

import itertools, json, subprocess, sys, os, collections
from concurrent.futures import ThreadPoolExecutor

binary, geometry, outdir = sys.argv[1:4]
extra = sys.argv[4:]
os.makedirs(outdir, exist_ok=True)
N = {"prompt": 3000, "backward": 3000, "displaced": 1500, "stress": 1500}
configs = list(
    itertools.product(
        ["eigen", "sympy"],
        ["const", "tilted", "varying", "zero"],
        [0, 1],
        ["prompt", "backward", "displaced", "stress"],
        [0, 1],
    )
)


def hits(path):
    c = collections.Counter()
    for l in open(path):
        c[l] += 1
    return c


def run(cfg):
    stepper, field, kick, category, fix = cfg
    tag = f"{stepper}-{field}-k{kick}-{category}-f{fix}"
    cmd = [
        binary,
        f"geometry={geometry}",
        f"stepper={stepper}",
        f"field={field}",
        f"kick={kick}",
        f"category={category}",
        f"fix={fix}",
        f"n={N[category]}",
        f"detail={outdir}/detail-{tag}.txt",
        f"hits={outdir}/hits-{tag}.txt",
        *extra,
    ]
    out = subprocess.run(cmd, capture_output=True, text=True)
    lines = [l for l in out.stdout.splitlines() if l.startswith("{")]
    if not lines:
        return cfg, {"error": out.returncode, "stderr": out.stderr[-300:]}
    return cfg, json.loads(lines[-1])


with ThreadPoolExecutor(12) as pool:
    results = dict(pool.map(run, configs))
json.dump(
    {"-".join(map(str, k)): v for k, v in results.items()},
    open(f"{outdir}/results.json", "w"),
    indent=1,
)

print(
    f"{'stepper':7s} {'field':8s} k {'category':10s} | {'expected':>8s} | "
    f"{'main: miss late spur':>21s} | {'fix: miss late spur':>20s} | "
    f"{'steps':>9s} {'d steps':>7s} | {'hits +/-':>9s} | fail retr"
)
tot = collections.Counter()
for stepper, field, kick, category in itertools.product(
    ["eigen", "sympy"],
    ["const", "tilted", "varying", "zero"],
    [0, 1],
    ["prompt", "backward", "displaced", "stress"],
):
    a = results[(stepper, field, kick, category, 0)]
    b = results[(stepper, field, kick, category, 1)]
    if "error" in a or "error" in b:
        print(stepper, field, kick, category, "ERROR", a if "error" in a else b)
        continue
    tag = f"{stepper}-{field}-k{kick}-{category}"
    ha, hb = hits(f"{outdir}/hits-{tag}-f0.txt"), hits(f"{outdir}/hits-{tag}-f1.txt")
    gained, lost = sum((hb - ha).values()), sum((ha - hb).values())
    if kick:
        gained = lost = -1  # trajectories differ once a surface is missed
    print(
        f"{stepper:7s} {field:8s} {kick} {category:10s} | {b['expected']:8d} | "
        f"{a['missing']:9d} {a['late']:5d} {a['spurious']:5d} | "
        f"{b['missing']:8d} {b['late']:5d} {b['spurious']:5d} | "
        f"{b['steps']:9d} {b['steps'] - a['steps']:7d} | {gained:4d} {lost:4d} | "
        f"{b['failed']:4d} {b['kickRetraced']:4d}"
    )
    for k in ("expected", "missing", "late", "spurious", "failed", "kickRetraced"):
        tot["fix_" + k] += b[k]
        tot["main_" + k] += a[k]
    if not kick:
        tot["gained"] += gained
        tot["lost"] += lost
print(dict(tot))
