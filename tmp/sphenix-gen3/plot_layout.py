"""r-z and x-y layout plots of the Gen3 sPHENIX geometry.

usage: plot_layout.py <dir>   (reads gen3_* and gen3_zgaps_* CSVs from <dir>)
"""

import csv
import math
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle

D = sys.argv[1]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#d8d7d2"
COL = {  # categorical slots 1-4, fixed order
    "MVTX": "#2a78d6",
    "Silicon": "#eb6834",
    "TPC": "#1baf7a",
    "MICROMEGAS": "#eda100",
}
LABEL = {"MVTX": "MVTX", "Silicon": "INTT", "TPC": "TPC", "MICROMEGAS": "Micromegas"}

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "font.size": 10,
        "axes.titlesize": 11,
    }
)


def read(name):
    return list(csv.DictReader(open(f"{D}/{name}")))


def subsystem_of(volname):
    for s in COL:
        if volname.startswith(s):
            return s
    return None


def load(prefix):
    vols = read(f"{prefix}_volumes.csv")
    sens = read(f"{prefix}_sensors.csv")
    mats = read(f"{prefix}_material_surfaces.csv")
    for s in sens:
        s["v"] = [
            tuple(float(s[f"{c}{i}"]) for c in "xyz") for i in range(4)
        ]
    for m in mats:
        for k in ("cx", "cy", "r", "z0", "z1"):
            m[k] = float(m[k])
    return vols, sens, mats


def legend_handles():
    h = [
        Line2D([], [], color=COL[s], lw=2, label=f"{LABEL[s]} sensors")
        for s in COL
    ]
    h += [
        Line2D([], [], color=INK, lw=1.6, ls=(0, (4, 2)), label="material carrier (passive)"),
        Line2D([], [], color=INK2, lw=0.9, label="material on portal"),
        Rectangle((0, 0), 1, 1, fc="none", ec=GRID, lw=0.8, label="volume boundary"),
    ]
    return h


# ---------------------------------------------------------------- r-z
def draw_rz(ax, vols, sens, mats, rmax, zmax, title):
    ax.set_title(title, loc="left")
    # volumes: light fill for subsystem volumes, outline for all
    for v in vols:
        if v["name"] == "World":
            continue
        rmin, rmx, hz, cz = (float(v[k]) for k in ("rmin", "rmax", "hz", "cz"))
        sub = subsystem_of(v["name"])
        is_gap = "Gap" in v["name"]
        fc = COL[sub] if (sub and not is_gap) else "none"
        ax.add_patch(
            Rectangle(
                (cz - hz, rmin), 2 * hz, rmx - rmin,
                fc=fc, alpha=0.08 if fc != "none" else 1, ec="none",
            )
        )
        ax.add_patch(
            Rectangle((cz - hz, rmin), 2 * hz, rmx - rmin, fc="none", ec=GRID, lw=0.6)
        )
    # sensors: (z, r) outline of each surface
    for sub in COL:
        polys = [
            [(p[2], math.hypot(p[0], p[1])) for p in (s["v"][0], s["v"][1], s["v"][3], s["v"][2])]
            for s in sens
            if s["subsystem"] == sub
        ]
        ax.add_collection(
            PolyCollection(polys, facecolors=COL[sub], edgecolors=COL[sub], lw=0.4, alpha=0.9)
        )
    # material
    for m in mats:
        if m["kind"] == "portal":
            ax.plot([m["z0"], m["z1"]], [m["r"]] * 2, color=INK2, lw=0.9, zorder=3)
        else:
            off = math.hypot(m["cx"], m["cy"])
            # off-axis carrier: r from the beam axis varies with phi -> band
            if off > 0.1:
                ax.fill_between(
                    [m["z0"], m["z1"]], m["r"] - off, m["r"] + off,
                    color=INK, alpha=0.07, lw=0, zorder=2,
                )
            ax.plot(
                [m["z0"], m["z1"]], [m["r"]] * 2,
                color=INK, lw=1.6, ls=(0, (4, 2)), zorder=4,
            )
    # eta guides
    for eta in (0.5, 1.0):
        th = 2 * math.atan(math.exp(-eta))
        for s in (1, -1):
            zz = min(zmax, rmax / math.tan(th))
            ax.plot([0, s * zz], [0, zz * math.tan(th)], ls=":", color=INK2, lw=0.7)
        ax.text(
            min(zmax, rmax / math.tan(th)) * 0.97, min(rmax, zmax * math.tan(th)) * 0.97,
            f"η={eta:g}", color=INK2, fontsize=8, ha="right", va="top",
        )
    ax.set_xlim(-zmax, zmax)
    ax.set_ylim(0, rmax)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("r [mm]")


def label_rz(ax, items):
    for text, z, r, sub in items:
        ax.text(z, r, text, color=INK, fontsize=9, ha="left", va="center",
                bbox=dict(fc=SURFACE, ec=COL[sub], lw=1, pad=2, boxstyle="round,pad=0.25"))


def rz_figure():
    fig, axs = plt.subplots(2, 2, figsize=(16, 11.5),
                            gridspec_kw=dict(width_ratios=[1.15, 1]))
    for row, (prefix, mode) in enumerate(
        [("gen3", "default: R-stack expands every volume to full z"),
         ("gen3_zgaps", "--zgaps: Gen1-like gap | barrel | gap in z")]
    ):
        vols, sens, mats = load(prefix)
        draw_rz(axs[row][0], vols, sens, mats, 900, 1150, f"Full detector — {mode}")
        draw_rz(axs[row][1], vols, sens, mats, 125, 260, "Inner silicon (MVTX + INTT)")
        label_rz(axs[row][0], [("Micromegas", -1120, 880, "MICROMEGAS"),
                               ("TPC (48 layer volumes)", -1000, 540, "TPC"),
                               ("MVTX + INTT", -1000, 70, "Silicon")])
        label_rz(axs[row][1], [("INTT: 1 volume, 4 carriers", -250, 117, "Silicon"),
                               ("MVTX: 1 volume, 3 carriers\n(off-axis: shaded = r range of carrier)",
                                -250, 58, "MVTX")])
    fig.legend(handles=legend_handles(), loc="lower center", ncol=7, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("sPHENIX Gen3 (TGeo BlueprintBuilder) — r–z layout", x=0.01, ha="left",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(f"{D}/sphenix_gen3_rz.png", dpi=120)


# ---------------------------------------------------------------- x-y
def draw_xy(ax, vols, sens, mats, lim, title, center=(0, 0), show_tpc_layers=True):
    ax.set_title(title, loc="left")
    cx0, cy0 = center
    # volume boundaries (coaxial) as circles
    radii = set()
    for v in vols:
        if v["name"] == "World":
            continue
        if not show_tpc_layers and v["name"].startswith("TPC"):
            continue
        radii.update((float(v["rmin"]), float(v["rmax"])))
    for r in radii:
        if r > 0:
            ax.add_patch(Circle((0, 0), r, fc="none", ec=GRID, lw=0.6))
    # sensors: xy footprint = segment between the farthest pair of vertices
    for sub in COL:
        segs = []
        for s in sens:
            if s["subsystem"] != sub:
                continue
            pts = [(p[0], p[1]) for p in s["v"]]
            a, b = max(
                ((p, q) for p in pts for q in pts),
                key=lambda pq: (pq[0][0] - pq[1][0]) ** 2 + (pq[0][1] - pq[1][1]) ** 2,
            )
            segs.append([a, b])
        ax.add_collection(LineCollection(segs, colors=COL[sub], lw=1.6 if sub != "TPC" else 0.5))
    for m in mats:
        if m["kind"] == "portal":
            ax.add_patch(Circle((m["cx"], m["cy"]), m["r"], fc="none", ec=INK2, lw=0.8))
        else:
            ax.add_patch(Circle((m["cx"], m["cy"]), m["r"], fc="none", ec=INK, lw=1.4,
                                ls=(0, (4, 2)), zorder=4))
    ax.plot([0], [0], marker="+", color=INK, ms=10, mew=1.2)
    ax.set_xlim(cx0 - lim, cx0 + lim)
    ax.set_ylim(cy0 - lim, cy0 + lim)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")


def xy_figure():
    vols, sens, mats = load("gen3")
    fig, axs = plt.subplots(1, 3, figsize=(20, 7.6))
    draw_xy(axs[0], vols, sens, mats, 900, "Full detector (z-projection)")
    draw_xy(axs[1], vols, sens, mats, 125, "MVTX + INTT")
    draw_xy(axs[2], vols, sens, mats, 50, "MVTX (axis offset from beam line)")
    # mark the MVTX axis
    for m in mats:
        if m["subsystem"] == "MVTX" and m["kind"] == "carrier":
            axs[2].plot([m["cx"]], [m["cy"]], marker="x", color=COL["MVTX"], ms=9, mew=2)
            axs[2].annotate(
                f"MVTX axis ({m['cx']:.1f}, {m['cy']:.1f}) mm", (m["cx"], m["cy"]),
                xytext=(m["cx"] + 4, m["cy"] - 9), fontsize=9,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8),
            )
            break
    axs[2].annotate("beam line", (0, 0), xytext=(-30, 6), fontsize=9,
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    axs[1].text(0, -118, "INTT: 2 staggered sub-layers per barrel → 4 carriers",
                ha="center", fontsize=9, color=INK)
    fig.legend(handles=legend_handles(), loc="lower center", ncol=7, frameon=False)
    fig.suptitle("sPHENIX Gen3 (TGeo BlueprintBuilder) — x–y layout (all z projected)",
                 x=0.01, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(f"{D}/sphenix_gen3_xy.png", dpi=120)


rz_figure()
xy_figure()
