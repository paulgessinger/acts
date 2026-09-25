"""r-z views of detray portal masks for one ITk volume, before/after the
converter fix. Volume extents come from the ACTS JSON export."""

import json
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

INSIDE = "#2a78d6"
OUTSIDE = "#e34948"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#c9c8c2"
FOCUS_FILL = "#e8eef8"
SURFACE = "#fcfcfb"

acts = json.load(open("/tmp/itk_gen3_noannulus.json"))
before = json.load(open("/tmp/itk_detray/itk_geometry.json"))["data"]["volumes"]
after = json.load(open("/tmp/itk_detray_new/itk_geometry.json"))["data"]["volumes"]


def extent(v):
    b = v["bounds"]["values"]
    z = (v["transform"]["translation"] or [0, 0, 0])[2]
    return b[0], b[1], z - b[2], z + b[2]


leaves = {v["name"]: extent(v) for v in acts["volumes"] if not v.get("children")}


def masks(vols, name):
    names = [v["name"] for v in vols]
    v = vols[names.index(name)]
    out = []
    for s in v["surfaces"]:
        if s["type"] != 0:
            continue
        tz = s["transform"]["translation"][2]
        for m in s["masks"]:
            b = m["boundaries"]
            tgt = names[m["volume_link"]] if m["volume_link"] < len(names) else "end of world"
            if m["shape"] == 4:  # cylinder portal: r, lower z, upper z
                out.append(("cyl", b[0], tz + b[1], tz + b[2], tgt))
            elif m["shape"] == 6:  # ring: inner r, outer r at z
                out.append(("disc", tz, b[0], b[1], tgt))
    return out


def draw(ax, vols, focus, window, title, labels):
    rmin, rmax, zmin, zmax = leaves[focus]
    zlo, zhi, rlo, rhi = window
    for n, (a, b, c, d) in leaves.items():
        if d < zlo or c > zhi or b < rlo or a > rhi or n == focus:
            continue
        ax.add_patch(Rectangle((c, a), d - c, b - a, fc="none", ec=GRID, lw=0.6))
    ax.add_patch(Rectangle((zmin, rmin), zmax - zmin, rmax - rmin, fc=FOCUS_FILL,
                           ec=INK, lw=1.4, zorder=2))
    ms = masks(vols, focus)
    nout = 0
    for kind, p, lo, hi, tgt in ms:
        lim = (zmin, zmax) if kind == "cyl" else (rmin, rmax)
        segs = [(max(lo, lim[0]), min(hi, lim[1]), INSIDE),
                (lo, min(hi, lim[0]), OUTSIDE), (max(lo, lim[1]), hi, OUTSIDE)]
        nout += lo < lim[0] - 1e-3 or hi > lim[1] + 1e-3
        for a, b, col in segs:
            if b <= a:
                continue
            if kind == "cyl":
                ax.plot([a, b], [p, p], color=col, lw=3, solid_capstyle="butt", zorder=3)
            else:
                ax.plot([p, p], [a, b], color=col, lw=3, solid_capstyle="butt", zorder=3)
        # end ticks mark where one mask stops and the next starts
        tick = 0.012 * ((rhi - rlo) if kind == "cyl" else (zhi - zlo))
        for e in (lo, hi):
            if kind == "cyl":
                ax.plot([e, e], [p - tick, p + tick], color=INK, lw=0.8, zorder=4)
            else:
                ax.plot([p - tick, p + tick], [e, e], color=INK, lw=0.8, zorder=4)
    for kind, p, lo, hi, tgt in ms:
        if labels and tgt in labels:
            x, y = ((lo + hi) / 2, p) if kind == "cyl" else (p, (lo + hi) / 2)
            ax.annotate(labels[tgt], (x, y), xytext=(0, 7), textcoords="offset points",
                        ha="center", fontsize=7.5, color=INK2, zorder=5)
    ax.set_xlim(zlo, zhi)
    ax.set_ylim(rlo, rhi)
    ax.set_title(f"{title}: {len(ms)} masks, {nout} leave the volume",
                 fontsize=9.5, color=INK, loc="left")
    ax.set_xlabel("z [mm]", fontsize=8.5, color=INK2)
    ax.set_ylabel("r [mm]", fontsize=8.5, color=INK2)
    ax.tick_params(labelsize=7.5, colors=INK2)
    for sp in ax.spines.values():
        sp.set_color(GRID)


def figure(fname, focus, window, heading, labels=None):
    fig, axes = plt.subplots(2, 1, figsize=(10, 6.4), sharex=True, facecolor=SURFACE)
    for ax in axes:
        ax.set_facecolor(SURFACE)
    draw(axes[0], before, focus, window, "Before (main)", labels)
    draw(axes[1], after, focus, window, "After (#6159)", labels)
    fig.suptitle(heading, fontsize=11, color=INK, x=0.01, ha="left")
    fig.legend(handles=[
        Rectangle((0, 0), 1, 1, fc=FOCUS_FILL, ec=INK, lw=1.4, label=f"{focus}"),
        Rectangle((0, 0), 1, 1, fc="none", ec=GRID, lw=0.6, label="other volumes"),
        Line2D([], [], color=INSIDE, lw=3, label="portal mask, inside the volume"),
        Line2D([], [], color=OUTSIDE, lw=3, label="portal mask, outside the volume"),
    ], loc="lower center", ncol=4, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(f"/tmp/pr_figs/{fname}.svg")
    fig.savefig(f"/tmp/pr_figs/{fname}.png", dpi=110)


figure("case_a_strip_barrel", "Strip_Brl_3", (-3150, 3150, 975, 1025),
       "A: a barrel layer borders a gap volume that spans the whole strip region",
       labels={"Strip::Gap2": "Strip::Gap2", "Strip_Brl::Gap3": "", "Strip_Brl::Gap4": ""})
figure("case_b_outer_pixel_gap", "OuterPixel_nEC_1::Gap17", (-3000, 0, 130, 335),
       "B: a thin endcap gap volume receives every piece of its neighbours' faces")
figure("case_c_staggered_rings", "OuterPixel_nEC_0::Gap3", (-2720, -2230, 140, 275),
       "C: staggered ring stacks, the neighbouring gaps only partly overlap",
       labels={"OuterPixel_nEC_1::Gap2": "nEC_1::Gap2", "OuterPixel_nEC_1::Gap3": "nEC_1::Gap3"})
