#!/usr/bin/env python3

from pathlib import Path
from typing import Annotated, Callable
from colliderml.util import HepMC3Meta
import typer

import ROOT
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mplhep
import numpy as np
import pandas as pd
import pyhepmc
from matplotlib.backends.backend_pdf import PdfPages
import hist
from hist import Hist
import numpy
import colliderml.logging
from rich.progress import (
    Progress,
    MofNCompleteColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from colliderml.root import TH1


def get_colors():
    """
    Returns a list of colors for plotting.
    """
    colors = [
        "C0",  # blue
        "C1",  # orange
        "C2",  # green
        "C3",  # red
        "C4",  # purple
        "C5",  # brown
        "C6",  # pink
        "C7",  # gray
        "C8",  # olive
        "C9",  # cyan
    ]
    return colors


def get_color(i):
    """
    Returns a color for plotting based on the index.
    """
    colors = get_colors()
    if i < 0 or i >= len(colors):
        raise ValueError(f"Index {i} is out of range for colors list.")
    return colors[i]


def get_markers():
    """
    Returns a list of markers for plotting.
    """
    markers = [
        "o",  # circle
        "s",  # square
        "^",  # triangle_up
        "v",  # triangle_down
        "D",  # diamond
        "*",  # star
        "X",  # x
        "+",  # plus
        "x",  # x
        "|",  # vertical line
        "_",  # horizontal line
    ]
    return markers


def get_marker(i):
    """
    Returns a marker for plotting based on the index.
    """
    markers = get_markers()
    if i < 0 or i >= len(markers):
        raise ValueError(f"Index {i} is out of range for markers list.")
    return markers[i]


app = typer.Typer()


@app.command()
def hepmc(
    input_file: Annotated[
        Path,
        typer.Argument(
            ...,
            exists=True,
            dir_okay=False,
            readable=True,
            help="Input HepMC3 file",
        ),
    ],
    output: Annotated[
        Path | None,
        typer.Option(
            ...,
            "-o",
            "--output",
            help="Output plot file ",
        ),
    ] = None,
    max_events: Annotated[
        int | None,
        typer.Option(
            ...,
            "-n",
            "--max-events",
            help="Maximum number of events to process",
        ),
    ] = None,
    show: Annotated[
        bool,
        typer.Option(
            ...,
            help="Show plot interactively (in addition to saving)",
        ),
    ] = False,
    peek_events: int = 100,
):
    """
    Plot HepMC3 file content: vertices, kinematics, and PDG distributions.

    This script reads HepMC3 files and creates comprehensive plots of the event data.
    """

    logger = colliderml.logging.get_logger(__name__)
    # Determine output filename
    if output is None:
        output = input_file.with_suffix(".pdf")

    if output.suffix != ".pdf":
        logger.warning(
            "Output file %s does not end with .pdf, changing to .pdf", output
        )

    logger.info("Input file: %s", input_file)
    logger.info("Output file: %s", output)

    um = 1e-3
    mm = 1

    statuses = set()
    pids = set()

    minima = {}
    maxima = {}

    with pyhepmc.open(input_file) as f:
        for _, event in zip(range(peek_events), f):
            is_beam = event.numpy.particles.status == 4
            statuses.update(event.numpy.particles.status[~is_beam].tolist())
            pids.update(event.numpy.particles.pid[~is_beam].tolist())

            for q in ["px", "py", "pz"]:
                minima.setdefault(q, float("inf"))
                minima[q] = min(
                    minima[q], getattr(event.numpy.particles, q)[~is_beam].min()
                )
                maxima.setdefault(q, float("-inf"))
                maxima[q] = max(
                    maxima[q], getattr(event.numpy.particles, q)[~is_beam].max()
                )

    print("Statuses:", sorted(statuses))
    print("PIDs:", sorted(pids))

    nbins = 100

    hists = dict(
        hist_particles_pid=hist.Hist(
            hist.axis.IntCategory(pids, name="pid", label="Particle ID (PDG)"),
        )
    )

    hists.update(
        {
            f"hist_particles_p{q}": hist.Hist(
                hist.axis.Regular(
                    nbins, -150, 150, name=f"p{q}", label=f"Particle $p_{{{q}}}$ [GeV]"
                ),
            )
            for q in "xy"
        }
    )

    hists["hist_particles_pz"] = hist.Hist(
        hist.axis.Regular(nbins, -600, 600, name="z", label="Particle $p_{z}$ [GeV]"),
    )

    hists["hist_particles_p"] = hist.Hist(
        hist.axis.Regular(nbins, 0, 1000, name="p", label="Particle $p$ [GeV]"),
    )

    hists["hist_particles_pt"] = hist.Hist(
        hist.axis.Regular(
            nbins, 0, 200, name="p", label="Particle $p_\\text{T}$ [GeV]"
        ),
    )

    hists["hist_particles_eta"] = hist.Hist(
        hist.axis.Regular(nbins, -6, 6, name="eta", label="Particle $\\eta$"),
    )

    hists.update(
        {
            f"hist_vertices_{q}": hist.Hist(
                hist.axis.Regular(
                    nbins, -50 * um, 50 * um, name=f"{q}", label=f"Vertex ${q}$ [mm]"
                ),
            )
            for q in "xy"
        }
    )

    hists["hist_vertices_z"] = hist.Hist(
        hist.axis.Regular(
            nbins, -200 * mm, 200 * mm, name="z", label="Vertex $z$ [mm]"
        ),
    )

    meta = HepMC3Meta.for_file(input_file)

    events = meta.num_events
    if max_events is not None:
        events = min(events, max_events)

    prog = Progress(
        *Progress.get_default_columns(),
        MofNCompleteColumn(),
    )
    with prog, pyhepmc.open(input_file) as f:
        print("File open go")
        for _, event in prog.track(
            zip(range(events), f),
            total=events,
            description="Processing events",
        ):
            status = event.numpy.particles.status

            is_beam = status == 4

            particles = {
                q: getattr(event.numpy.particles, q)[~is_beam]
                for q in [
                    "px",
                    "py",
                    "pz",
                ]
            }

            hists["hist_particles_pid"].fill(event.numpy.particles.pid[~is_beam])
            for q in "xyz":
                hists[f"hist_particles_p{q}"].fill(particles[f"p{q}"])
                hists[f"hist_vertices_{q}"].fill(getattr(event.numpy.vertices, q))

            px, py, pz = [particles[f"p{q}"] for q in "xyz"]
            particles_p = numpy.sqrt(px**2 + py**2 + pz**2)
            hists["hist_particles_p"].fill(particles_p)
            particles_pt = numpy.sqrt(px**2 + py**2)
            hists["hist_particles_pt"].fill(particles_pt)

            with numpy.errstate(divide="ignore"):
                particles_eta = numpy.arctanh(
                    pz[particles_p != 0] / particles_p[particles_p != 0]
                )

            hists["hist_particles_eta"].fill(particles_eta)

    for key, histo in hists.items():
        if key != "hist_particles_pid":
            histo = histo[hist.rebin(10)]
        print(key)
        print(histo)

    with PdfPages(output) as pdf:
        for key, histo in hists.items():
            fig, ax = plt.subplots(figsize=(8, 6))

            ax.set_title(key.replace("_", " ").title())
            histo.plot(ax=ax)

            # ax.set_yscale("log")

            pdf.savefig(fig)
            plt.close(fig)

    logger.info("Wrote output to %s", output)


reco_app = typer.Typer()


def common_label(**kwargs):
    kwargs.setdefault("loc", 1)
    kwargs.setdefault("italic", (True, False, False))
    kwargs.setdefault("exp", "ColliderML")
    mplhep.label.exp_text(**kwargs)


def enlarge_top(ax, factor: float = 1.2):
    """
    Enlarge the top of the y-axis by a given factor.
    """
    ymin, ymax = ax.get_ylim()
    delta = ymax - ymin
    ymax = ymin + delta * factor
    ax.set_ylim(ymin, ymax)


@reco_app.command()
def comparison_finding(
    finding_perf: list[Path],
    labels: Annotated[list[str], typer.Option("-l", "--label")],
    output: Path | None = None,
    show: bool = False,
    title: str | None = None,
):

    if len(finding_perf) != len(labels):
        raise ValueError("Number of finding_perf files must match number of labels")

    finding_perf = [ROOT.TFile.Open(str(p.absolute())) for p in finding_perf]

    mplhep.style.use("ATLAS")

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ax.set_xlabel(r"true $\eta$")
    ax.set_ylabel("Technical efficiency")

    ax.set_xlim(-3, 3)

    ax.hlines(1, -3, 3, linestyles="--", color="gray")

    for i, (label, perf) in enumerate(zip(labels, finding_perf)):
        eff_vs_eta = TH1(perf.Get("trackeff_vs_eta"), xrange=(-3, 3))
        eff_vs_eta.errorbar(
            ax,
            label=label,
            marker=get_marker(i),
            linestyle="",
            color=get_color(i),
        )

    ax.legend(loc="upper right")

    enlarge_top(ax=ax, factor=1.15)
    common_label(ax=ax, text="Simulation", supp=title)

    fig.tight_layout()

    if output is not None:
        fig.savefig(output)

    if output is None or show:
        plt.show()


@reco_app.command()
def comparison_fitting(
    fitting_perf: list[Path],
    labels: Annotated[list[str], typer.Option("-l", "--label")],
    output_base: Annotated[Path, typer.Option()],
    show: bool = False,
    title: str | None = None,
):

    if len(fitting_perf) != len(labels):
        raise ValueError("Number of fitting_perf files must match number of labels")

    fitting_perf = [ROOT.TFile.Open(str(p.absolute())) for p in fitting_perf]

    mplhep.style.use("ATLAS")

    params = ["d0", "z0", "phi", "theta", "qop", "t"]
    ylabels = ["d_0", "z_0", r"\phi", r"\theta", "q/p", "t"]
    units = ["[mm]", "[mm]", "", "", "[1/GeV]", "[ns]"]

    for metric in ["pullmean", "pullwidth"]:

        fig, axs = plt.subplots(6, 1, figsize=(8, 8), sharex=True)

        for i, (ax, param, ylabel) in enumerate(zip(axs, params, ylabels)):
            if metric == "pullmean":
                ax.hlines(0, -3, 3, linestyles="--", color="gray")
            elif metric == "pullwidth":
                ax.hlines(1, -3, 3, linestyles="--", color="gray")

            ax.set_xlim(-3, 3)

            if i == 5:
                ax.set_xlabel(r"true $\eta$")
            ax.set_ylabel(f"${ylabel}$")

            for j, (label, perf) in enumerate(zip(labels, fitting_perf)):
                pull_mean = TH1(perf.Get(f"{metric}_{param}_vs_eta"), xrange=(-3, 3))

                pull_mean.errorbar(
                    ax,
                    label=label,
                    marker=get_marker(j),
                    linestyle="",
                    color=get_color(j),
                )

            if metric == "pullmean":
                low, high = ax.get_ylim()
                bound = max(abs(low), abs(high))
                ax.set_ylim(-bound, bound)

            if i == 0:
                ax.legend(bbox_to_anchor=(1.0, 2.0), loc="upper right", ncol=2)
                common_label(ax=ax, text="Simulation", supp=title, loc=0)

            else:
                atlasify.atlasify(axes=ax, outside=True, atlas=False, offset=0)
                ax.get_legend().remove()

        fig.tight_layout(h_pad=-0.01)

        output_file = (
            output_base.parent / f"{output_base.stem}_{metric}{output_base.suffix}"
        )
        print("Saving to", output_file)

        fig.savefig(output_file)

        if show:
            plt.show()

    for param, label, unit in zip(params, ylabels, units):

        fig, ax = plt.subplots(1, 1, figsize=(8, 4))

        ax.set_xlabel(r"true $\eta$")
        ax.set_ylabel(rf"$\sigma_{{{label}}}$ {unit}")

        ax.set_xlim(-3, 3)

        for i, (label, perf) in enumerate(zip(labels, fitting_perf)):
            res_vs_eta = TH1(perf.Get(f"reswidth_{param}_vs_eta"), xrange=(-3, 3))
            res_vs_eta.errorbar(
                ax,
                label=label,
                marker=get_marker(i),
                linestyle="",
                color=get_color(i),
            )

        common_label(ax=ax, text="Simulation", supp=title)
        enlarge_top(ax=ax, factor=1.15)
        ax.legend()

        fig.tight_layout()

        output_file = (
            output_base.parent / f"{output_base.stem}_res_{param}{output_base.suffix}"
        )

        fig.savefig(output_file)


@reco_app.command()
def diagnostics(
    fitting_perf: list[Path],
    labels: Annotated[list[str], typer.Option("-l", "--label")],
    output_base: Annotated[Path, typer.Option()],
    show: bool = False,
    title: str | None = None,
):
    if len(fitting_perf) != len(labels):
        raise ValueError("Number of fitting_perf files must match number of labels")

    fitting_perf = [ROOT.TFile.Open(str(p.absolute())) for p in fitting_perf]

    mplhep.style.use("ATLAS")

    latex = {
        "eta": r"$\eta$",
        "pT": r"$p_T$",
        "phi": r"$\phi$",
        "prodR": "$R$",
        "d0": "$d_0$",
        "z0": "$z_0$",
        "theta": r"$\theta$",
        "qop": "$q/p$",
        "t": "$t$",
    }

    units = {
        "pT": "GeV",
        "prodR": "mm",
        "d0": "mm",
        "z0": "mm",
        "qop": "1/GeV",
        "t": "ns",
    }
    xranges = {"eta": (-3, 3), "phi": (-numpy.pi, numpy.pi)}

    keys = {"trackeff": ["eta", "pT", "phi", "prodR"], "nStates": ["eta", "pT"]}
    for k in ["nHoles", "nOutliers", "nSharedHits", "nMeasurements"]:
        keys[k] = keys["nStates"]

    ylabels = {
        "trackeff": "Technical efficiency",
        "nStates": "Number of states",
        "nHoles": "Number of holes",
        "nOutliers": "Number of outliers",
        "nSharedHits": "Number of shared hits",
        "nMeasurements": "Number of measurements",
    }

    categories = {
        m: "finding"
        for m in [
            "trackeff",
            "nStates",
            "nHoles",
            "nOutliers",
            "nSharedHits",
            "nMeasurements",
        ]
    }

    for param in ["d0", "z0", "phi", "theta", "qop", "t"]:
        keys[f"pullmean_{param}"] = keys["nStates"]
        ylabel = f"Pull mean of {latex[param]}"
        if unit := units.get(param):
            ylabel = f"{ylabel} [{unit}]"
        ylabels[f"pullmean_{param}"] = ylabel
        categories[f"pullmean_{param}"] = "pulls/mean"

        keys[f"pullwidth_{param}"] = keys["nStates"]
        ylabel = f"Pull width of {latex[param]}"
        if unit := units.get(param):
            ylabel = f"{ylabel} [{unit}]"
        ylabels[f"pullwidth_{param}"] = ylabel
        categories[f"pullwidth_{param}"] = "pulls/width"

    for key, quantities in keys.items():
        ylabel = ylabels[key]

        for qty in quantities:

            xlabel = latex[qty]
            if unit := units.get(qty):
                xlabel = f"{xlabel} [{unit}]"

            xrange = xranges.get(qty)

            fig, ax = plt.subplots(1, 1, figsize=(8, 4))

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

            if xrange is not None:
                ax.set_xlim(*xrange)

            for i, (label, perf) in enumerate(zip(labels, fitting_perf)):
                res_vs_eta = TH1(perf.Get(f"{key}_vs_{qty}"), xrange=xrange)
                res_vs_eta.errorbar(
                    ax,
                    label=label,
                    marker=get_marker(i),
                    linestyle="",
                    color=get_color(i),
                )

            if key.startswith("pull"):
                ax.axhline(1 if "width" in key else 0, ls="--", color="gray")

            common_label(ax=ax, text="Simulation", supp=title)
            enlarge_top(ax=ax, factor=1.2)

            if len(fitting_perf) > 1:
                # don't need a legend for single sample
                ax.legend()

            fig.tight_layout()

            output_dir = output_base.parent / categories[key]
            output_dir.mkdir(parents=True, exist_ok=True)

            output_file = (
                output_dir / f"{output_base.stem}_{key}_vs_{qty}{output_base.suffix}"
            )

            print("Saving to", output_file)
            fig.savefig(output_file)

            if show:
                plt.show()
            plt.close(fig)


app.add_typer(reco_app, name="reco")
