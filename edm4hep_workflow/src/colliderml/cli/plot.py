#!/usr/bin/env python3

from pathlib import Path
from typing import Annotated, Callable
from colliderml.util import HepMC3Meta
import typer

import ROOT
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mplhep as hep
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
