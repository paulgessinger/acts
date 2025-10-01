#!/usr/bin/env python3


# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "matplotlib",
#     "mplhep",
#     "numpy",
#     "pyhepmc",
#     "rich",
#     "typer",
# ]
# ///


"""
Plot HepMC3 file content: vertices, particle kinematics, and PDG distributions.

This script reads HepMC3 files and creates comprehensive plots of the event data.
"""

import sys
from pathlib import Path
from typing import Annotated, Optional

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pyhepmc
import typer
from matplotlib.backends.backend_pdf import PdfPages


def analyze_hepmc_file(filename, max_events=None):
    """
    Read HepMC3 file and extract particle and vertex information.

    Returns:
        dict with arrays of vertex positions, particle properties, etc.
    """
    print(f"Reading {filename}...")

    data = {
        "vtx_x": [],
        "vtx_y": [],
        "vtx_z": [],
        "vtx_t": [],
        "particle_x": [],
        "particle_y": [],
        "particle_z": [],
        "particle_px": [],
        "particle_py": [],
        "particle_pz": [],
        "particle_e": [],
        "particle_pdg": [],
    }

    with pyhepmc.open(filename) as f:
        for event_num, event in enumerate(f):
            if max_events is not None and event_num >= max_events:
                break

            # Extract vertex positions
            for vtx in event.vertices:
                pos = vtx.position
                data["vtx_x"].append(pos.x)
                data["vtx_y"].append(pos.y)
                data["vtx_z"].append(pos.z)
                data["vtx_t"].append(pos.t)

            # Extract particle information
            for particle in event.particles:
                mom = particle.momentum
                data["particle_px"].append(mom.px)
                data["particle_py"].append(mom.py)
                data["particle_pz"].append(mom.pz)
                data["particle_e"].append(mom.e)
                data["particle_pdg"].append(particle.pid)

                # Get production vertex position
                prod_vtx = particle.production_vertex
                if prod_vtx:
                    pos = prod_vtx.position
                    data["particle_x"].append(pos.x)
                    data["particle_y"].append(pos.y)
                    data["particle_z"].append(pos.z)
                else:
                    data["particle_x"].append(0.0)
                    data["particle_y"].append(0.0)
                    data["particle_z"].append(0.0)

            if (event_num + 1) % 100 == 0:
                print(f"  Processed {event_num + 1} events...")

    print(f"  Total events processed: {event_num + 1}")

    # Convert to numpy arrays
    for key in data:
        data[key] = np.array(data[key])

    return data


def compute_kinematics(data):
    """Compute derived kinematic quantities."""
    px = data["particle_px"]
    py = data["particle_py"]
    pz = data["particle_pz"]
    e = data["particle_e"]

    # Transverse momentum
    pt = np.sqrt(px**2 + py**2)

    # Total momentum
    p = np.sqrt(px**2 + py**2 + pz**2)

    # Pseudorapidity (eta)
    theta = np.arctan2(pt, pz)
    eta = -np.log(np.tan(theta / 2.0))
    eta = np.where(np.isfinite(eta), eta, 0.0)  # Handle infinities

    # Azimuthal angle (phi)
    phi = np.arctan2(py, px)

    return {
        "pt": pt,
        "p": p,
        "eta": eta,
        "phi": phi,
    }


def plot_analysis(data, kinematics, output_file):
    """Create comprehensive multi-page PDF with individual plots."""
    print(f"Creating plots...")

    # Set ATLAS style
    plt.style.use(hep.style.ATLAS)

    pdg_codes, pdg_counts = np.unique(data["particle_pdg"], return_counts=True)

    # Statistics for title pages
    stats = {
        "total_particles": len(data["particle_pdg"]),
        "total_vertices": len(data["vtx_x"]),
        "unique_pdg": len(pdg_codes),
    }

    with PdfPages(output_file) as pdf:
        # Page 1: Vertex Positions (3 plots)
        fig1, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8.5, 11))

        hep.histplot(
            *np.histogram(data["vtx_x"], bins=100),
            ax=ax1,
            label="Vertex X",
        )
        ax1.set_xlabel("Vertex X [mm]")
        ax1.set_ylabel("Events")
        ax1.set_title("Vertex X Distribution")
        ax1.legend()

        hep.histplot(
            *np.histogram(data["vtx_y"], bins=100),
            ax=ax2,
            label="Vertex Y",
        )
        ax2.set_xlabel("Vertex Y [mm]")
        ax2.set_ylabel("Events")
        ax2.set_title("Vertex Y Distribution")
        ax2.legend()

        hep.histplot(
            *np.histogram(data["vtx_z"], bins=100),
            ax=ax3,
            label="Vertex Z",
        )
        ax3.set_xlabel("Vertex Z [mm]")
        ax3.set_ylabel("Events")
        ax3.set_title("Vertex Z Distribution")
        ax3.legend()

        fig1.suptitle("Vertex Position Distributions", fontsize=14, fontweight="bold")
        fig1.tight_layout()
        pdf.savefig(fig1)
        plt.close(fig1)

        # Page 2: Vertex Time
        fig2, ax = plt.subplots(1, 1, figsize=(8.5, 11))

        hep.histplot(
            *np.histogram(data["vtx_t"], bins=100),
            ax=ax,
            label="Vertex Time",
        )
        ax.set_xlabel("Vertex T [mm/c]")
        ax.set_ylabel("Events")
        ax.set_title("Vertex Time Distribution")
        ax.legend()

        fig2.tight_layout()
        pdf.savefig(fig2)
        plt.close(fig2)

        # Page 3: Particle Production Vertices
        fig3, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8.5, 11))

        hep.histplot(
            *np.histogram(data["particle_x"], bins=100),
            ax=ax1,
            label="Particle Source X",
        )
        ax1.set_xlabel("Particle Source X [mm]")
        ax1.set_ylabel("Particles")
        ax1.set_title("Particle Production X Distribution")
        ax1.legend()

        hep.histplot(
            *np.histogram(data["particle_y"], bins=100),
            ax=ax2,
            label="Particle Source Y",
        )
        ax2.set_xlabel("Particle Source Y [mm]")
        ax2.set_ylabel("Particles")
        ax2.set_title("Particle Production Y Distribution")
        ax2.legend()

        hep.histplot(
            *np.histogram(data["particle_z"], bins=100),
            ax=ax3,
            label="Particle Source Z",
        )
        ax3.set_xlabel("Particle Source Z [mm]")
        ax3.set_ylabel("Particles")
        ax3.set_title("Particle Production Z Distribution")
        ax3.legend()

        fig3.suptitle("Particle Production Vertices", fontsize=14, fontweight="bold")
        fig3.tight_layout()
        pdf.savefig(fig3)
        plt.close(fig3)

        # Page 4: Momentum Distributions
        fig4, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 11))

        hep.histplot(
            *np.histogram(kinematics["p"], bins=100),
            ax=ax1,
            label="Total Momentum",
        )
        ax1.set_xlabel("Total Momentum [GeV/c]")
        ax1.set_ylabel("Particles")
        ax1.set_yscale("log")
        ax1.set_title("Total Particle Momentum")
        ax1.legend()

        hep.histplot(
            *np.histogram(kinematics["pt"], bins=100),
            ax=ax2,
            label="Transverse Momentum",
        )
        ax2.set_xlabel("Transverse Momentum pT [GeV/c]")
        ax2.set_ylabel("Particles")
        ax2.set_yscale("log")
        ax2.set_title("Particle pT Distribution")
        ax2.legend()

        fig4.suptitle("Momentum Distributions", fontsize=14, fontweight="bold")
        fig4.tight_layout()
        pdf.savefig(fig4)
        plt.close(fig4)

        # Page 5: Angular Distributions
        fig5, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 11))

        hep.histplot(
            *np.histogram(kinematics["eta"], bins=100),
            ax=ax1,
            label="Pseudorapidity",
        )
        ax1.set_xlabel("Pseudorapidity η")
        ax1.set_ylabel("Particles")
        ax1.set_title("Eta Distribution")
        ax1.legend()

        hep.histplot(
            *np.histogram(kinematics["phi"], bins=100),
            ax=ax2,
            label="Azimuthal Angle",
        )
        ax2.set_xlabel("Azimuthal Angle φ [rad]")
        ax2.set_ylabel("Particles")
        ax2.set_title("Phi Distribution")
        ax2.legend()

        fig5.suptitle("Angular Distributions", fontsize=14, fontweight="bold")
        fig5.tight_layout()
        pdf.savefig(fig5)
        plt.close(fig5)

        # Page 6: PDG Particle Types
        fig6 = plt.figure(figsize=(8.5, 11))
        ax = fig6.add_subplot(111)

        # Sort by count and take top 30
        sorted_indices = np.argsort(pdg_counts)[::-1][:30]
        top_pdg = pdg_codes[sorted_indices]
        top_counts = pdg_counts[sorted_indices]

        y_pos = np.arange(len(top_pdg))
        ax.barh(y_pos, top_counts, alpha=0.7, edgecolor="black")
        ax.set_yticks(y_pos)
        ax.set_yticklabels([str(int(pdg)) for pdg in top_pdg])
        ax.set_xlabel("Count")
        ax.set_ylabel("PDG Code")
        ax.set_title("Top 30 Particle Types (PDG)", fontsize=14, fontweight="bold")
        ax.invert_yaxis()
        ax.grid(True, alpha=0.3, axis="x")

        # Add statistics
        stats_text = f"Total particles: {stats['total_particles']}\n"
        stats_text += f"Total vertices: {stats['total_vertices']}\n"
        stats_text += f"Unique PDG codes: {stats['unique_pdg']}"

        ax.text(
            0.98,
            0.02,
            stats_text,
            transform=ax.transAxes,
            fontsize=10,
            family="monospace",
            verticalalignment="bottom",
            horizontalalignment="right",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

        fig6.tight_layout()
        pdf.savefig(fig6)
        plt.close(fig6)

        # Set PDF metadata
        d = pdf.infodict()
        d["Title"] = "HepMC3 Event Analysis"
        d["Author"] = "plot_hepmc.py"
        d["Subject"] = "Particle physics event analysis"
        d["Keywords"] = "HepMC3, particle physics, event display"

    print(f"Multi-page PDF saved to {output_file}")
    print(f"  - Page 1: Vertex Position Distributions (X, Y, Z)")
    print(f"  - Page 2: Vertex Time Distribution")
    print(f"  - Page 3: Particle Production Vertices (X, Y, Z)")
    print(f"  - Page 4: Momentum Distributions (p, pT)")
    print(f"  - Page 5: Angular Distributions (η, φ)")
    print(f"  - Page 6: PDG Particle Types")

    return None


def main(
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
):
    """
    Plot HepMC3 file content: vertices, kinematics, and PDG distributions.

    This script reads HepMC3 files and creates comprehensive plots of the event data.
    """
    # Determine output filename
    if output is None:
        output = input_file.with_suffix(".pdf")

    # Analyze file
    data = analyze_hepmc_file(input_file, max_events=max_events)

    # Compute kinematics
    kinematics = compute_kinematics(data)

    # Create plots (multi-page PDF)
    plot_analysis(data, kinematics, output)

    # Note: --show doesn't work with multi-page PDFs
    # User can open the PDF directly to view


if __name__ == "__main__":
    typer.run(main)
