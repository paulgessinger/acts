#!/usr/bin/env python3

from pathlib import Path
from typing import Annotated, Callable
import typer

import ROOT


def split(delim: str) -> Callable[[str], list[str]]:
    def fn(arg: str):
        return arg.split(delim)

    return fn


def plot_from_file(
    file_path: Path,
    plots: list[str],
    output_dir: Path,
    output_file: Path,
    canvas: ROOT.TCanvas,
    label: str | None,
):
    """Open ROOT file, find keys, and draw specified plots."""
    root_file = ROOT.TFile.Open(str(file_path))
    keys = [k.GetName() for k in root_file.GetListOfKeys()]

    if label is not None:
        text = ROOT.TText(0.12, 0.95, label)
        text.SetNDC()
    for key in sorted(plots):
        canvas.Clear()

        obj = root_file.Get(key)
        assert key in keys, f"Object {key} not found in file"
        obj.Draw()

        if label is not None:
            text.Draw("same")

        canvas.Print(f"{output_file}", obj.GetTitle())  # add page
        canvas.Print(str(output_dir / f"{key}.pdf"))


def main(
    finding: Annotated[Path, typer.Argument(dir_okay=False, exists=True)],
    fitting: Annotated[Path, typer.Argument(dir_okay=False, exists=True)],
    output: Annotated[Path, typer.Argument(file_okay=False)],
    finding_plots: Annotated[str, typer.Option(..., parser=split(","))],
    fitting_plots: Annotated[str, typer.Option(..., parser=split(","))],
    label: str | None = None,
):
    print(finding, fitting)

    output.mkdir(parents=True, exist_ok=True)

    ROOT.gStyle.SetOptStat(0)

    outfile = output / f"plots_{output.name}.pdf"

    canvas = ROOT.TCanvas("c", "c", 800, 600)
    canvas.Print(f"{outfile}[")  # open pdf

    plot_from_file(
        finding, finding_plots, output, output_file=outfile, canvas=canvas, label=label
    )
    plot_from_file(
        fitting, fitting_plots, output, output_file=outfile, canvas=canvas, label=label
    )

    canvas.Print(f"{outfile}]")  # close pdf


typer.run(main)
