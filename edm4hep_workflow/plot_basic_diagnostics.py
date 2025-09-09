#!/usr/bin/env python3

from pathlib import Path
from typing import Annotated, Callable
import typer

import ROOT


def split(delim: str) -> Callable[[str], list[str]]:
    def fn(arg: str):
        return arg.split(delim)

    return fn


def main(
    finding: Annotated[Path, typer.Argument(dir_okay=False, exists=True)],
    fitting: Annotated[Path, typer.Argument(dir_okay=False, exists=True)],
    output: Annotated[Path, typer.Argument(file_okay=False)],
    finding_plots: Annotated[str, typer.Option(..., parser=split(","))],
    fitting_plots: Annotated[str, typer.Option(..., parser=split(","))],
):
    print(finding, fitting)

    output.mkdir(parents=True, exist_ok=True)

    finding = ROOT.TFile.Open(str(finding))
    finding_keys = [k.GetName() for k in finding.GetListOfKeys()]

    c = ROOT.TCanvas("c", "c", 800, 600)

    for key in finding_plots:
        obj = finding.Get(key)
        assert key in finding_keys, f"Object {key} not found in file"
        obj.Draw()
        c.SaveAs(str(output / f"{key}.pdf"))

    fitting = ROOT.TFile.Open(str(fitting))
    fitting_keys = [k.GetName() for k in fitting.GetListOfKeys()]

    print(fitting_keys)

    for key in fitting_plots:
        obj = fitting.Get(key)
        assert key in fitting_keys, f"Object {key} not found in file"
        obj.Draw()
        c.SaveAs(str(output / f"{key}.pdf"))


typer.run(main)
