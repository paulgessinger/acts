#!/usr/bin/env python3

# /// script
# dependencies = [
#   "typer",
#   "rich",
#   "pandas",
# ]
# ///

import typer
from typing import Annotated
from pathlib import Path
import pandas
from rich.console import Console
from rich.table import Table


def main(file: Annotated[Path, typer.Argument(dir_okay=True, exists=True)]):
    df = pandas.read_csv(file)
    df = df.sort_values(by="time_total_s", ascending=False)

    df["rel"] = df.time_total_s / df.time_total_s.sum()

    table = Table(title="Timing breakdown", expand=True)

    table.add_column("Component", justify="right", style="cyan")
    table.add_column(r"Total time \[s]", justify="right", style="magenta")
    table.add_column(r"% of total", justify="right", style="magenta")
    table.add_column(r"Time per event \[ms]", justify="right", style="green")

    for row in df.itertuples():
        table.add_row(
            row.identifier,
            f"{row.time_total_s:> 8.2f}",
            f"{row.rel*100:> 3.2f}",
            f"{row.time_perevent_s*1000.:> 8.2f}",
        )

    table.add_section()
    table.add_row(
        "TOTAL",
        f"{df.time_total_s.sum():> 8.2f}",
        "100",
        f"{df.time_perevent_s.sum()*1000:> 8.2f}",
        style="green bold",
    )

    console = Console()
    console.print(table)


typer.run(main)
