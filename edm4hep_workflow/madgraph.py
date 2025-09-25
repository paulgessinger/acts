#!/usr/bin/env python3

from pathlib import Path
from typing import Annotated

import typer

from config import SampleConfig

app = typer.Typer()


@app.command()
def init(
    sample_file: Annotated[Path, typer.Argument(..., exists=True, dir_okay=False)]
):
    sample_config = SampleConfig.load(sample_file)
    print(sample_config)


@app.command()
def generate():
    print("hi")


if __name__ == "__main__":
    app()
