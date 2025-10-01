from typing import Annotated

import typer

EVENTS = Annotated[int, typer.Option("--events", "-n")]
JOBS = Annotated[int, typer.Option("--jobs", "-j")]
