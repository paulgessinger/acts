from typing import Annotated
import logging

import typer

import colliderml.cli.madgraph
from colliderml.logging import configure_logging

app = typer.Typer()


@app.callback()
def main(verbose: Annotated[int, typer.Option("--verbose", "-v", count=True)] = 0):
    level = logging.INFO
    if verbose >= 2:
        level = logging.DEBUG
    configure_logging(level)
    logger = logging.getLogger(__name__)
    logger.debug("Logging configured, level=%s", level)


app.add_typer(colliderml.cli.madgraph.app, name="madgraph")
