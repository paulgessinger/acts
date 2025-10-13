from typing import Annotated
import logging

import typer

import colliderml.cli.madgraph
import colliderml.cli.particle_gun
import colliderml.cli.pythia8
import colliderml.cli.simulation
import colliderml.cli.digitization
import colliderml.cli.reconstruction
import colliderml.cli.tree_size
import colliderml.cli.plot_basic_diagnostics
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
app.command("particle-gun")(colliderml.cli.particle_gun.main)
app.command("pythia8")(colliderml.cli.pythia8.main)
app.command("simulation")(colliderml.cli.simulation.main)
app.command("digitization")(colliderml.cli.digitization.main)
app.command("tree-size")(colliderml.cli.tree_size.main)
app.command("reconstruction")(colliderml.cli.reconstruction.main)

app.add_typer(colliderml.cli.plot_basic_diagnostics.app, name="plot")
