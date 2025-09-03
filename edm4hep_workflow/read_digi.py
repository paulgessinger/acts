#!/usr/bin/env python3
from pathlib import Path
import typer
import acts


def main(
    input: Path,
    jobs: int = -1,
    logLevel: str = "INFO",
):
    from acts import UnitConstants as u
    import acts.examples
    from acts.examples import Sequencer
    from acts.examples.edm4hep import (
        PodioReader,
    )

    logLevel = getattr(acts.logging, logLevel.upper(), acts.logging.INFO)

    s = Sequencer(numThreads=jobs)
    s.config.logLevel = acts.logging.DEBUG

    podioReader = PodioReader(
        level=logLevel,
        inputPath=str(input),
        outputFrame="events",
        category="events",
    )
    s.addReader(podioReader)

    s.run()


typer.run(main)
