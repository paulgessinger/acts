#!/usr/bin/env python3

import contextlib
import os
import logging
import textwrap
from pathlib import Path
from typing import Annotated, Generator
import subprocess
import tempfile
from string import Template
import tarfile
import shutil
import re
import collections

import typer

from colliderml.config import RunMode, MadgraphSampleConfig as SampleConfig
from colliderml import constants
import colliderml.logging
import colliderml.util

app = typer.Typer()

logger = colliderml.logging.get_logger(__name__)

MG_OUTPUT_DIR = "proc_output_mg"


def format_dict_for_logging(d: dict[str, str]) -> str:
    """Format dictionary for rich logging with bold keys and bold blue values."""
    return ", ".join(
        f"[bold blue]{k}[/bold blue]=[bold green]{v}[/bold green]" for k, v in d.items()
    )


class Madgraph:
    _exe: Path

    def __init__(self):
        self._exe = Path(
            subprocess.run(
                ["which", "mg5_aMC"], capture_output=True, text=True, check=True
            ).stdout.strip()
        )
        self._exe = self._exe.resolve()
        logger.debug("Found Madgraph executable at %s", self._exe)

    def run(self, script: Path, **kwargs):
        subprocess.run([self._exe, script], check=True, **kwargs)

    def run_iter(self, script: Path, **kwargs) -> Generator[str, None, None]:
        yield from colliderml.util.stream_subprocess(
            [str(self._exe), str(script)], **kwargs
        )


def parse_madgraph_card(path_or_string: Path | str) -> dict[str, str]:
    result = {}
    if isinstance(path_or_string, str):
        raw = path_or_string.splitlines()
    elif isinstance(path_or_string, Path):
        raw = path_or_string.read_text().splitlines()
    else:
        raise TypeError("path_or_string must be a Path or str")

    for line in raw:
        m = re.match(
            r"\s+(?P<val>.*?)\s+=\s+(?P<key>[A-Za-z0-9_:]+).*",
            line,
        )
        if m is None:
            continue
        key = m.group("key")
        val = m.group("val")
        result[key] = val

    return result


def update_madgraph_card(file_path: Path, updates: dict[str, str], log_base: Path):
    logger.info(
        "Updating Madgraph card [b]%s[/b] with overrides %s",
        file_path.relative_to(log_base),
        format_dict_for_logging(updates),
        extra={"highlighter": False},
    )

    raw = file_path.read_text()

    existing = set()

    def repl(m: re.Match) -> str:
        prefix, val, infix, key, suffix = m.groups()

        existing.add(key)
        if key in updates and val != updates[key]:
            logger.debug("Updating key %s: %s -> %s", key, val, updates[key])
            val = updates[key]
        return prefix + str(val) + infix + key + suffix

    updated = re.sub(r"(\s+)(.*?)(\s+=\s+)([A-Za-z0-9_:]+)(.*)", repl, raw)

    for key in sorted(set(updates.keys()) - existing):
        logger.debug("Adding new key %s = %s", key, updates[key])
        updated += f"\n{updates[key]} = {key}\n"

    file_path.write_text(updated)


@contextlib.contextmanager
def with_madgraph_card_customization(
    file_path: Path, updates: dict[str, str], log_base: Path
):
    logger.info(
        "Applying temporary Madgraph card customizations to [b]%s[/b]",
        file_path.relative_to(log_base),
    )
    original = file_path.read_text()
    update_madgraph_card(file_path, updates, log_base)

    try:
        yield
    finally:
        logger.info(
            "Restoring original Madgraph card [b]%s[/b]",
            file_path.relative_to(log_base),
        )
        file_path.write_text(original)


def parse_pythia8_card(path_or_string: Path | str) -> dict[str, str]:
    if isinstance(path_or_string, str):
        raw = path_or_string.splitlines()
    elif isinstance(path_or_string, Path):
        raw = path_or_string.read_text().splitlines()

    result = {}
    for line_num, line in enumerate(raw, 1):
        # skip comment lines
        if re.match(r"^\s*[!#].*$", line):
            continue

        m = re.match(r"^\s*?(?P<key>[\w:]+)\s*=\s*(?P<val>[^#\n]+?)(?:\s*#.*)?$", line)
        if m is None:
            continue
        key = m.group("key")
        val = m.group("val").strip()

        # Check for duplicate keys
        if key in result:
            file_info = (
                f" in {path_or_string}" if isinstance(path_or_string, Path) else ""
            )
            raise ValueError(
                f"Duplicate key '{key}' found on line {line_num}{file_info}. "
                f"Previous value: '{result[key]}', new value: '{val}'"
            )

        result[key] = val

    return result


def update_pythia8_card(file_path: Path, updates: dict[str, str], log_base: Path):
    logger.info(
        "Updating Pythia8 card [b]%s[/b] with overrides %s",
        file_path.relative_to(log_base),
        format_dict_for_logging(updates),
        extra={"highlighter": False},
    )

    raw = file_path.read_text()

    existing = set()

    def repl(m: re.Match) -> str:
        prefix, key, infix, val, suffix = m.groups()

        existing.add(key)
        if key in updates and val != updates[key]:
            logger.debug("Updating key %s: %s -> %s", key, val, updates[key])
            val = updates[key]

        if val is None:
            val = ""

        if suffix is None:
            suffix = ""

        return prefix + key + infix + str(val) + suffix

    ex = r"^(\s*?)([\w:]+)(\s*=\s*)(.*?)(\s*#.*)?$"
    updated = re.sub(ex, repl, raw, flags=re.MULTILINE)

    for key in sorted(set(updates.keys()) - existing):
        logger.debug("Adding new key %s = %s", key, updates[key])
        updated += f"\n{key} = {updates[key]}\n"

    file_path.write_text(updated)


def backup_file(file_path: Path, suffix: str, log_base: Path):
    new_path = file_path.with_suffix(file_path.suffix + suffix)

    # If target backup already exists, find the next available numbered backup
    if new_path.exists():
        counter = 1
        while True:
            numbered_path = file_path.with_suffix(
                file_path.suffix + suffix + f".{counter}"
            )
            if not numbered_path.exists():
                new_path = numbered_path
                break
            counter += 1

    logger.info(
        "Backing up [i]%s[/i] to [i]%s[/i]",
        file_path.relative_to(log_base),
        new_path.relative_to(log_base),
        extra={"highlighter": False},
    )
    shutil.copy(file_path, new_path)


def fix_fortran_dollar_syntax(file_path: Path, log_base: Path):
    """
    Fix non-standard Fortran $ format descriptor with standard advance='no'.

    Replaces patterns like:
    write(26,'(x,i1,a2$)') with write(26,'(x,i1,a2)',advance='no')
    """
    if not file_path.exists():
        logger.warning("Fortran file not found: %s", file_path.relative_to(log_base))
        return

    logger.info("Fixing Fortran syntax in [b]%s[/b]", file_path.relative_to(log_base))

    content = file_path.read_text()
    original_content = content

    # Pattern to match write statements with $ format descriptor
    # Matches: write(unit,'(format$)') args
    pattern = r"write\(([^,]+),\s*'\(([^']*)\$\)'\)\s*([^\n]*)"

    def replace_dollar_format(match):
        unit = match.group(1)
        format_str = match.group(2)  # format without the $
        args = match.group(3)

        # Construct the replacement with advance='no'
        return f"write({unit},'({format_str})',advance='no') {args}"

    content = re.sub(pattern, replace_dollar_format, content)

    if content != original_content:
        # Create backup before modifying
        backup_file(file_path, ".orig", log_base)
        file_path.write_text(content)
        logger.debug("Fixed non-standard Fortran $ syntax in %s", file_path.name)
    else:
        logger.debug("No Fortran $ syntax found in %s", file_path.name)


def fix_madgraph_fortran(output_dir: Path):
    """Fix non-standard Fortran syntax in MadGraph generated files."""
    logger.info("Fixing non-standard Fortran syntax in MadGraph files")

    # Files that need Fortran fixes
    fortran_files = [
        "SubProcesses/symmetry_fks_v3.f",
        "SubProcesses/write_ajob.f",
    ]

    for file_rel_path in fortran_files:
        file_path = output_dir / file_rel_path
        if not file_path.exists():
            continue
        fix_fortran_dollar_syntax(file_path, output_dir)


def apply_card_customizations(sample_config: SampleConfig, process_dir: Path):
    cards_dir = process_dir / "Cards"
    if not cards_dir.exists():
        logger.error("Could not find cards directory at %s", cards_dir)
        raise typer.Exit(1)

    run_card_path = cards_dir / "run_card.dat"
    run_data_content = parse_madgraph_card(run_card_path)
    logger.debug("run data settings: %s", run_data_content)

    backup_file(run_card_path, ".orig", log_base=process_dir)
    update_madgraph_card(
        run_card_path,
        sample_config.run_card,
        log_base=process_dir,
    )

    if sample_config.run_mode == RunMode.nlo_fxfx:
        shower_card_path = cards_dir / "shower_card.dat"
        logger.info(
            "Applying shower card customizations for run mode [b]%s[/b]",
            sample_config.run_mode,
        )
        if not shower_card_path.exists():
            raise FileNotFoundError(f"Could not find shower card at {shower_card_path}")
        shower_data_content = parse_pythia8_card(shower_card_path)
        logger.debug("shower data settings: %s", shower_data_content)

        backup_file(shower_card_path, ".orig", log_base=process_dir)

        update_pythia8_card(
            shower_card_path,
            sample_config.card_customizations.shower_card,
            log_base=process_dir,
        )
    elif sample_config.run_mode == RunMode.lo_mlm:
        logger.info(
            "Applying Pythia8 card customizations for run mode [b]%s[/b]",
            sample_config.run_mode,
        )
        pythia8_card_path: Path | None = None
        for candidate in ["pythia8_card.dat", "pythia8_card_default.dat"]:
            candidate_path = cards_dir / candidate
            if candidate_path.exists():
                pythia8_card_path = candidate_path
                break
        if pythia8_card_path is None:
            raise FileNotFoundError("Found no pythia8_card.dat in Cards/")

        pythia8_data_content = parse_pythia8_card(pythia8_card_path)
        logger.debug("pythia8 data settings: %s", pythia8_data_content)

        backup_file(pythia8_card_path, ".orig", log_base=process_dir)

        update_pythia8_card(
            pythia8_card_path,
            sample_config.card_customizations.pythia8_card,
            log_base=process_dir,
        )


def apply_run_configuration(
    events: int,
    seed: int,
    process_dir: Path,
    margin: float = constants.MG_GENERATION_REL_MARGIN,
):
    logger.info(
        "Setting up run configuration: %s",
        format_dict_for_logging({"events": str(events), "seed": str(seed)}),
    )

    cards_dir = process_dir / "Cards"
    if not cards_dir.exists():
        logger.error("Could not find cards directory at %s", cards_dir)
        raise typer.Exit(1)

    run_card_path = cards_dir / "run_card.dat"
    if not run_card_path.exists():
        logger.error("Could not find run card at %s", run_card_path)
        raise typer.Exit(1)

    backup_file(run_card_path, ".orig", log_base=process_dir)
    update_madgraph_card(
        run_card_path,
        {"nevents": str(int(events * (1 + margin))), "iseed": str(seed)},
        log_base=process_dir,
    )

    pythia8_cards = [
        cards_dir / "pythia8_card.dat",
        cards_dir / "pythia8_card_default.dat",
    ]

    for pythia8_card in pythia8_cards:
        if not pythia8_card.exists():
            logger.debug("No Pythia8 card found at %s", pythia8_card)
            continue

        backup_file(pythia8_card, ".orig", log_base=process_dir)
        update_pythia8_card(
            pythia8_card,
            {"Random:seed": str(seed), "Main:numberOfEvents": str(events)},
            log_base=process_dir,
        )

    shower_card_path = cards_dir / "shower_card.dat"
    if shower_card_path.exists():
        backup_file(shower_card_path, ".orig", log_base=process_dir)
        update_pythia8_card(
            shower_card_path,
            {"nevents": str(events), "rnd_seed": str(seed)},
            log_base=process_dir,
        )


def parse_number_of_merged_events(output: str) -> int:
    m = re.search(r".*Nb of events after merging :  (\d+)", output)
    if m is None:
        logger.error(
            "Could not determine number of events after merging from Madgraph output"
        )
        raise typer.Exit(1)

    return int(m.group(1))


def locate_single_output(process_dir: Path, name: str) -> Path:
    events_dir = process_dir / "Events"
    if not events_dir.exists():
        logger.error("Could not find Events/ directory in %s", process_dir)
        raise typer.Exit(1)

    run_dir = events_dir / name
    if not run_dir.exists():
        logger.error("Could not find %s/ directory in Events/", name)
        raise typer.Exit(1)

    events_file = list(run_dir.glob("*.hepmc.gz"))

    if len(events_file) == 0:
        logger.error("Could not find any *.hepmc.gz files in %s", run_dir)
        raise typer.Exit(1)

    if len(events_file) > 1:
        logger.error("Found multiple *.hepmc.gz files in %s", run_dir)
        raise typer.Exit(1)

    return events_file[0]


def run_madgraph_via_script(
    process_dir: Path,
    events: int,
    seed: int,
    name: str,
    jobs: int,
    margin: float = constants.MG_GENERATION_REL_MARGIN,
) -> tuple[int, Path]:
    mg = Madgraph()

    scratch_dir = process_dir.parent

    apply_run_configuration(events=events, seed=seed, process_dir=process_dir)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".mg5") as script_file:
        script_file_path = Path(script_file.name).resolve()

        script_lines = [
            "help launch",
            f"launch {process_dir} --name={name}",
            "shower=PYTHIA8",
            # f"set nb_cores {jobs}",
            f"set nevents {int(events * (1+margin))}",
            f"set iseed {seed}",
        ]

        script = "\n".join(script_lines)

        logger.debug(
            "Writing madgraph script %s to execute:\n%s",
            script_file.name,
            script,
        )

        script_file.write(script)
        script_file.flush()

        lines = collections.deque(maxlen=20)
        for line in mg.run_iter(script_file_path, cwd=scratch_dir):
            print(line, end="")
            lines.append(line)

        merged_events = parse_number_of_merged_events("".join(lines))
        logger.info("Number of events after merging: %d", merged_events)

        events_file = locate_single_output(process_dir, name)

        return (merged_events, events_file)


@app.command()
def init(
    sample_file: Annotated[Path, typer.Argument(..., exists=True, dir_okay=False)],
    output: Path,
    scratch_dir: Path | None = None,
):
    mg = Madgraph()
    logger.info("Loading sample config from %s", sample_file)
    sample_config = SampleConfig.load(sample_file)
    print(sample_config)
    logger.info("Label: %s", sample_config.label)
    logger.info("Model: %s", sample_config.model)
    logger.info("Process:\n%s", sample_config.generate_command)
    logger.info("Run mode: %s", sample_config.run_mode)

    with contextlib.ExitStack() as stack:
        if scratch_dir is None:
            scratch_dir = Path(stack.enter_context(tempfile.TemporaryDirectory()))

        scratch_dir = scratch_dir.resolve()

        if not scratch_dir.exists():
            scratch_dir.mkdir(parents=True)

        logger.debug("Working in scratch directory %s", scratch_dir)

        tpl = Template(
            textwrap.dedent(
                """
                import model ${model}
                ${definitions}
                ${generate_command}
                output ${output_dir} -f
                exit
                """.strip()
            )
        )

        mg_script = tpl.substitute(
            {**sample_config.model_dump(), "output_dir": MG_OUTPUT_DIR}
        )
        logger.debug("Generating Madgraph script:\n%s", mg_script)

        script_file = scratch_dir / "process_script.mg5"
        script_file.write_text(mg_script)

        logger.info("Running Madgraph...")
        mg.run(script_file, cwd=scratch_dir)
        logger.info(
            "Applying run card and shower card customizations for run mode %s",
            sample_config.run_mode,
        )
        output_dir = scratch_dir / MG_OUTPUT_DIR

        # Fix non-standard Fortran syntax
        fix_madgraph_fortran(output_dir)

        apply_card_customizations(sample_config, output_dir)

        generate_events_exe = output_dir / "bin" / "generate_events"
        if not generate_events_exe.exists():
            raise FileNotFoundError(generate_events_exe)

        if sample_config.run_mode == RunMode.nlo_fxfx:
            logger.info("Build Grids/Envelopes (zero events)")

            with with_madgraph_card_customization(
                output_dir / "Cards" / "run_card.dat",
                {"nevents": "0", "req_acc": "0.001"},
                log_base=output_dir,
            ):
                cmd = [str(generate_events_exe), "-f", "--name", "run_build"]
                logger.info(f"Running grid/envelope compilation via {' '.join(cmd)}")
                subprocess.run(cmd, cwd=output_dir, check=True)

        elif sample_config.run_mode == RunMode.lo_mlm:
            # @TODO: Check if there's any way to do this
            logger.info("Not precompiling anything for LO MLM mode")

        # @TODO: madgraph might already be creating a tarball in amcatnlo.tar.gz
        logger.info("Packing output to %s", output)
        with tarfile.open(output, "w:gz") as tar:
            tar.add(scratch_dir, arcname="")


@app.command()
def generate(
    sample_file: Annotated[Path, typer.Argument(..., exists=True, dir_okay=False)],
    tarball: Annotated[Path, typer.Argument(..., exists=True, dir_okay=False)],
    output: Path,
    events_per_file: int | None = None,
    scratch_dir: Path | None = None,
    events: Annotated[int, typer.Option("--events", "-n")] = 10,
    seed: int = constants.SEED_DEFAULT,
    force: Annotated[bool, typer.Option("--force", "-f")] = False,
    max_iterations: int = constants.GENERATION_ITER_MAX,
    jobs: int = os.cpu_count() or 1,
):
    import acts.examples.hepmc3

    logger.info("Loading sample config from %s", sample_file)
    sample_config = SampleConfig.load(sample_file)
    print(sample_config)
    logger.info("Label: %s", sample_config.label)
    logger.info("Model: %s", sample_config.model)
    logger.info(
        "Process:\n[italic]%s[/italic]",
        sample_config.generate_command,
        extra={"highlighter": False},
    )
    logger.info("Run mode: %s", sample_config.run_mode)

    if output.exists():
        if force:
            logger.warning(
                "Output file %s already exists, overwriting due to --force", output
            )
            output.unlink()
        else:
            logger.error(
                "Output file %s already exists, use --force to overwrite", output
            )
            raise typer.Exit(1)

    if seed == constants.SEED_DEFAULT:
        logger.warning(
            "Default seed %d given. You will likely use a different seed!",
            constants.SEED_DEFAULT,
        )

    with contextlib.ExitStack() as stack:
        if scratch_dir is None:
            scratch_dir = Path(stack.enter_context(tempfile.TemporaryDirectory()))

        scratch_dir = scratch_dir.resolve()

        if not scratch_dir.exists():
            scratch_dir.mkdir(parents=True)

        logger.debug("Working in scratch directory %s", scratch_dir)

        logger.info("Unpacking Madgraph tarball %s", tarball)

        with tarfile.open(tarball, "r:gz") as tar:
            tar.extractall(path=scratch_dir)

        process_dir = scratch_dir / MG_OUTPUT_DIR
        if not process_dir.exists():
            logger.error(
                "Could not find extracted Madgraph output directory %s", process_dir
            )
            raise typer.Exit(1)

        apply_card_customizations(sample_config, process_dir)

        if sample_config.run_mode == RunMode.nlo_fxfx:
            logger.info("Launching event generation for [b]nlo_fxfx[/b] (precompiled)")

            apply_run_configuration(
                events=events,
                seed=seed,
                process_dir=process_dir,
                margin=0,  # we don't seem to need margin here
            )

            generate_events_exe = (process_dir / "bin" / "generate_events").resolve()
            if not generate_events_exe.exists():
                logger.error("Could not find executable %s", generate_events_exe)
                raise typer.Exit(1)

            run_name = "run_01"
            cmd = [
                str(generate_events_exe),
                "--only_generation",  # do not reproduce grids/envelopes
                "--nocompile",  # do not recompile
                "--force",
                "--name",
                run_name,
            ]
            logger.info("Running event generation with command %s", " ".join(cmd))
            subprocess.run(cmd, cwd=process_dir, check=True)

            events_file = locate_single_output(process_dir, run_name)

            logger.info("Normalizing output to %s", output)

            acts.examples.hepmc3.normalizeFiles(
                inputFiles=[events_file],
                singleOutputPath=output,
                maxEvents=events,
                format=acts.examples.hepmc3.formatFromFilename(output),
                compression=acts.examples.hepmc3.compressionFromFilename(output),
                compressionLevel=9,
                verbose=True,
            )

        elif sample_config.run_mode == RunMode.lo_mlm:
            logger.info(
                "Launching event generation for [b]lo_mlm[/b]",
            )

            requested_events = 0
            produced_events = 0

            events_collect_dir = process_dir / "Events_collect"
            events_collect_dir.mkdir(exist_ok=True)

            events_dir = process_dir / "Events"

            if events_dir.exists():
                logger.debug("Cleaning existing Events/ directory at %s", events_dir)
                shutil.rmtree(events_dir)
                events_dir.mkdir()

            run_files = []

            for i in range(1, max_iterations + 1):
                logger.info("Iteration %d to produce events", i)

                run_name = f"run_iter_{i:02d}"

                if i == 1:
                    # first run, let's use the target event count
                    if sample_config.nevents_scale_factor != 1.0:
                        logger.info(
                            "Using scale factor %f to speculatively generate more events in the first iteration",
                            sample_config.nevents_scale_factor,
                        )
                    run_events = int(events * sample_config.nevents_scale_factor)
                else:
                    # subsequent runs, estimate veto rate from existing
                    veto_pass_rate = produced_events / requested_events
                    logger.info(
                        f"Estimated veto pass rate: {veto_pass_rate * 100:.2f}%"
                    )
                    run_events = int(
                        (events - produced_events)
                        / veto_pass_rate
                        * (1 + constants.GENERATED_EVENT_MARGIN_REL)
                    )

                logger.info("Requesting %d events in this iteration", run_events)

                merged_events, output_file = run_madgraph_via_script(
                    process_dir, run_events, seed * i, run_name, jobs
                )

                run_files.append(output_file)

                requested_events += run_events
                produced_events += merged_events

                logger.info(
                    "Have produced %d events so far (%d / %d)",
                    produced_events,
                    i,
                    constants.GENERATION_ITER_MAX,
                )

                if produced_events >= events:
                    logger.info("Have produced enough events, stopping")
                    break

            logger.info(
                "Combining %d output files into %s (events per file: %s)",
                len(run_files),
                output,
                events_per_file,
            )
            logger.debug("Files to combine:\n%s", "\n".join(map(str, run_files)))

            acts.examples.hepmc3.normalizeFiles(
                inputFiles=run_files,
                singleOutputPath=output,
                maxEvents=events,
                eventsPerFile=events_per_file,
                format=acts.examples.hepmc3.formatFromFilename(output),
                compression=acts.examples.hepmc3.compressionFromFilename(output),
                compressionLevel=9,
                verbose=True,
            )
