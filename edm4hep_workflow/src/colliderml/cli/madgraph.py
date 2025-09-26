#!/usr/bin/env python3

import contextlib
import os
import logging
from pathlib import Path
from typing import Annotated
import subprocess
import tempfile
from string import Template
import tarfile
from rich import print
import shutil
import re

import typer

from colliderml.config import RunMode, SampleConfig
import colliderml.logging

app = typer.Typer()

logger = colliderml.logging.get_logger(__name__)

MG_OUTPUT_DIR = "proc_output_mg"


def format_dict_for_logging(d: dict[str, str]) -> str:
    """Format dictionary for rich logging with bold keys and bold blue values."""
    return ", ".join(f"[b]{k}[/b]=[bold blue]{v}[/bold blue]" for k, v in d.items())


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
        sample_config.card_customizations.run_card,
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


def apply_run_configuration(events: int, seed: int, process_dir: Path):
    logger.info(
        "Setting up run configuration: [b]events=%d[/b], [b]seed=%d[/b]", events, seed
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
        {"nevents": str(events), "iseed": str(seed)},
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
            """
import model ${model}
${definitions}
${generate_command}
output ${output_dir} -f
exit
                       """.strip()
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


SEED_DEFAULT = 42


@app.command()
def generate(
    sample_file: Annotated[Path, typer.Argument(..., exists=True, dir_okay=False)],
    tarball: Annotated[Path, typer.Argument(..., exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Argument(..., file_okay=False)] = Path.cwd(),
    scratch_dir: Path | None = None,
    events: Annotated[int, typer.Option("--events", "-n")] = 10,
    seed: int = SEED_DEFAULT,
):
    logger.info("Loading sample config from %s", sample_file)
    sample_config = SampleConfig.load(sample_file)
    print(sample_config)
    logger.info("Label: %s", sample_config.label)
    logger.info("Model: %s", sample_config.model)
    logger.info("Process:\n%s", sample_config.generate_command)
    logger.info("Run mode: %s", sample_config.run_mode)

    if seed == SEED_DEFAULT:
        logger.warning(
            "Default seed %d given. You will likely use a different seed!", SEED_DEFAULT
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

        apply_run_configuration(events=events, seed=seed, process_dir=process_dir)

        if sample_config.run_mode == RunMode.nlo_fxfx:
            logger.info("Launching event generation for [b]nlo_fxfx[/b] (precompiled)")

            generate_events_exe = (process_dir / "bin" / "generate_events").resolve()
            if not generate_events_exe.exists():
                logger.error("Could not find executable %s", generate_events_exe)
                raise typer.Exit(1)

            cmd = [
                str(generate_events_exe),
                "--only_generation",  # do not reproduce grids/envelopes
                "--nocompile",  # do not recompile
                "--force",
            ]
            logger.info("Running event generation with command %s", " ".join(cmd))
            subprocess.run(cmd, cwd=process_dir, check=True)

        elif sample_config.run_mode == RunMode.lo_mlm:
            raise NotImplementedError("LO MLM event generation not implemented yet")
