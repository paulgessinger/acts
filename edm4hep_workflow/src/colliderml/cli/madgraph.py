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
            r"\s+(?P<val>.*?)\s+=\s+(?P<key>[A-Za-z0-9_]+).*",
            line,
        )
        if m is None:
            continue
        key = m.group("key")
        val = m.group("val")
        result[key] = val

    return result


def update_madgraph_card(file_path: Path, updates: dict[str, str], log_base: Path):
    logger.info("Updating Madgraph card [b]%s[/b] with overrides %s", file_path.relative_to(log_base), updates)

    raw = file_path.read_text()

    existing = set()

    def repl(m: re.Match) -> str:
        prefix, val, infix, key, suffix = m.groups()

        existing.add(key)
        if key in updates:
            logger.debug("Updating key %s: %s -> %s", key, val, updates[key])
            val = updates[key]
        return prefix + str(val) + infix + key + suffix

    updated = re.sub(r"(\s+)(.*?)(\s+=\s+)([A-Za-z0-9_]+)(.*)", repl, raw)

    for key in set(updates.keys()) - existing:
        logger.debug("Adding new key %s: %s", key, updates[key])
        updated += f"\n{updates[key]} = {key} ! added after the fact\n"

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
    for line in raw:
        # skip comment lines
        if re.match(r"^\s*#.*$", line):
            continue

        m = re.match(r"^\s*?(?P<key>[\w:]+)\s*=\s*(?P<val>[^#\n]+?)(?:\s*#.*)?$", line)
        if m is None:
            continue
        key = m.group("key")
        val = m.group("val")
        result[key] = val.strip()

    return result


def update_pythia8_card(file_path: Path, updates: dict[str, str], log_base: Path):
    logger.info("Updating Pythia8 card [b]%s[/b] with overrides %s", file_path.relative_to(log_base), updates)

    raw = file_path.read_text()

    existing = set()

    def repl(m: re.Match) -> str:
        prefix, key, infix, val, suffix = m.groups()

        existing.add(key)
        logger.debug("matched %s = %s", key, val)
        if key in updates :
            logger.debug("Updating key %s: %s -> %s", key, val, updates[key])
            if val is None:
                suffix = " " + suffix

            val = updates[key]

        if val is None:
            val = ""

        return prefix + key + infix + str(val) + suffix

    ex = r"^(\s*?)([\w:]+)(\s*=\s*)(.*?)(\s?#\s.*)"
    updated = re.sub(ex, repl, raw, flags=re.MULTILINE)

    for key in set(updates.keys()) - existing:
        logger.debug("Adding new key %s: %s", key, updates[key])
        updated += f"\n{key} = {updates[key]} # added after the fact\n"

    file_path.write_text(updated)


def backup_file(file_path: Path, suffix: str, log_base: Path):
    new_path = file_path.with_suffix(file_path.suffix + suffix)
    logger.info(
        "Backing up [i]%s[/i] to [i]%s[/i]",
        file_path.relative_to(log_base),
        new_path.relative_to(log_base),
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


def apply_card_customizations(sample_config: SampleConfig, output_dir: Path):
    cards_dir = output_dir / "Cards"
    if not cards_dir.exists():
        raise FileNotFoundError(f"Could not find cards directory at {cards_dir}")

    run_card_path = cards_dir / "run_card.dat"
    run_data_content = parse_madgraph_card(run_card_path)
    logger.debug("run data settings: %s", run_data_content)

    backup_file(run_card_path, ".orig", log_base=output_dir)
    update_madgraph_card(
        run_card_path,
        sample_config.card_customizations.run_card,
        log_base=output_dir,
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

        backup_file(shower_card_path, ".orig", log_base=output_dir)

        update_pythia8_card(
            shower_card_path,
            sample_config.card_customizations.shower_card,
            log_base=output_dir,
        )
    elif sample_config.run_mode == RunMode.lo_mlm:
        logger.info("Applying Pythia8 card customizations for run mode [b]%s[/b]", sample_config.run_mode)
        pythia8_card_path: Path|None = None
        for candidate in ["pythia8_card.dat", "pythia8_card_default.dat"]:
            candidate_path = cards_dir / candidate
            if candidate_path.exists():
                pythia8_card_path = candidate_path
                break
        if pythia8_card_path is None:
            raise FileNotFoundError("Found no pythia8_card.dat in Cards/")

        pythia8_data_content = parse_pythia8_card(pythia8_card_path)
        logger.debug("pythia8 data settings: %s", pythia8_data_content)

        backup_file(pythia8_card_path, ".orig", log_base=output_dir)

        update_pythia8_card(
            pythia8_card_path,
            sample_config.card_customizations.pythia8_card,
            log_base=output_dir,
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

        return

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
                logger.info(
                    f"Running grid/envelope compilation via {" ".join(cmd)}"
                )
                subprocess.run(cmd, cwd=output_dir, check=True)

        elif sample_config.run_mode == RunMode.lo_mlm:
            # Daniel reports it's not possible to precompile this, but lets see
            logger.info("Build Grids/Envelopes (zero events)")

            with with_madgraph_card_customization(
                output_dir / "Cards" / "run_card.dat",
                {"nevents": "0", "req_acc": "0.001"},
                log_base=output_dir,
            ):

                cmd = [str(generate_events_exe), "-f", "--name", "run_build"]
                logger.info(
                    f"Running grid/envelope compilation via {" ".join(cmd)}"
                )
                subprocess.run(cmd, cwd=output_dir, check=True)

        # @TODO: madgraph might already be creating a tarball in amcatnlo.tar.gz
        logger.info("Packing output to %s", output)
        with tarfile.open(output, "w:gz") as tar:
            tar.add(scratch_dir, arcname="")


@app.command()
def generate():
    print("hi")
