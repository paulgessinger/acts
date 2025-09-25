#!/usr/bin/env python3

from contextlib import ExitStack
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
    logger.info("Updating Madgraph card [b]%s[/b]", file_path.relative_to(log_base))

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

        m = re.match(r"^\s*?(?P<key>\w+)\s*=\s*(?P<val>[^#\n]+)\b", line)
        if m is None:
            continue
        key = m.group("key")
        val = m.group("val")
        result[key] = val.strip()

    return result


def update_pythia8_card(file_path: Path, updates: dict[str, str], log_base: Path):
    logger.info("Updating Pythia8 card [b]%s[/b]", file_path.relative_to(log_base))

    raw = file_path.read_text()

    existing = set()

    def repl(m: re.Match) -> str:
        prefix, key, infix, val, suffix = m.groups()

        existing.add(key)
        if key in updates:
            logger.debug("Updating key %s: %s -> %s", key, val, updates[key])
            if val is None:
                suffix = " " + suffix

            val = updates[key]

        if val is None:
            val = ""

        return prefix + key + infix + str(val) + suffix

    ex = r"^(\s*?)(\w+)(\s*=\s*)(.*?)(\s?#\s.*)"
    print(re.findall(ex, raw))
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
        raise NotImplementedError("LO MLM not implemented yet")


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

    with ExitStack() as stack:
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

        apply_card_customizations(sample_config, output_dir)

        # @TODO: madgraph might already be creating a tarball in amcatnlo.tar.gz
        logger.info("Packing output to %s", output)
        with tarfile.open(output, "w:gz") as tar:
            tar.add(scratch_dir, arcname="")


@app.command()
def generate():
    print("hi")
