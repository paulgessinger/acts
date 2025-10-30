#!/usr/bin/env python3
# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "rich",
#     "typer",
# ]
# ///

"""Display ROOT TTree branch sizes in a sorted, human-readable table."""

import os
import re
import sys
import tempfile
from dataclasses import dataclass
from enum import Enum
from typing import Annotated
from tabulate import tabulate

import typer
import numpy

import colliderml.logging

try:
    import ROOT
except ImportError:
    print("Error: PyROOT is not available", file=sys.stderr)
    sys.exit(1)

from rich.console import Console
from rich.table import Table

from colliderml.util import human_readable_size


class SortBy(str, Enum):
    """Sort options for branch display."""

    file_size = "file_size"
    total_size = "total_size"
    compression = "compression"
    name = "name"


@dataclass
class BranchInfo:
    """Information about a TTree branch."""

    name: str
    type_info: str
    entries: int
    total_size: int
    file_size: int
    compression: float
    avg_size_per_entry: float = 0.0


def parse_bytes(size_str: str) -> int:
    """Parse byte size from string like '829584 bytes'."""
    match = re.search(r"(\d+)\s*bytes", size_str)
    return int(match.group(1)) if match else 0


def parse_compression(comp_str: str) -> float:
    """Parse compression ratio from string."""
    match = re.search(r"Compression=\s*([\d.]+)", comp_str)
    return float(match.group(1)) if match else 1.0


def parse_entries(entries_str: str) -> int:
    """Parse entries count from string like 'Entries :      100'."""
    match = re.search(r"Entries\s*:\s*(\d+)", entries_str)
    return int(match.group(1)) if match else 0


def human_readable_size_float(size_bytes: float) -> str:
    """Convert bytes to human-readable format (float version)."""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"


def parse_tree_print(tree) -> list[BranchInfo]:
    """Parse TTree::Print() output and extract branch information."""
    # Redirect ROOT output to capture it
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".txt", delete=False) as tmp:
        tmp_path = tmp.name

    ROOT.gROOT.ProcessLine(f".> {tmp_path}")
    tree.Print()
    ROOT.gROOT.ProcessLine(".>")

    branches = []
    with open(tmp_path) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Look for branch lines starting with *Br
        if line.startswith("*Br"):
            # Extract branch name and type
            match = re.match(r"\*Br\s+\d+\s*:\s*(.+?)\s*:\s*(.*)\*", line)
            if match:
                branch_name = match.group(1).strip()
                type_info = match.group(2).strip() if match.group(2) else ""

                # Next line might have type information
                if i + 1 < len(lines) and lines[i + 1].strip().startswith("*"):
                    next_line = lines[i + 1].strip()
                    if "|" in next_line:
                        type_match = re.search(r"\|\s*(.+?)\s*\*", next_line)
                        if type_match:
                            type_info = type_match.group(1).strip()
                        i += 1

                # Look for the Entries line with size info
                if i + 1 < len(lines):
                    entries_line = lines[i + 1].strip()
                    if "Entries" in entries_line and "Total" in entries_line:
                        entries = parse_entries(entries_line)
                        total_size = parse_bytes(
                            re.search(
                                r"Total\s+Size=\s*(.+?)\s+File", entries_line
                            ).group(1)
                            if re.search(r"Total\s+Size=\s*(.+?)\s+File", entries_line)
                            else ""
                        )
                        file_size = parse_bytes(entries_line)
                        compression = parse_compression(entries_line)

                        branches.append(
                            BranchInfo(
                                name=branch_name,
                                type_info=type_info,
                                entries=entries,
                                total_size=total_size,
                                file_size=file_size,
                                compression=compression,
                            )
                        )
                        i += 1

        i += 1

    return branches


def compute_average_sizes(tree, branches: list[BranchInfo]) -> None:
    """Compute average size per entry for each branch by iterating through the tree."""
    num_entries = tree.GetEntries()
    if num_entries == 0:
        return

    # Get all branch names
    branch_dict = {b.name: b for b in branches}

    # Iterate through all entries to compute average sizes
    for branch_info in branches:
        branch = tree.GetBranch(branch_info.name)
        if branch:
            # Sum up compressed bytes for this branch across all entries
            total_bytes = 0
            for entry in range(num_entries):
                total_bytes += branch.GetEntry(entry)

            # Calculate average
            branch_info.avg_size_per_entry = (
                total_bytes / num_entries if num_entries > 0 else 0.0
            )


def display_branch_table(
    branches: list[BranchInfo], sort_by: str = "file_size", show_avg_size: bool = False
):
    """Display branch information in a rich table."""
    # Sort branches
    if sort_by == "file_size":
        branches = sorted(branches, key=lambda b: b.file_size, reverse=True)
    elif sort_by == "total_size":
        branches = sorted(branches, key=lambda b: b.total_size, reverse=True)
    elif sort_by == "compression":
        branches = sorted(branches, key=lambda b: b.compression, reverse=True)
    elif sort_by == "name":
        branches = sorted(branches, key=lambda b: b.name)

    # Calculate total file size for percentage calculation
    total_file_size = sum(b.file_size for b in branches)

    table = Table(title="TTree Branch Sizes")
    table.add_column("Branch Name", style="cyan", no_wrap=False)
    table.add_column("Type", style="magenta")
    table.add_column("Entries", justify="right", style="blue")
    table.add_column("Total Size", justify="right", style="green")
    table.add_column("File Size", justify="right", style="yellow")
    table.add_column("% of Total", justify="right", style="bright_magenta")

    if show_avg_size:
        table.add_column("Avg/Entry", justify="right", style="bright_cyan")

    for branch in branches:
        percentage = (
            (branch.file_size / total_file_size * 100) if total_file_size > 0 else 0.0
        )
        row = [
            branch.name,
            (
                branch.type_info[:50] + "..."
                if len(branch.type_info) > 50
                else branch.type_info
            ),
            f"{branch.entries:,}",
            human_readable_size(branch.total_size),
            human_readable_size(branch.file_size),
            f"{percentage:.1f}%",
        ]

        if show_avg_size:
            row.append(human_readable_size_float(branch.avg_size_per_entry))

        table.add_row(*row)

    # Add totals row
    total_size = sum(b.total_size for b in branches)
    total_entries = sum(b.entries for b in branches)

    table.add_section()
    total_row = [
        f"TOTAL ({len(branches)} branches)",
        "",
        f"{total_entries:,}",
        human_readable_size(total_size),
        human_readable_size(total_file_size),
        "100.0%",
    ]

    if show_avg_size:
        avg_total = sum(b.avg_size_per_entry for b in branches)
        total_row.append(human_readable_size_float(avg_total))

    table.add_row(*total_row, style="bold")

    console = Console()
    console.print(table)


app = typer.Typer(
    help="Display ROOT TTree branch sizes in a sorted, human-readable table"
)


@app.command()
def size(
    file: Annotated[str, typer.Argument(help="ROOT file path")],
    tree: Annotated[str, typer.Argument(help="TTree name")],
    sort: Annotated[SortBy, typer.Option(help="Sort by")] = SortBy.file_size,
    avg_size: Annotated[
        bool, typer.Option("--avg-size", help="Compute average size per entry (slow)")
    ] = False,
):
    """Display ROOT TTree branch sizes in a sorted, human-readable table."""
    # Open ROOT file and get tree
    root_file = ROOT.TFile.Open(file)
    if not root_file or root_file.IsZombie():
        typer.echo(f"Error: Cannot open file {file}", err=True)
        raise typer.Exit(1)

    ttree = root_file.Get(tree)
    if not ttree:
        typer.echo(f"Error: Cannot find tree {tree} in {file}", err=True)
        raise typer.Exit(1)

    # Parse and display
    branches = parse_tree_print(ttree)

    if avg_size:
        typer.echo("Computing average sizes per entry (this may take a while)...")
        compute_average_sizes(ttree, branches)

    display_branch_table(branches, sort_by=sort.value, show_avg_size=avg_size)

    root_file.Close()


@app.command()
def entries(
    files: Annotated[list[str], typer.Argument(help="ROOT file path")],
    markdown: bool = False,
):
    logger = colliderml.logging.get_logger(__name__)

    # figure out tree if not given
    first_file = ROOT.TFile.Open(files[0])
    keys = [k.GetName() for k in first_file.GetListOfKeys()]
    logger.info("Found keys: %s", keys)
    first_file.Close()

    table = Table(title="ROOT File Entries")
    table.add_column("File", style="yellow", no_wrap=False)
    for key in keys:
        table.add_column(key, style="cyan", no_wrap=False)
    table.add_column("file size", style="red", no_wrap=False)
    table.add_column("size / entry", style="red", no_wrap=False)

    column_names = ["File"] + keys + ["file size", "size / entry"]

    rows = []
    for file in files:
        # Open ROOT file and get tree
        root_file = ROOT.TFile.Open(file)
        if not root_file or root_file.IsZombie():
            typer.echo(f"Error: Cannot open file {file}", err=True)
            raise typer.Exit(1)

        row = []
        for tree in keys:
            ttree = root_file.Get(tree)
            if not ttree:
                typer.echo(f"Error: Cannot find tree {tree} in {file}", err=True)
                raise typer.Exit(1)

            row.append(ttree.GetEntries())

        max_entries = max(row)
        max_tree = keys[numpy.argmax(row)]
        file_size = os.path.getsize(file)
        size_per_entry = file_size / max_entries if max_entries > 0 else 0
        row = [f"{c:,}" for c in row]
        row.append(human_readable_size(file_size))
        row.append(f"{human_readable_size(size_per_entry)} ({max_tree})")

        table.add_row(file, *row)
        rows.append([file, *row])

        root_file.Close()

    if markdown:
        print(tabulate(rows, headers=column_names, tablefmt="github"))
        return
    console = Console()
    console.print(table)
