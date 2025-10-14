import dataclasses
import subprocess
from pathlib import Path
from typing import Generator, Self
import re
import json
import sys
import pydantic


def stream_subprocess(
    args: list[str],
    cwd: Path | str | None = None,
    check: bool = True,
    **kwargs,
) -> Generator[str, None, None]:
    """
    Run a command and yield output line by line as it's produced.

    Stderr is merged into stdout. Lines are yielded in real-time.
    Raises CalledProcessError if process fails and check=True.

    Args:
        args: Command and arguments to run
        cwd: Working directory for command
        check: If True, raise CalledProcessError on non-zero exit
        **kwargs: Additional arguments passed to subprocess.Popen

    Yields:
        Tuple of ('stdout', line) for each output line

    Raises:
        CalledProcessError: If check=True and process returns non-zero exit code

    Example:
        for stream, line in stream_subprocess(["ls", "-la"]):
            if stream == "stdout":
                print(line, end="")
    """
    # Remove conflicting kwargs that we handle ourselves
    kwargs.pop("stdout", None)
    kwargs.pop("stderr", None)
    kwargs.pop("text", None)
    kwargs.pop("capture_output", None)

    # If stdin not specified, inherit from parent for interactive prompts
    if "stdin" not in kwargs:
        kwargs["stdin"] = None

    process = subprocess.Popen(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # Merge stderr into stdout
        text=True,
        cwd=cwd,
        **kwargs,
    )

    # Collect output for potential error reporting
    output_lines = []

    # Stream output line by line
    if process.stdout:
        for line in process.stdout:
            output_lines.append(line)
            yield line

    # Wait for process to complete
    process.wait()

    # Check return code and raise if needed
    if check and process.returncode != 0:
        raise subprocess.CalledProcessError(
            process.returncode,
            args,
            output="".join(output_lines),
        )


def which(name: str) -> Path:
    exe = Path(
        subprocess.run(
            ["which", name],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    )

    if not exe.exists():
        raise FileNotFoundError(f"Executable '{name}' not found in PATH")

    return exe


@dataclasses.dataclass
class Hepmc3File:
    prefix: str
    compression: str

    @property
    def full_name(self) -> str:
        return f"{self.prefix}.hepmc3{self.compression_suffix}"

    @property
    def compression_suffix(self) -> str:
        suffix = ""
        if self.compression == "zlib":
            suffix = ".gz"
        elif self.compression == "lzma":
            suffix = ".xz"
        return suffix


def parse_hepmc3_file(file: Path):
    m = re.match(r"(.*)\.hepmc3(\.\w+)?", file.name)
    if m is None:
        raise ValueError("Unable to determine output prefix from output name")

    compression = "none"
    comp_suffix = m.group(2)
    if comp_suffix is not None:
        if comp_suffix == ".gz":
            compression = "zlib"
        elif comp_suffix == ".xz":
            compression = "lzma"
        elif comp_suffix == ".zst":
            compression = "zstd"
        else:
            raise ValueError(f"Unsupported compression suffix {comp_suffix}")

    prefix = m.group(1)

    return Hepmc3File(prefix=prefix, compression=compression)


def hepmc3_get_meta_path(file: Path) -> Path:
    return file.parent / f"{file.name}.json"


class HepMC3Meta(pydantic.BaseModel):
    num_events: int

    @classmethod
    def for_file(cls, file: Path) -> Self:
        with hepmc3_get_meta_path(file).open("r") as meta_f:
            meta = json.load(meta_f)
            return cls.model_validate(meta)


def hadd(input_files: list[Path], output_file: Path):
    hadd_exe = which("hadd")

    subprocess.run(
        [hadd_exe, str(output_file)] + [str(f) for f in input_files], check=True
    )


def human_readable_size(size_bytes: int) -> str:
    """Convert bytes to human-readable format."""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"
