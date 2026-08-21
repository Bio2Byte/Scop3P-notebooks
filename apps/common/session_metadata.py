from __future__ import annotations

from importlib import metadata
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from typing import Any


_DEPENDENCIES = (
    "shiny",
    "pandas",
    "requests",
    "numpy",
    "scipy",
    "networkx",
    "biopython",
    "py3Dmol",
    "pyvis",
    "bokeh",
    "b2bTools",
)
_TOOL_PROBES = {
    "TM-align": ("TM-align",),
    "hmmer": ("hmmsearch", "-h"),
    "t_coffee": ("t_coffee", "-version"),
}


def write_metadata(
    *,
    metadata_path: Path,
    log_file_path: Path,
    log_dir: Path,
    session_started_at: str,
    trail_file_path: Path | None = None,
) -> Path:
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    payload = build_metadata(
        log_file_path=log_file_path,
        log_dir=log_dir,
        session_started_at=session_started_at,
        trail_file_path=trail_file_path,
    )
    metadata_path.write_text(_to_yaml(payload), encoding="utf-8")
    return metadata_path


def build_metadata(
    *,
    log_file_path: Path,
    log_dir: Path,
    session_started_at: str,
    trail_file_path: Path | None = None,
) -> dict[str, Any]:
    return {
        "schema_version": "1.0",
        "application": {
            "name": os.getenv("SCOP3P_APP_NAME", "unknown"),
            "title": "Scop3P-Toolkit",
            "description": "Tools for exploring and extending Scop3P",
        },
        "session": {
            "started_at_utc": session_started_at,
            "working_directory": str(Path.cwd()),
            "log_directory": str(log_dir),
            "log_file": str(log_file_path),
            # The step-by-step experiment record, recorded separately so a run can be
            # handed over as a standalone document.
            "trail_file": str(trail_file_path) if trail_file_path else None,
        },
        "image": {
            "title": os.getenv("SCOP3P_IMAGE_TITLE", "Scop3P-Toolkit"),
            "version": os.getenv("SCOP3P_IMAGE_VERSION", "unknown"),
            "revision": os.getenv("SCOP3P_IMAGE_REVISION", "unknown"),
            "created": os.getenv("SCOP3P_IMAGE_CREATED", "unknown"),
            "source": os.getenv("SCOP3P_IMAGE_SOURCE", "https://github.com/Bio2Byte/Scop3P-notebooks"),
            "license": os.getenv("SCOP3P_IMAGE_LICENSE", "Apache-2.0"),
        },
        "runtime": {
            "python": sys.version.split()[0],
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "machine": platform.machine(),
        },
        "dependencies": _dependency_versions(),
        "external_tools": _external_tool_versions(),
    }


def _dependency_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for package in _DEPENDENCIES:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def _external_tool_versions() -> dict[str, dict[str, str | None]]:
    return {name: _probe_tool(command) for name, command in _TOOL_PROBES.items()}


def _probe_tool(command: tuple[str, ...]) -> dict[str, str | None]:
    executable = shutil.which(command[0])
    if executable is None:
        return {"path": None, "version": None}

    try:
        process = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=2,
        )
    except (OSError, subprocess.TimeoutExpired):
        return {"path": executable, "version": None}

    output = "\n".join(part.strip() for part in (process.stdout, process.stderr) if part.strip())
    first_line = output.splitlines()[0] if output else None
    return {"path": executable, "version": first_line if process.returncode == 0 else None}


def _to_yaml(value: Any, *, indent: int = 0) -> str:
    lines = list(_yaml_lines(value, indent=indent))
    return "\n".join(lines) + "\n"


def _yaml_lines(value: Any, *, indent: int) -> list[str]:
    prefix = " " * indent
    if isinstance(value, dict):
        lines: list[str] = []
        for key, item in value.items():
            if isinstance(item, dict):
                lines.append(f"{prefix}{key}:")
                lines.extend(_yaml_lines(item, indent=indent + 2))
            elif isinstance(item, list):
                lines.append(f"{prefix}{key}:")
                lines.extend(_yaml_lines(item, indent=indent + 2))
            else:
                lines.append(f"{prefix}{key}: {_yaml_scalar(item)}")
        return lines
    if isinstance(value, list):
        lines = []
        for item in value:
            if isinstance(item, dict):
                lines.append(f"{prefix}-")
                lines.extend(_yaml_lines(item, indent=indent + 2))
            else:
                lines.append(f"{prefix}- {_yaml_scalar(item)}")
        return lines
    return [f"{prefix}{_yaml_scalar(value)}"]


def _yaml_scalar(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int | float):
        return str(value)
    text = str(value)
    escaped = text.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'
