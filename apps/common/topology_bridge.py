"""Locate and import the topology viewer package.

The package lives at ``notebooks/topology_viewer/topology/`` and stays there: it is
the single source of truth shared by the Voila notebook, its standalone test script
(``notebooks/topology_viewer/test_topology.py``) and the Shiny app in
``apps/topology_viewer/``. Because ``docker/Dockerfile`` only copies ``./apps`` into
the image, the package is copied separately to ``/opt/scop3p/topology_viewer`` and
this module finds it in whichever environment it is running:

    1. ``$SCOP3P_TOPOLOGY_PATH``      operator override
    2. ``/opt/scop3p/topology_viewer`` the in-container copy
    3. ``<repo>/notebooks/topology_viewer``  a source checkout: ``shiny run``, pytest

Two deliberate choices worth keeping:

``sys.path.append``, never ``insert(0)``. ``notebooks/topology_viewer/`` also holds
``preview.py``, ``make_fixtures.py``, ``test_topology.py`` and ``dssp_topology_app.py``.
Appending makes them importable as top-level names but loses every name race against
real repository modules, which is the outcome we want.

**This module never raises on import.** ``apps/portal/main.py`` imports every app at
module scope, so an exception here would take down the other four apps and fail four
unrelated smoke tests. Callers read :data:`TOPOLOGY_ERROR` and degrade; only the
explicit :func:`load_topology` call raises, and it raises :class:`TopologyUnavailable`
with the full candidate list rather than a bare ``ImportError`` traceback.
"""

from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Iterable, Mapping

TOPOLOGY_PACKAGE = "topology"
ENV_VAR = "SCOP3P_TOPOLOGY_PATH"
CONTAINER_DIR = Path("/opt/scop3p/topology_viewer")
_REPO_DIR = Path(__file__).resolve().parents[2] / "notebooks" / "topology_viewer"


class TopologyUnavailable(RuntimeError):
    """The topology package could not be found on any known path."""


def candidate_dirs(env: Mapping[str, str] | None = None) -> list[Path]:
    """Directories that may contain the ``topology`` package, in priority order.

    Pure: touches no filesystem, so the precedence rules are testable on their own.
    """
    environment = os.environ if env is None else env
    candidates: list[Path] = []
    override = (environment.get(ENV_VAR) or "").strip()
    if override:
        candidates.append(Path(override))
    candidates.append(CONTAINER_DIR)
    candidates.append(_REPO_DIR)
    return candidates


def resolve_topology_dir(candidates: Iterable[Path]) -> Path | None:
    """First candidate that actually holds the package, or ``None``.

    The test is for ``<candidate>/topology/__init__.py`` rather than for the
    candidate directory itself. An empty ``/opt/scop3p/topology_viewer`` left behind
    by a partial image build would otherwise shadow a working source checkout.
    """
    for candidate in candidates:
        try:
            if (candidate / TOPOLOGY_PACKAGE / "__init__.py").is_file():
                return candidate
        except OSError:
            # An unreadable or malformed path is just not a candidate.
            continue
    return None


def describe_failure(candidates: Iterable[Path]) -> str:
    """The message a user gets when nothing resolved. Names every path tried."""
    tried = " ".join(
        f"({index}) {candidate}" for index, candidate in enumerate(candidates, start=1)
    )
    return (
        f"Topology package not found. Looked for '<dir>/{TOPOLOGY_PACKAGE}/__init__.py' "
        f"in: {tried or '(no candidates)'}. Set {ENV_VAR} to the directory that "
        f"contains the '{TOPOLOGY_PACKAGE}' package."
    )


_CACHED: ModuleType | None = None


def load_topology() -> ModuleType:
    """Import the topology package, memoised.

    Raises :class:`TopologyUnavailable` when no candidate holds the package.
    """
    global _CACHED
    if _CACHED is not None:
        return _CACHED

    candidates = candidate_dirs()
    resolved = resolve_topology_dir(candidates)
    if resolved is None:
        raise TopologyUnavailable(describe_failure(candidates))

    if str(resolved) not in sys.path:
        sys.path.append(str(resolved))

    try:
        module = importlib.import_module(TOPOLOGY_PACKAGE)
    except Exception as error:  # noqa: BLE001 - reported, not swallowed
        raise TopologyUnavailable(
            f"Found the topology package at {resolved} but importing it failed: "
            f"{type(error).__name__}: {error}"
        ) from error

    loaded_from = Path(getattr(module, "__file__", "") or "").resolve().parent.parent
    if loaded_from != resolved:
        # Something else called 'topology' was already importable. Say so rather
        # than silently drawing diagrams with code nobody meant to ship.
        raise TopologyUnavailable(
            f"Resolved the topology package to {resolved} but 'import topology' "
            f"loaded {loaded_from} instead. Remove the conflicting package from "
            f"sys.path, or set {ENV_VAR} explicitly."
        )

    _CACHED = module
    return module


topology: ModuleType | None
TOPOLOGY_ERROR: str | None
TOPOLOGY_DIR: Path | None

try:
    topology = load_topology()
    TOPOLOGY_ERROR = None
    TOPOLOGY_DIR = resolve_topology_dir(candidate_dirs())
except TopologyUnavailable as error:
    topology = None
    TOPOLOGY_ERROR = str(error)
    TOPOLOGY_DIR = None


def _attr(dotted: str):
    """Resolve a possibly dotted attribute path, or ``None``.

    ``fetch_alphafold`` is not re-exported by ``topology/__init__.py``, so it has
    to be reached through the ``topology.app`` submodule.
    """
    if topology is None:
        return None
    current: object = topology
    for part in dotted.split("."):
        current = getattr(current, part, None)
        if current is None:
            return None
    return current


# The names the Shiny app needs, and nothing more. All are ``None`` when the
# package is unavailable; the app checks TOPOLOGY_ERROR before touching them.
build_view = _attr("build_view")
fetch_alphafold = _attr("app.fetch_alphafold")
save_html = _attr("save_html")
load_structure = _attr("load_structure")
Structure = _attr("Structure")
annotations_module = _attr("annotations")
logo = _attr("logo")
__topology_version__ = _attr("__version__")

__all__ = [
    "TOPOLOGY_PACKAGE",
    "ENV_VAR",
    "CONTAINER_DIR",
    "TopologyUnavailable",
    "candidate_dirs",
    "resolve_topology_dir",
    "describe_failure",
    "load_topology",
    "topology",
    "TOPOLOGY_ERROR",
    "TOPOLOGY_DIR",
    "build_view",
    "fetch_alphafold",
    "save_html",
    "load_structure",
    "Structure",
    "annotations_module",
    "logo",
    "__topology_version__",
]
