from __future__ import annotations

from pathlib import Path

import pytest

from common import topology_bridge as bridge


def _make_package(root: Path, name: str = "topology") -> Path:
    package = root / name
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("VALUE = 1\n")
    return root


def test_candidate_order_puts_env_override_first(tmp_path: Path) -> None:
    candidates = bridge.candidate_dirs({bridge.ENV_VAR: str(tmp_path / "override")})
    assert candidates[0] == tmp_path / "override"
    assert candidates[1] == bridge.CONTAINER_DIR
    assert candidates[2].name == "topology_viewer"


def test_candidate_order_without_override() -> None:
    candidates = bridge.candidate_dirs({})
    assert candidates[0] == bridge.CONTAINER_DIR
    assert candidates[1].parts[-2:] == ("notebooks", "topology_viewer")


def test_blank_override_is_ignored() -> None:
    assert bridge.candidate_dirs({bridge.ENV_VAR: "   "})[0] == bridge.CONTAINER_DIR


def test_resolve_prefers_the_earlier_candidate(tmp_path: Path) -> None:
    first = _make_package(tmp_path / "first")
    second = _make_package(tmp_path / "second")
    assert bridge.resolve_topology_dir([first, second]) == first
    assert bridge.resolve_topology_dir([second, first]) == second


def test_resolve_skips_a_directory_without_the_package(tmp_path: Path) -> None:
    """A stale, empty container directory must not shadow a real checkout.

    This is why the existence test targets ``<dir>/topology/__init__.py`` rather
    than ``<dir>`` itself.
    """
    empty = tmp_path / "empty"
    (empty / "topology").mkdir(parents=True)  # directory present, no __init__.py
    real = _make_package(tmp_path / "real")
    assert bridge.resolve_topology_dir([empty, real]) == real


def test_resolve_returns_none_when_nothing_matches(tmp_path: Path) -> None:
    assert bridge.resolve_topology_dir([tmp_path / "nope"]) is None
    assert bridge.resolve_topology_dir([]) is None


def test_failure_message_names_every_candidate_and_the_env_var(tmp_path: Path) -> None:
    candidates = [tmp_path / "a", bridge.CONTAINER_DIR, tmp_path / "c"]
    message = bridge.describe_failure(candidates)
    for candidate in candidates:
        assert str(candidate) in message
    assert bridge.ENV_VAR in message
    assert "topology/__init__.py" in message


def test_unavailable_is_not_an_import_error() -> None:
    """The app degrades on TopologyUnavailable; an ImportError would be a traceback."""
    assert issubclass(bridge.TopologyUnavailable, RuntimeError)
    assert not issubclass(bridge.TopologyUnavailable, ImportError)


def test_bridge_resolves_in_a_source_checkout() -> None:
    """Guards the ``parents[2]`` arithmetic and the notebooks/ layout.

    If the topology package is ever moved, or this file is relocated within
    ``apps/common/``, this is the test that says so.
    """
    assert bridge.TOPOLOGY_ERROR is None, bridge.TOPOLOGY_ERROR
    assert bridge.topology is not None
    assert bridge.TOPOLOGY_DIR is not None
    assert (bridge.TOPOLOGY_DIR / "topology" / "__init__.py").is_file()
    assert bridge.__topology_version__


@pytest.mark.parametrize(
    "name",
    ["build_view", "fetch_alphafold", "save_html", "load_structure", "Structure"],
)
def test_reexported_callables_are_present(name: str) -> None:
    assert getattr(bridge, name) is not None


def test_fetch_alphafold_is_reached_through_the_app_submodule() -> None:
    """It is not re-exported by ``topology/__init__.py``, hence the dotted lookup."""
    assert not hasattr(bridge.topology, "fetch_alphafold")
    assert bridge.fetch_alphafold is bridge.topology.app.fetch_alphafold


@pytest.mark.parametrize(
    "name",
    [
        "fetch_structures",
        "fetch_ptms",
        "fetch_variants",
        "fetch_numbering",
        "fetch_structure_file",
        "structure_file_url",
    ],
)
def test_annotation_helpers_the_app_calls_exist(name: str) -> None:
    assert hasattr(bridge.annotations_module, name)
