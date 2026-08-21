"""Repo-wide logging conventions, enforced structurally.

These are the rules that keep the log readable across five protocols. Reviewing them by
eye does not scale to 150-odd call sites, and a single stray print or a mis-levelled
line is invisible in a passing test suite -- so they are asserted against the source.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

APPS_DIR = Path(__file__).resolve().parents[2] / "apps"
PY_FILES = sorted(path for path in APPS_DIR.rglob("*.py") if "__pycache__" not in str(path))
APP_FILES = sorted(APPS_DIR.glob("*/app.py"))


def test_there_are_files_to_check() -> None:
    assert len(PY_FILES) >= 10
    assert len(APP_FILES) >= 5


@pytest.mark.parametrize("path", PY_FILES, ids=lambda p: str(p.relative_to(APPS_DIR)))
def test_no_print_statements(path: Path) -> None:
    """Output goes through logging, so it carries a level, a timestamp and a name.

    A print also cannot be silenced by SCOP3P_LOG_LEVEL and never reaches the log file,
    which is how debug output ends up in a user's console.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    offenders = [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "print"
    ]
    assert not offenders, f"print() at {path.name}:{offenders} -- use the logger instead"


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_every_app_opens_an_experiment_trail(path: Path) -> None:
    """Without this the run has no record of which protocol was opened."""
    source = path.read_text(encoding="utf-8")
    assert "new_trail()" in source, f"{path.parent.name} never creates a trail"
    assert "trail.opened(" in source, f"{path.parent.name} never records the protocol opening"


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_the_trail_is_per_session_not_per_process(path: Path) -> None:
    """A module-level trail would interleave every user's steps into one sequence."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in tree.body:  # module level only
        for child in ast.walk(node):
            if (
                isinstance(child, ast.Call)
                and isinstance(child.func, ast.Name)
                and child.func.id == "new_trail"
            ):
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    continue
                raise AssertionError(
                    f"{path.parent.name}: new_trail() at module scope (line {child.lineno})"
                )


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_no_app_still_uses_the_retired_click_helper(path: Path) -> None:
    """One way to log a click, so the vocabulary cannot drift back apart."""
    assert "log_action_button_click" not in path.read_text(encoding="utf-8")


def test_b2btools_is_always_called_inside_the_quiet_wrapper() -> None:
    """b2bTools prints its progress and its dependencies warn on every prediction.

    Both must be captured to DEBUG. This is asserted rather than remembered because the
    symptom -- a console full of third-party chatter -- shows up for the user, not in CI.
    """
    predict_sites = []
    for path in PY_FILES:
        source = path.read_text(encoding="utf-8")
        if "SingleSeq(" not in source:
            continue
        tree = ast.parse(source)
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            body_source = ast.get_source_segment(source, node) or ""
            if "SingleSeq(" not in body_source:
                continue
            predict_sites.append((path.name, node.name))
            assert "quiet_third_party" in body_source, (
                f"{path.name}:{node.name} calls b2bTools without quiet_third_party"
            )
    assert predict_sites, "no b2bTools call sites found -- this test would be vacuous"


@pytest.mark.parametrize(
    "path",
    sorted(p for p in (APPS_DIR.parent / "tests").rglob("*.py")),
    ids=lambda p: p.name,
)
def test_tests_import_project_code_under_one_name(path: Path) -> None:
    """``common.x`` and ``apps.common.x`` are two different modules.

    pytest.ini puts both ``.`` and ``apps`` on the path, so the same file can be imported
    twice under two names -- giving two module objects, two loggers, and two copies of
    every module-level object. That broke cache isolation in a way that looked like a
    flaky test: the fixture cleared one registry while the service used the other, so a
    test passed alone and failed in its file.
    """
    if path.name == Path(__file__).name:
        pytest.skip("this module names the pattern in order to check for it")
    source = path.read_text(encoding="utf-8")
    assert "apps.common." not in source, (
        f"{path.name} refers to apps.common.*; use common.* so there is exactly one "
        "module identity for project code"
    )
