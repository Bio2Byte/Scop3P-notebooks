"""Structural guard against Shiny's silently-suspended outputs.

Shiny decides whether to suspend an output from the client-reported
``.clientdata_output_<id>_hidden`` value, and ``Session._is_hidden()`` treats "never
reported" as hidden. An output that lives in a ``nav_panel`` which is not the initially
active tab is therefore suspended at page load and never woken when the user opens its
tab: it sits at "recalculating" forever, with nothing logged anywhere. The fix is to
declare it ``@output(suspend_when_hidden=False)``.

That failure is invisible in a smoke test -- the app constructs, serves, and returns 200 --
so it has to be caught structurally. This walks each app's AST rather than trusting anyone
to remember the decorator.

Verified against shiny 1.7.0, which requirements-shiny.txt permits (shiny>=1.1,<2).
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

APPS_DIR = Path(__file__).resolve().parents[2] / "apps"
APP_FILES = sorted(APPS_DIR.glob("*/app.py"))


def _is_call_to(node: ast.AST, prefix: str) -> bool:
    """True for ``ui.<prefix>...(...)`` calls."""
    return (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr.startswith(prefix)
    )


def _output_ids(node: ast.AST) -> set[str]:
    """Every ``ui.output_*("id")`` id anywhere below this node."""
    found: set[str] = set()
    for child in ast.walk(node):
        if _is_call_to(child, "output_") and child.args:
            first = child.args[0]
            if isinstance(first, ast.Constant) and isinstance(first.value, str):
                found.add(first.value)
    return found


def _ids_in_non_initial_panels(tree: ast.Module) -> set[str]:
    """Output ids that render inside a tab which is not the one shown at page load.

    No navset in this repo pins ``selected=``, so the first ``nav_panel`` is the active
    one and every later panel starts hidden. A navset that grows a ``selected=`` kwarg
    would make this reasoning wrong, so that is asserted rather than assumed.
    """
    hidden: set[str] = set()
    for node in ast.walk(tree):
        if not _is_call_to(node, "navset"):
            continue
        assert not any(
            keyword.arg == "selected" for keyword in node.keywords
        ), "a navset pins selected= -- this test's 'first panel is active' rule no longer holds"
        panels = [arg for arg in node.args if _is_call_to(arg, "nav_panel")]
        for panel in panels[1:]:
            hidden |= _output_ids(panel)
    return hidden


def _render_functions(tree: ast.Module) -> dict[str, bool]:
    """Every ``@render.*`` function, mapped to whether it carries the guard."""
    functions: dict[str, bool] = {}
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        renders = guarded = False
        for decorator in node.decorator_list:
            if isinstance(decorator, ast.Attribute) and _name_of(decorator.value) == "render":
                renders = True
            elif isinstance(decorator, ast.Call) and _name_of(decorator.func) == "output":
                guarded = any(
                    keyword.arg == "suspend_when_hidden"
                    and isinstance(keyword.value, ast.Constant)
                    and keyword.value.value is False
                    for keyword in decorator.keywords
                )
        if renders:
            functions[node.name] = guarded
    return functions


def _name_of(node: ast.AST) -> str | None:
    return node.id if isinstance(node, ast.Name) else None


def test_there_are_apps_to_check() -> None:
    """Guards against the glob silently matching nothing."""
    assert len(APP_FILES) >= 5


@pytest.mark.parametrize("app_file", APP_FILES, ids=lambda path: path.parent.name)
def test_outputs_on_hidden_tabs_disable_suspension(app_file: Path) -> None:
    tree = ast.parse(app_file.read_text(encoding="utf-8"))
    hidden_ids = _ids_in_non_initial_panels(tree)
    renders = _render_functions(tree)

    unguarded = sorted(
        output_id
        for output_id in hidden_ids
        if output_id in renders and not renders[output_id]
    )
    assert not unguarded, (
        f"{app_file.parent.name}: {unguarded} render on a tab that is not active at page "
        "load, so Shiny suspends them and never wakes them. Add "
        "@output(suspend_when_hidden=False) above @render.*"
    )


@pytest.mark.parametrize("app_file", APP_FILES, ids=lambda path: path.parent.name)
def test_the_hidden_tab_scan_finds_something_to_check(app_file: Path) -> None:
    """A positive control: an app with tabs must yield ids on non-initial tabs.

    Without this, a refactor that stopped matching ``nav_panel`` would make the test above
    pass by finding nothing at all.
    """
    tree = ast.parse(app_file.read_text(encoding="utf-8"))
    # Test for a real navset call, not the substring: several apps mention "navset" only
    # in a comment explaining why they deliberately avoid one.
    if not any(_is_call_to(node, "navset") for node in ast.walk(tree)):
        pytest.skip("single-page app, no tabs to hide")
    assert _ids_in_non_initial_panels(tree), "tabs present but no output ids found on them"
