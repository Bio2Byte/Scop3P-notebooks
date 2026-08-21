"""An input must be updated with the function that matches its widget.

``ui.update_select`` against a selectize widget, or ``ui.update_selectize`` against a
plain select, is accepted by Shiny and then does nothing. There is no error and no log
line -- the control simply keeps its old options, which looks exactly like an upstream
returning no results. Converting the structure pickers to searchable selects created 7
new chances to get this wrong, so it is checked against the source rather than trusted.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

APPS_DIR = Path(__file__).resolve().parents[2] / "apps"
APP_FILES = sorted(APPS_DIR.glob("*/app.py"))

#: Helpers that create a widget of a given kind, beyond the ui.* constructors.
SELECTIZE_FACTORIES = {"scop3p_structure_picker"}


def _calls(tree: ast.AST) -> list[ast.Call]:
    return [node for node in ast.walk(tree) if isinstance(node, ast.Call)]


def _callee(node: ast.Call) -> str:
    if isinstance(node.func, ast.Attribute):
        return node.func.attr
    if isinstance(node.func, ast.Name):
        return node.func.id
    return ""


def _first_string_arg(node: ast.Call) -> str | None:
    if node.args and isinstance(node.args[0], ast.Constant):
        value = node.args[0].value
        if isinstance(value, str):
            return value
    return None


def _widget_kinds(tree: ast.AST) -> dict[str, str]:
    """input id -> "select" | "selectize", from how the widget was created."""
    kinds: dict[str, str] = {}
    for call in _calls(tree):
        name = _callee(call)
        input_id = _first_string_arg(call)
        if input_id is None:
            continue
        if name in SELECTIZE_FACTORIES or name == "input_selectize":
            kinds[input_id] = "selectize"
        elif name == "input_select":
            # input_select(selectize=True) is a selectize widget wearing another name.
            selectize = any(
                kw.arg == "selectize"
                and isinstance(kw.value, ast.Constant)
                and kw.value.value is True
                for kw in call.keywords
            )
            kinds[input_id] = "selectize" if selectize else "select"
    return kinds


def _update_calls(tree: ast.AST) -> list[tuple[str, str, int]]:
    """(input id, "select" | "selectize", line) for every update call with a literal id."""
    found: list[tuple[str, str, int]] = []
    for call in _calls(tree):
        name = _callee(call)
        if name not in {"update_select", "update_selectize"}:
            continue
        input_id = _first_string_arg(call)
        if input_id is None:
            continue  # a loop variable; covered by the dedicated test below
        found.append((input_id, name.removeprefix("update_"), call.lineno))
    return found


def test_there_are_apps_to_check() -> None:
    assert len(APP_FILES) >= 5


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_updates_match_the_widget_they_target(path: Path) -> None:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    kinds = _widget_kinds(tree)

    mismatches = [
        f"{path.parent.name}:{line} update_{used}(\"{input_id}\") "
        f"but it was created as a {kinds[input_id]}"
        for input_id, used, line in _update_calls(tree)
        if input_id in kinds and kinds[input_id] != used
    ]
    assert not mismatches, "; ".join(mismatches)


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_a_loop_over_input_ids_updates_one_widget_kind(path: Path) -> None:
    """structure_viz updates four pickers in a loop, so they must agree on their kind.

    If one of them were left a plain select, the single update call in the loop could not
    be right for all four.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    kinds = _widget_kinds(tree)

    for node in ast.walk(tree):
        if not isinstance(node, ast.For) or not isinstance(node.iter, ast.Tuple):
            continue
        looped = [
            element.value
            for element in node.iter.elts
            if isinstance(element, ast.Constant) and isinstance(element.value, str)
        ]
        looped = [input_id for input_id in looped if input_id in kinds]
        if len(looped) < 2:
            continue
        distinct = {kinds[input_id] for input_id in looped}
        assert len(distinct) == 1, (
            f"{path.parent.name}:{node.lineno} loops over {looped} but they are "
            f"{distinct} -- one update call cannot serve both"
        )


def test_the_searchable_pickers_are_actually_searchable() -> None:
    """The point of the conversion: P04637 offers 629 chains in RIN Alignment.

    Pinned by id, because reverting any one of these to a plain select would restore an
    unscannable dropdown without failing anything else.
    """
    expected = {
        "structure_viz": {"pdb_id", "rin_pdb_id", "tm_pdb1_id", "tm_pdb2_id"},
        "rinalign": {"left_structure", "right_structure"},
        "topology_viewer": {"structure_choice"},
    }
    for app_name, ids in expected.items():
        path = APPS_DIR / app_name / "app.py"
        kinds = _widget_kinds(ast.parse(path.read_text(encoding="utf-8")))
        for input_id in ids:
            assert kinds.get(input_id) == "selectize", f"{app_name}.{input_id} is not searchable"


def test_small_pickers_stay_plain_selects() -> None:
    """A search box for two chains is more clicks, not fewer."""
    kinds = _widget_kinds(
        ast.parse((APPS_DIR / "structure_viz" / "app.py").read_text(encoding="utf-8"))
    )
    for input_id in ("chain", "rin_chain", "rin_color_mode", "rin_b2b_metric"):
        assert kinds.get(input_id) == "select", f"{input_id} became searchable needlessly"
