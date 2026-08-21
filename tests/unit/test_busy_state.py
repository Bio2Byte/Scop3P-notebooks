"""Waiting-time handling: the contracts that a passing app can still violate.

A task button disables itself in the browser the moment it is clicked. That is what stops
a second click queuing the same slow work again -- and it is also the hazard: if the server
never releases the button, the control stays dead for the rest of the session with nothing
on screen to explain why. Nothing raises, no test fails, and the app looks fine until you
click the button.

So the pairing is checked against the source: every ``background_task_button`` must have a
matching ``finish_task``, and no ``task_button`` may have one.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

APPS_DIR = Path(__file__).resolve().parents[2] / "apps"
APP_FILES = sorted(APPS_DIR.glob("*/app.py"))

#: Buttons whose handler is instant and local. A spinner that flashes for 20ms is noise.
INSTANT_BUTTONS = {"load_example", "load_example_mut", "clear_annotations", "reset_b2b"}


def _calls(tree: ast.AST, name: str) -> list[ast.Call]:
    return [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and (
            (isinstance(node.func, ast.Name) and node.func.id == name)
            or (isinstance(node.func, ast.Attribute) and node.func.attr == name)
        )
    ]


def _first_str(node: ast.Call) -> str | None:
    if node.args and isinstance(node.args[0], ast.Constant):
        value = node.args[0].value
        if isinstance(value, str):
            return value
    return None


def _ids(tree: ast.AST, name: str) -> set[str]:
    return {found for node in _calls(tree, name) if (found := _first_str(node))}


def test_there_are_apps_to_check() -> None:
    assert len(APP_FILES) >= 5


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_every_background_button_is_released(path: Path) -> None:
    """The stuck-button bug: disabled on click, never re-enabled."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    background = _ids(tree, "background_task_button")
    released = _ids(tree, "finish_task")
    missing = sorted(background - released)
    assert not missing, (
        f"{path.parent.name}: {missing} use background_task_button (auto_reset=False) but "
        "never call finish_task, so the button stays disabled forever"
    )


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_a_self_resetting_button_is_not_also_released_manually(path: Path) -> None:
    """``task_button`` resets itself; calling finish_task for one signals a mix-up."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    overlap = sorted(_ids(tree, "task_button") & _ids(tree, "finish_task"))
    assert not overlap, (
        f"{path.parent.name}: {overlap} are self-resetting task_buttons but also call "
        "finish_task -- they were probably meant to be background_task_buttons"
    )


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_finish_task_targets_a_real_button(path: Path) -> None:
    """A typo'd id in finish_task silently releases nothing."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    known = _ids(tree, "background_task_button") | _ids(tree, "task_button")
    if not known:
        pytest.skip("no task buttons in this app")
    unknown = sorted(_ids(tree, "finish_task") - known)
    assert not unknown, f"{path.parent.name}: finish_task({unknown}) matches no button"


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_a_background_task_exists_where_one_is_claimed(path: Path) -> None:
    """``background_task_button`` promises the work is off the event loop.

    Without a ``background(...)`` task the button would be held busy by a handler that is
    blocking the loop anyway, which is the worst of both: no live message *and* a control
    that needs manual release.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    if not _ids(tree, "background_task_button"):
        pytest.skip("no background buttons in this app")
    assert _calls(tree, "background"), (
        f"{path.parent.name}: claims a background button but never calls background()"
    )


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_instant_actions_are_not_dressed_as_slow_ones(path: Path) -> None:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    busy = _ids(tree, "task_button") | _ids(tree, "background_task_button")
    wrong = sorted(busy & INSTANT_BUTTONS)
    assert not wrong, f"{path.parent.name}: {wrong} are instant; a spinner there is noise"


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_an_app_with_busy_buttons_shows_spinners(path: Path) -> None:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    if not (_ids(tree, "task_button") | _ids(tree, "background_task_button")):
        pytest.skip("no task buttons in this app")
    assert _calls(tree, "busy_indicators"), (
        f"{path.parent.name}: has task buttons but never includes busy_indicators()"
    )


def test_the_slow_actions_across_every_protocol_are_guarded() -> None:
    """Pinned by id, so a slow action cannot quietly go back to an unguarded button.

    These are the handlers that make the user wait: a network round trip, a structure
    download, a KD-tree, a subprocess, or a Bio2Byte prediction.
    """
    expected = {
        "structure_viz": {
            "set_accession", "fetch_ptm", "fetch_variants", "fetch_af", "render_structure",
            "fetch_seq", "run_b2b", "render_b2b_3d", "rin_dl_af", "build_rin", "show_rin",
            "run_tmalign", "load_tmalign_structures",
        },
        "mutation_effect": {"run_wt", "run_mut", "run_inf"},
        "peptide_mapper": {"load_btn", "load_upload", "build_mapping", "map_all", "export_html"},
        "topology_viewer": {"fetch_btn", "fetch_ptms", "fetch_variants"},
        "rinalign": {"fetch", "generate", "compare"},
    }
    for app_name, ids in expected.items():
        tree = ast.parse((APPS_DIR / app_name / "app.py").read_text(encoding="utf-8"))
        guarded = _ids(tree, "task_button") | _ids(tree, "background_task_button")
        missing = sorted(ids - guarded)
        assert not missing, f"{app_name}: {missing} are unguarded slow actions"


def test_the_heaviest_work_runs_off_the_event_loop() -> None:
    """Bio2Byte is the longest wait in the suite, so it earns a live status message.

    A synchronous handler cannot show one: the message is queued until the handler
    returns. These two must therefore stay on the background path.
    """
    for app_name, input_id in (("structure_viz", "run_b2b"), ("mutation_effect", "run_wt")):
        tree = ast.parse((APPS_DIR / app_name / "app.py").read_text(encoding="utf-8"))
        assert input_id in _ids(tree, "background_task_button"), (
            f"{app_name}.{input_id} no longer runs off the event loop, so its "
            "'working...' message can never reach the browser"
        )


# --------------------------------------------------------------------------------------
# The helpers must accept what the call sites actually pass
# --------------------------------------------------------------------------------------


def _keywords_used(tree: ast.AST, name: str) -> set[str]:
    used: set[str] = set()
    for call in _calls(tree, name):
        used.update(keyword.arg for keyword in call.keywords if keyword.arg)
    return used


@pytest.mark.parametrize("helper", ["task_button", "background_task_button"])
@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_the_button_helpers_accept_every_keyword_passed_to_them(path: Path, helper: str) -> None:
    """A swallowed keyword raises at *render* time, not at construction.

    Several of these buttons are built inside a ``@render.ui`` output and pass ``disabled=``
    to reflect whether their prerequisite is met. When the helper did not accept it, the
    TypeError replaced the entire output with an error message and the button vanished --
    on the TM-align tab and on RIN Alignment's Compare / Align. Constructing the app UI
    does not evaluate those outputs, so every existing test still passed.
    """
    import inspect

    from common import busy

    signature = inspect.signature(getattr(busy, helper))
    accepts_extra = any(
        parameter.kind is inspect.Parameter.VAR_KEYWORD
        for parameter in signature.parameters.values()
    )
    tree = ast.parse(path.read_text(encoding="utf-8"))
    used = _keywords_used(tree, helper)
    if not used:
        pytest.skip(f"{path.parent.name} does not call {helper} with keywords")

    unsupported = sorted(
        keyword for keyword in used if keyword not in signature.parameters
    )
    assert accepts_extra or not unsupported, (
        f"{path.parent.name} passes {unsupported} to {helper}, which does not accept it; "
        "if the button is built in a @render.ui output this raises at render time and the "
        "control disappears"
    )


@pytest.mark.parametrize("helper", ["task_button", "background_task_button"])
def test_the_button_helpers_really_forward_a_disabled_flag(helper: str) -> None:
    """Accepting the keyword is not enough -- it has to reach the rendered markup.

    ``**kwargs`` that is accepted and then dropped would leave a button that looks enabled
    while its prerequisite is unmet, which is worse than the crash it replaced.
    """
    from common import busy

    build = getattr(busy, helper)
    assert "disabled" in str(build("b", "Go", disabled=True))
    assert "disabled" not in str(build("b", "Go", disabled=False))


def test_conditionally_enabled_buttons_still_carry_their_gate() -> None:
    """These two are only clickable once their prerequisite exists.

    They must be gated with ``gate(...)``, not ``disabled=``: a task button owns its own
    ``disabled`` property and drops the attribute on initialisation, so the server's flag
    is discarded and the button renders fully clickable. Verified in the browser -- the
    markup carried no ``disabled`` attribute at all. Losing the gate is silent: the button
    looks fine and simply lets the user start work that cannot succeed.
    """
    expected = {
        "structure_viz": "run_tmalign",
        "rinalign": "compare",
    }
    for app_name, input_id in expected.items():
        source = (APPS_DIR / app_name / "app.py").read_text(encoding="utf-8")
        tree = ast.parse(source)

        builds = [
            call
            for name in ("task_button", "background_task_button")
            for call in _calls(tree, name)
            if _first_str(call) == input_id
        ]
        assert builds, f"{app_name}.{input_id} is not a task button"

        # No task button may rely on disabled=, because it does not work.
        for call in builds:
            assert not any(keyword.arg == "disabled" for keyword in call.keywords), (
                f"{app_name}.{input_id} passes disabled= to a task button, which the "
                "bslib component discards -- wrap it in gate() instead"
            )

        # And the working gate must be present, wrapping that button.
        wrapped = [
            call
            for call in _calls(tree, "gate")
            if any(
                isinstance(inner, ast.Call)
                and _first_str(inner) == input_id
                for inner in ast.walk(call)
            )
        ]
        assert wrapped, f"{app_name}.{input_id} is no longer wrapped in gate()"
        assert any(
            keyword.arg == "ready" for call in wrapped for keyword in call.keywords
        ), f"{app_name}.{input_id} is gated without a ready= condition"


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_no_task_button_anywhere_relies_on_disabled(path: Path) -> None:
    """The whole class of bug, not just the two known sites."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    offenders = sorted(
        _first_str(call) or "?"
        for name in ("task_button", "background_task_button")
        for call in _calls(tree, name)
        if any(keyword.arg == "disabled" for keyword in call.keywords)
    )
    assert not offenders, (
        f"{path.parent.name}: {offenders} pass disabled= to a task button; the bslib "
        "component discards it, so the control renders clickable. Use gate()."
    )


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.parent.name)
def test_an_app_that_gates_ships_the_stylesheet(path: Path) -> None:
    """gate() relies on a class; without the CSS the control looks and acts enabled."""
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    if not _calls(tree, "gate"):
        pytest.skip("this app gates nothing")
    assert "INERT_CSS" in source, (
        f"{path.parent.name} calls gate() but never includes INERT_CSS, so the inert "
        "class has no styling and the control stays clickable"
    )
