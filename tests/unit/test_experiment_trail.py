"""The experiment trail: an ordered, per-session record of a protocol run.

The trail is a deliverable, not debug output -- it is what makes a figure's provenance
readable afterwards. So the properties that make it trustworthy are pinned here: steps
are ordered, sessions do not interleave, the vocabulary is closed, and a guess never
appears at a level that implies success.
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import pytest

from common.logging_utils import (
    TRAIL_ACTIONS,
    ExperimentTrail,
    _reset_logging_for_tests,
    configure_logging,
    get_trail_file_path,
    new_trail,
    quiet_third_party,
)


@pytest.fixture
def trail() -> ExperimentTrail:
    return ExperimentTrail("structure-viz", session_id="sess0001")


# --------------------------------------------------------------------------------------
# Ordering and attribution
# --------------------------------------------------------------------------------------


def test_steps_are_numbered_in_order(trail, caplog) -> None:
    with caplog.at_level(logging.INFO, logger="scop3p.trail"):
        trail.opened("Structure Visualisation")
        trail.entered("UniProtKB accession", "P07949")
        trail.clicked("Set protein")
    steps = [record.getMessage() for record in caplog.records]
    assert "step=1" in steps[0] and "action=open" in steps[0]
    assert "step=2" in steps[1] and "action=input" in steps[1]
    assert "step=3" in steps[2] and "action=click" in steps[2]


def test_step_numbering_does_not_interleave_across_sessions(caplog) -> None:
    """One process serves many browser sessions; two tabs must stay separable."""
    first = ExperimentTrail("structure-viz", session_id="aaaa1111")
    second = ExperimentTrail("structure-viz", session_id="bbbb2222")
    with caplog.at_level(logging.INFO, logger="scop3p.trail"):
        first.clicked("Set protein")
        second.clicked("Set protein")
        first.clicked("Fetch PTMs")
    messages = [record.getMessage() for record in caplog.records]
    assert "session=aaaa1111 step=1" in messages[0]
    assert "session=bbbb2222 step=1" in messages[1]
    assert "session=aaaa1111 step=2" in messages[2]


def test_each_trail_gets_its_own_session_id() -> None:
    assert new_trail("x").session_id != new_trail("x").session_id


def test_step_returns_its_own_number(trail) -> None:
    assert [trail.clicked("a"), trail.clicked("b")] == [1, 2]
    assert trail.steps_recorded == 2


# --------------------------------------------------------------------------------------
# Levels: the record must not flatten success and failure together
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "call, expected_level",
    [
        (lambda t: t.opened("P"), logging.INFO),
        (lambda t: t.entered("f", "v"), logging.INFO),
        (lambda t: t.selected("f", "v"), logging.INFO),
        (lambda t: t.clicked("b"), logging.INFO),
        (lambda t: t.produced("r"), logging.INFO),
        (lambda t: t.exported("f"), logging.INFO),
        (lambda t: t.blocked("why"), logging.WARNING),
        (lambda t: t.failed("what"), logging.ERROR),
    ],
)
def test_levels_match_the_meaning(trail, caplog, call, expected_level) -> None:
    with caplog.at_level(logging.DEBUG, logger="scop3p.trail"):
        call(trail)
    assert caplog.records[-1].levelno == expected_level


def test_a_blocked_step_is_not_reported_as_success(trail, caplog) -> None:
    """A blocked action at INFO would read as though it had run."""
    with caplog.at_level(logging.DEBUG, logger="scop3p.trail"):
        trail.blocked("no cached graph")
    assert caplog.records[-1].levelno > logging.INFO


# --------------------------------------------------------------------------------------
# Vocabulary
# --------------------------------------------------------------------------------------


def test_the_named_wrappers_only_emit_known_actions(trail, caplog) -> None:
    """Every convenience method must stay inside the closed vocabulary."""
    with caplog.at_level(logging.DEBUG, logger="scop3p.trail"):
        trail.opened("P")
        trail.entered("f", "v")
        trail.selected("f", "v")
        trail.clicked("b")
        trail.produced("r")
        trail.blocked("w")
        trail.failed("w")
        trail.exported("f")
    emitted = {
        message.split("action=")[1].split()[0]
        for message in (record.getMessage() for record in caplog.records)
    }
    assert emitted <= set(TRAIL_ACTIONS)
    assert emitted == set(TRAIL_ACTIONS), "a wrapper for some action is missing or unused"


def test_an_unknown_action_is_reported_but_does_not_raise(trail, caplog) -> None:
    """Losing a log line must never break a protocol -- but it must be visible."""
    with caplog.at_level(logging.DEBUG, logger="scop3p.trail"):
        assert trail.step("teleport", "something odd") == 1
    assert any(record.levelno == logging.WARNING for record in caplog.records)
    assert any("teleport" in record.getMessage() for record in caplog.records)


# --------------------------------------------------------------------------------------
# Formatting: readable by a person, parseable by a script
# --------------------------------------------------------------------------------------


def test_values_with_spaces_are_quoted(trail, caplog) -> None:
    with caplog.at_level(logging.INFO, logger="scop3p.trail"):
        trail.produced("RIN built for chain A", nodes=283)
    message = caplog.records[-1].getMessage()
    assert 'detail="RIN built for chain A"' in message
    assert "nodes=283" in message


def test_simple_values_are_not_quoted(trail, caplog) -> None:
    with caplog.at_level(logging.INFO, logger="scop3p.trail"):
        trail.selected("PDB entry", "2IVT")
    assert "detail=" in caplog.records[-1].getMessage()
    assert 'numbering="' not in caplog.records[-1].getMessage()


def test_extra_fields_are_ordered_so_lines_diff_cleanly(trail, caplog) -> None:
    with caplog.at_level(logging.INFO, logger="scop3p.trail"):
        trail.produced("built", nodes=1, cutoff=8.0, edges=2)
    message = caplog.records[-1].getMessage()
    assert message.index("cutoff=") < message.index("edges=") < message.index("nodes=")


@pytest.mark.parametrize("value", [None, ""])
def test_an_empty_value_renders_as_a_placeholder(trail, caplog, value) -> None:
    """An empty value must not produce a dangling "field=" that breaks parsing."""
    with caplog.at_level(logging.INFO, logger="scop3p.trail"):
        trail.entered("accession", value)
    assert "= -" in caplog.records[-1].getMessage()


# --------------------------------------------------------------------------------------
# The standalone record file
# --------------------------------------------------------------------------------------


def test_the_trail_is_written_to_its_own_file(tmp_path, monkeypatch) -> None:
    _reset_logging_for_tests()
    monkeypatch.setenv("SCOP3P_LOG_DIR", str(tmp_path))
    try:
        configure_logging()
        trail = new_trail("structure-viz", "sess0001")
        trail.opened("Structure Visualisation")
        trail.entered("UniProtKB accession", "P07949")
        logging.getLogger("scop3p.structure_viz").info(
            "a diagnostic", extra={"event": "diag"}
        )

        trail_path = get_trail_file_path()
        assert trail_path is not None
        recorded = Path(trail_path).read_text(encoding="utf-8")
        assert "step=1" in recorded and "step=2" in recorded
        # The record is the narrative only; diagnostics stay in the combined log.
        assert "a diagnostic" not in recorded
    finally:
        _reset_logging_for_tests()


def test_the_trail_survives_the_log_level_being_raised(tmp_path, monkeypatch) -> None:
    """SCOP3P_LOG_LEVEL is for quieting diagnostics, not for losing the record."""
    _reset_logging_for_tests()
    monkeypatch.setenv("SCOP3P_LOG_DIR", str(tmp_path))
    monkeypatch.setenv("SCOP3P_LOG_LEVEL", "WARNING")
    try:
        configure_logging()
        new_trail("structure-viz", "sess0001").clicked("Set protein")
        recorded = Path(get_trail_file_path()).read_text(encoding="utf-8")
        assert "Set protein" in recorded
    finally:
        _reset_logging_for_tests()


# --------------------------------------------------------------------------------------
# Keeping third-party chatter out
# --------------------------------------------------------------------------------------


def test_stdout_is_captured_to_debug(caplog) -> None:
    logger = logging.getLogger("scop3p.test.quiet")
    with caplog.at_level(logging.DEBUG, logger="scop3p.test.quiet"):
        with quiet_third_party(logger):
            print("Reading input fasta...")
    assert any("Reading input fasta" in record.getMessage() for record in caplog.records)
    assert all(record.levelno == logging.DEBUG for record in caplog.records)


def test_warnings_are_captured_to_debug(caplog) -> None:
    logger = logging.getLogger("scop3p.test.quiet")
    with caplog.at_level(logging.DEBUG, logger="scop3p.test.quiet"):
        with quiet_third_party(logger):
            warnings.warn("InconsistentVersionWarning: SVC 1.0.2 vs 1.4.2")
    assert any("SVC 1.0.2" in record.getMessage() for record in caplog.records)


def test_nothing_reaches_the_log_at_info(caplog) -> None:
    """The whole point: at the default level this output is invisible."""
    logger = logging.getLogger("scop3p.test.quiet")
    with caplog.at_level(logging.INFO, logger="scop3p.test.quiet"):
        with quiet_third_party(logger):
            print("starting crunch 1")
            warnings.warn("noisy")
    assert caplog.records == []


def test_stdout_is_restored_afterwards() -> None:
    import sys

    before = sys.stdout
    with quiet_third_party(logging.getLogger("scop3p.test.quiet")):
        pass
    assert sys.stdout is before


def test_an_exception_propagates_and_still_flushes_what_was_printed(caplog) -> None:
    """A failing prediction must not lose the output that explains why."""
    logger = logging.getLogger("scop3p.test.quiet")
    with caplog.at_level(logging.DEBUG, logger="scop3p.test.quiet"):
        with pytest.raises(RuntimeError):
            with quiet_third_party(logger):
                print("got as far as loading the model")
                raise RuntimeError("predictor exploded")
    assert any("loading the model" in record.getMessage() for record in caplog.records)


def test_stdout_is_restored_after_an_exception() -> None:
    import sys

    before = sys.stdout
    with pytest.raises(RuntimeError):
        with quiet_third_party(logging.getLogger("scop3p.test.quiet")):
            raise RuntimeError("boom")
    assert sys.stdout is before


def test_blank_lines_are_not_recorded(caplog) -> None:
    logger = logging.getLogger("scop3p.test.quiet")
    with caplog.at_level(logging.DEBUG, logger="scop3p.test.quiet"):
        with quiet_third_party(logger):
            print("\n\n   \n")
    assert caplog.records == []
