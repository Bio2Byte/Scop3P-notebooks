"""The PTM-summary policy of the RINAlign app layer.

The service keeps its contract -- fetch_scop3p_ptm_positions raises Scop3PApiError on a
broken endpoint and returns an empty set for an uncovered protein (see
test_rinalign_service). What the app does with those two outcomes is decided by
``_ptm_summary``: a failure degrades to the UniProt positions with an "unavailable"
notice, while an empty set is reported as an answer about the protein.
"""

from __future__ import annotations

from rinalign.app import _ptm_summary


def test_scop3p_failure_degrades_to_uniprot_positions() -> None:
    combined, status_line, lines = _ptm_summary(None, {5, 9})
    assert combined == frozenset({5, 9})
    assert "Scop3P: unavailable" in status_line
    assert "UniProt: 2" in status_line
    assert "total unique: 2" in status_line
    # The detail box tells the user the degradation happened and that retrying helps.
    assert any("UniProt only" in line for line in lines)
    assert any("retry" in line.lower() for line in lines)


def test_uncovered_protein_is_an_answer_not_an_error() -> None:
    combined, status_line, lines = _ptm_summary(set(), {7})
    assert combined == frozenset({7})
    assert "Scop3P: 0" in status_line
    assert "not in Scop3P" in status_line
    assert "unavailable" not in status_line
    assert not any("unavailable" in line for line in lines)


def test_both_sources_union_and_count() -> None:
    combined, status_line, lines = _ptm_summary({1, 2, 3}, {3, 4})
    assert combined == frozenset({1, 2, 3, 4})
    assert "Scop3P: 3" in status_line
    assert "UniProt: 2" in status_line
    assert "total unique: 4" in status_line
    # Something was found, so the user is pointed at the overlay step.
    assert any("Compare" in line for line in lines)


def test_uniprot_toggled_off_reads_as_off_not_zero() -> None:
    combined, status_line, _lines = _ptm_summary({10}, None)
    assert combined == frozenset({10})
    assert "UniProt: off" in status_line


def test_nothing_anywhere_makes_no_overlay_promise() -> None:
    combined, status_line, lines = _ptm_summary(set(), set())
    assert combined == frozenset()
    assert "total unique: 0" in status_line
    assert not any("Compare" in line for line in lines)
