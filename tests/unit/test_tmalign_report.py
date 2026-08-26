"""Parsing TM-align's report.

The alignment always ran and the 3D view always rendered; what was missing was the
answer. The app showed ``report.splitlines()[0]`` -- TM-align's blank first line -- plus a
path to a temp file, so the TM-score, RMSD and aligned length were all discarded. A
structure alignment with no score does not tell the user anything, which is why this
reads to them as "TM-align does not work".

The sample below is a real report from TM-align 20210520 (2IVS chain A against 2IVT
chain A), kept verbatim including the leading blank line that caused the bug.
"""

from __future__ import annotations

import pytest

from common.structure_viz import TMAlignResult, parse_tmalign_report

REPORT = """
 **********************************************************************
 * TM-align (Version 20210520): protein and RNA structure alignment   *
 * References: Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005)  *
 **********************************************************************

Name of Structure_1: /tmp/x/s1.pdb (to be superimposed onto Structure_2)
Name of Structure_2: /tmp/x/s2.pdb
Length of Structure_1: 284 residues
Length of Structure_2: 288 residues

Aligned length= 278, RMSD=   0.84, Seq_ID=n_identical/n_aligned= 0.982
TM-score= 0.96803 (normalized by length of Structure_1: L=284, d0=6.20)
TM-score= 0.95468 (normalized by length of Structure_2: L=288, d0=6.24)
(You should use TM-score normalized by length of the reference structure)

Total CPU time is  0.09 seconds
"""


def test_the_blank_first_line_is_not_the_result() -> None:
    """The regression, stated directly."""
    assert REPORT.splitlines()[0] == "", "the sample no longer reproduces the bug"
    assert parse_tmalign_report(REPORT).summary().strip() != ""


def test_every_field_is_read() -> None:
    result = parse_tmalign_report(REPORT)
    assert result.tm_score_1 == pytest.approx(0.96803)
    assert result.tm_score_2 == pytest.approx(0.95468)
    assert result.aligned_length == 278
    assert result.rmsd == pytest.approx(0.84)
    assert result.sequence_identity == pytest.approx(0.982)
    assert result.length_1 == 284
    assert result.length_2 == 288


def test_the_two_scores_are_not_confused() -> None:
    """They are normalised by different structures and are not interchangeable."""
    result = parse_tmalign_report(REPORT)
    assert result.tm_score_1 != result.tm_score_2
    assert result.tm_score_1 > result.tm_score_2  # as in the sample


def test_the_summary_carries_the_numbers_a_reader_needs() -> None:
    summary = parse_tmalign_report(REPORT).summary()
    for expected in ("0.96803", "278", "0.84", "98.2%"):
        assert expected in summary, f"{expected} missing from the summary"


def test_the_summary_leads_with_the_score() -> None:
    """It is the answer; a reader should not have to hunt for it."""
    assert parse_tmalign_report(REPORT).summary().splitlines()[0].startswith("TM-score")


@pytest.mark.parametrize(
    "score, expected",
    [
        (0.05, "no structural similarity"),
        (0.16, "no structural similarity"),
        (0.30, "some similarity, probably a different fold"),
        (0.49, "some similarity, probably a different fold"),
        (0.50, "same fold"),
        (0.97, "same fold"),
    ],
)
def test_the_score_is_interpreted_against_the_published_thresholds(score, expected) -> None:
    """Xu & Zhang (2010): 0.17 is random, 0.5 is the same fold.

    Stated in the UI because a bare number invites the reader to invent a threshold.
    """
    assert TMAlignResult(tm_score_1=score).interpretation() == expected


def test_a_report_with_no_scores_says_so_rather_than_looking_empty() -> None:
    result = parse_tmalign_report("some unrelated output")
    assert result.tm_score_1 is None
    assert "no parseable scores" in result.summary()


def test_an_empty_report_does_not_raise() -> None:
    for text in ("", None, "\n\n"):
        assert parse_tmalign_report(text).summary()


def test_a_partial_report_yields_what_it_can() -> None:
    """A future version renaming one field should cost that field, not the summary."""
    partial = "Aligned length= 100, RMSD=   1.50, Seq_ID=n_identical/n_aligned= 0.500"
    result = parse_tmalign_report(partial)
    assert result.aligned_length == 100
    assert result.rmsd == pytest.approx(1.50)
    assert result.tm_score_1 is None
    summary = result.summary()
    assert "100" in summary and "1.50" in summary


def test_the_summary_omits_what_it_does_not_have() -> None:
    """No placeholder values: a missing field is absent, not reported as zero."""
    summary = TMAlignResult(tm_score_1=0.8).summary()
    assert "RMSD" not in summary
    assert "Aligned length" not in summary


def test_interpretation_without_any_score() -> None:
    assert TMAlignResult().interpretation() == "no score"


def test_the_better_normalisation_is_the_one_reported() -> None:
    """TM-align tells the reader to use the reference structure's normalisation.

    With no way to know which the user considers the reference, the higher of the two is
    reported and both are shown underneath, so nothing is hidden.
    """
    result = TMAlignResult(tm_score_1=0.40, tm_score_2=0.80)
    first = result.summary().splitlines()[0]
    assert "0.80000" in first
    assert "0.40000" in result.summary()


# --------------------------------------------------------------------------------------
# The wiring, not just the parser
# --------------------------------------------------------------------------------------


def test_the_app_parses_the_whole_report_not_one_line() -> None:
    """A correct parser fed one line is still the original bug.

    The parser tests above pass whether the app hands over the full report or just
    ``report.splitlines()[0]``, so the call site needs pinning separately -- that split was
    exactly the defect.
    """
    import ast
    from pathlib import Path

    source = (
        Path(__file__).resolve().parents[2] / "apps" / "structure_viz" / "app.py"
    ).read_text(encoding="utf-8")
    tree = ast.parse(source)

    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "parse_tmalign_report"
    ]
    assert calls, "the app no longer parses the TM-align report at all"
    for call in calls:
        assert call.args, "parse_tmalign_report called with no argument"
        argument = call.args[0]
        assert isinstance(argument, ast.Name) and argument.id == "report", (
            "the app passes a slice of the report rather than the whole thing; the "
            "scores live on lines 13-15 and the first line is blank"
        )


def test_the_app_shows_the_summary_to_the_user() -> None:
    """Parsing it and not displaying it would be the same bug with more code."""
    from pathlib import Path

    source = (
        Path(__file__).resolve().parents[2] / "apps" / "structure_viz" / "app.py"
    ).read_text(encoding="utf-8")
    assert "result.summary()" in source


# --------------------------------------------------------------------------------------
# The superposition view
# --------------------------------------------------------------------------------------
# A superposition is two structures. TM-align writes aligned.pdb (structure 1 rotated onto
# structure 2) and leaves structure 2 as the unmoved input -- its own PyMOL script loads
# both. Rendering only aligned.pdb drew a single shape, which cannot show where two
# structures differ.

SUPERPOSED = "ATOM      1  CA  GLU A 713      12.260  -4.242   0.406  1.00  0.00           C\n"
REFERENCE = "ATOM      1  CA  GLU A 713      11.100  -3.900   0.900  1.00  0.00           C\n"


def _view() -> str:
    from common.structure_viz import StructureViewerBuilder

    return StructureViewerBuilder.superposition_html(
        SUPERPOSED, REFERENCE, "2IVS chain A", "7NZN chain A"
    )


def test_both_structures_reach_the_viewer() -> None:
    """The defect, stated as a test: one structure is not a superposition."""
    html = _view()
    assert "12.260" in html, "the superposed structure is missing"
    assert "11.100" in html, "the reference structure is missing"


def test_the_two_structures_are_drawn_in_different_colours() -> None:
    """Both grey would be one indistinguishable shape, which defeats the purpose.

    The colour is passed per structure into the loader, so that is where to look -- not in
    addRepresentation, which receives it as a variable.
    """
    import re

    html = _view()
    colours = set(re.findall(r"load\(\w+, '(#[0-9a-fA-F]{6})'", html))
    assert len(colours) == 2, f"expected two distinct cartoon colours, got {colours}"


def test_the_legend_names_both_structures() -> None:
    """A two-colour picture with no key does not say which is which."""
    html = _view()
    assert "2IVS chain A" in html
    assert "7NZN chain A" in html


def test_the_view_is_framed_only_after_both_have_loaded() -> None:
    """autoView on the first arrival leaves the other outside the camera.

    Matched on the call rather than the bare word: the comment above it in the generated
    source also says "autoView", and a substring search finds that first.
    """
    html = _view()
    assert "Promise.all" in html
    assert html.index("Promise.all") < html.index("stage.autoView()")


def test_the_coordinates_are_json_encoded() -> None:
    """PDB text goes into a JS string literal; a raw splice would break on quotes."""
    html = _view()
    assert "\\n" in html  # newlines escaped rather than literal inside the literal


def test_labels_are_escaped() -> None:
    from common.structure_viz import StructureViewerBuilder

    html = StructureViewerBuilder.superposition_html(
        SUPERPOSED, REFERENCE, "<script>alert(1)</script>", "ok"
    )
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;" in html


def test_region_selection_limits_both_cartoons() -> None:
    from common.structure_viz import StructureViewerBuilder

    html = StructureViewerBuilder.superposition_html(
        SUPERPOSED, REFERENCE, region_sele="700-720"
    )
    assert '"700-720"' in html
    # The default stays the whole protein.
    assert '"protein"' in StructureViewerBuilder.superposition_html(SUPERPOSED, REFERENCE)


def test_site_overlay_reaches_the_viewer_with_style_and_targets() -> None:
    from common.structure_viz import StructureViewerBuilder

    html = StructureViewerBuilder.superposition_html(
        SUPERPOSED,
        REFERENCE,
        site_sele="713 or 720-722",
        site_rep="spacefill",
        sites_on_superposed=True,
        sites_on_reference=False,
        settings_note="Sites: PTMs | region: aligned",
    )
    assert '"713 or 720-722"' in html
    assert '"spacefill"' in html
    # Per-component gating is data the script reads, not duplicated code paths.
    assert "load(superposedText, '#d62728', 0.85, true)" in html
    assert "load(referenceText, '#1f77b4', 0.6, false)" in html
    assert "selected sites" in html
    assert "Sites: PTMs | region: aligned" in html


def test_no_site_selection_means_no_site_legend() -> None:
    from common.structure_viz import StructureViewerBuilder

    html = StructureViewerBuilder.superposition_html(SUPERPOSED, REFERENCE)
    assert "selected sites" not in html


def test_colour_roles_match_the_notebook() -> None:
    """Red = structure 1 (superposed), blue = structure 2 (reference)."""
    from common.structure_viz import StructureViewerBuilder

    html = StructureViewerBuilder.superposition_html(SUPERPOSED, REFERENCE)
    assert "load(superposedText, '#d62728'" in html
    assert "load(referenceText, '#1f77b4'" in html


def test_the_app_renders_a_superposition_not_a_single_structure() -> None:
    """Pins the call sites: the single-structure viewer cannot express this.

    Rendering lives in _render_tm_view (shared by the run handler and the
    highlight-setting redraw), so that is where superposition_html must appear; the
    run handler must delegate to it.
    """
    import ast
    from pathlib import Path

    source = (
        Path(__file__).resolve().parents[2] / "apps" / "structure_viz" / "app.py"
    ).read_text(encoding="utf-8")
    tree = ast.parse(source)

    def _body(name: str) -> str:
        node = next(
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef) and node.name == name
        )
        return ast.get_source_segment(source, node) or ""

    renderer = _body("_render_tm_view")
    handler = _body("_run_tmalign")
    assert "superposition_html" in renderer, "TM-align no longer renders a superposition"
    assert "_render_tm_view" in handler, "the run handler no longer draws its result"
    for body in (renderer, handler):
        assert "ptm_html" not in body, (
            "the single-structure viewer is back in the TM-align path; it paints one "
            "uniform-grey structure and cannot show where two differ"
        )
    assert "alignment.reference" in body, "the reference structure is not being rendered"
