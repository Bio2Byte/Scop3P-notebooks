"""Integrity tests for the RINAlign view builders.

Cell 6 of the notebook is ~50 kB of JavaScript held in Python f-strings with
hand-doubled ``{{``/``}}``. These tests are the guard against that being reformatted,
"tidied", or partially re-typed. They deliberately assert on structure and on the
2D/3D bridge rather than on exact bytes, so legitimate edits stay possible while
silent corruption does not.
"""

from __future__ import annotations

import re
from pathlib import Path

import networkx as nx
import pytest

from common import rinalign_views as views
from common.rinalign import align_rins, build_rin, diff_rins

FIXTURES = Path(__file__).resolve().parents[2] / "notebooks" / "topology_viewer" / "fixtures"

# Every token the two-way network <-> 3D bridge is made of. All six must land in ONE
# document: the force graph and the NGL stage share a single ``window``, and that is
# the whole mechanism. Splitting the linked view across two iframes breaks it
# silently -- nothing errors, clicks simply stop crossing.
BRIDGE_TOKENS = (
    "__RIN_HL",        # NGL click -> network highlight
    "__RIN_ONSELECT",  # network click -> NGL highlight
    "LVgo",            # dropdown -> both
    "LVSEARCH",        # the dropdown itself
    "NGL.Stage",       # the 3D half
    "d3.forceSimulation",  # the network half
)


@pytest.fixture(scope="module")
def diff_result():
    """Two networks over the same coordinates at different cutoffs.

    Same positions, so everything matches; the tighter cutoff drops contacts, which
    yields a real set of lost edges without needing network access.
    """
    left, _ = build_rin(FIXTURES / "annotated.pdb", cutoff=8.0)
    right, _ = build_rin(FIXTURES / "annotated.pdb", cutoff=7.0)
    return diff_rins(left, right), left, right


def _balanced(text: str, opener: str, closer: str) -> bool:
    depth = 0
    for character in text:
        if character == opener:
            depth += 1
        elif character == closer:
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


def _scripts(html: str) -> list[str]:
    return re.findall(r"<script[^>]*>(.*?)</script>", html, re.DOTALL)


# ---------------------------------------------------------------------------
# Every builder renders, and renders something structural
# ---------------------------------------------------------------------------


def _all_views(diff, left, right, **kwargs) -> dict[str, str]:
    pdb_text = (FIXTURES / "annotated.pdb").read_text()
    return {
        "summary": views.summary_html(diff, "L", "R"),
        "contact_map": views.contact_map_html(diff, "L", "R", **kwargs),
        "aligned": views.aligned_network_html(diff, left, right, "L", "R", **kwargs),
        "force": views.force_network_html(diff, left, right, "L", "R", **kwargs),
        "linked": views.linked_view_html(
            diff, left, right, "L", "R", pdb_text, "pdb", "A", **kwargs
        ),
    }


ANCHORS = {
    "summary": ["Jaccard", "Conserved", "Lost", "Gained"],
    "contact_map": ["getContext('2d')", "addEventListener('wheel'", "<canvas"],
    "aligned": ["<svg", "viewBox"],
    "force": ["d3.forceSimulation", "forceLink", "cdnjs.cloudflare.com/ajax/libs/d3"],
    "linked": list(BRIDGE_TOKENS),
}


@pytest.mark.parametrize("name", sorted(ANCHORS))
def test_each_view_renders_with_its_structural_anchors(diff_result, name: str) -> None:
    """Catches an f-string whose interpolation silently stopped producing markup.

    A NameError from a broken ``{{`` escape surfaces here rather than in the browser.
    """
    diff, left, right = diff_result
    html = _all_views(diff, left, right)[name]
    assert html
    for anchor in ANCHORS[name]:
        assert anchor in html, f"{name} lost its {anchor!r} anchor"


@pytest.mark.parametrize("name", sorted(ANCHORS))
def test_view_braces_and_parentheses_balance(diff_result, name: str) -> None:
    """Brace-escape damage shows up as unbalanced JavaScript.

    Note that ``}}`` legitimately appears in the rendered output -- adjacent block
    closes such as ``{d.fx=null;d.fy=null;}}`` -- so counting doubled braces would
    produce false alarms. Balance is the meaningful invariant.
    """
    diff, left, right = diff_result
    html = _all_views(diff, left, right)[name]
    for script in _scripts(html):
        assert _balanced(script, "{", "}"), f"{name}: unbalanced braces in a script block"
        assert _balanced(script, "(", ")"), f"{name}: unbalanced parens in a script block"


@pytest.mark.parametrize("name", sorted(ANCHORS))
def test_views_render_with_and_without_annotation_overlays(diff_result, name: str) -> None:
    diff, left, right = diff_result
    bare = _all_views(diff, left, right)[name]
    overlaid = _all_views(diff, left, right, ptm_pos={10, 25}, var_pos={33})[name]
    assert bare and overlaid

    # summary takes no overlay arguments, and the contact map accepts them but has
    # no overlay layer -- a documented gap inherited from the notebook, which
    # computed the position lists and then never interpolated them. The other three
    # views must visibly change.
    if name in {"summary", "contact_map"}:
        assert overlaid == bare
    else:
        assert overlaid != bare, f"{name} ignored the PTM/variant overlay"


# ---------------------------------------------------------------------------
# The bridge
# ---------------------------------------------------------------------------


def test_linked_view_keeps_the_whole_bridge_in_one_document(diff_result) -> None:
    """The executable form of the migration's central claim.

    ``linked_view_html`` returns ``head + force_network_html(...) + mid + script``.
    One string becomes one ``srcdoc`` becomes one ``window``, which is why the
    ``window.__RIN_HL`` / ``__RIN_ONSELECT`` handshake survives the port untouched
    and needs no ``postMessage`` or ``Shiny.setInputValue`` plumbing.
    """
    diff, left, right = diff_result
    html = views.linked_view_html(
        diff, left, right, "L", "R", (FIXTURES / "annotated.pdb").read_text(), "pdb", "A"
    )
    for token in BRIDGE_TOKENS:
        assert token in html, f"the bridge lost {token!r}"


def test_force_and_linked_views_must_not_share_a_document(diff_result) -> None:
    """They compute the same uid from the same data, so their ids and
    ``window['<uid>*']`` functions collide. Each needs its own iframe, and the app
    puts them in separate nav panels for that reason."""
    diff, left, right = diff_result
    force = views.force_network_html(diff, left, right, "L", "R")
    linked = views.linked_view_html(
        diff, left, right, "L", "R", (FIXTURES / "annotated.pdb").read_text(), "pdb", "A"
    )
    force_uid = re.search(r"window\['(fn[0-9a-f]{6})HL'\]", force)
    linked_uid = re.search(r"window\['(fn[0-9a-f]{6})HL'\]", linked)
    assert force_uid and linked_uid
    assert force_uid.group(1) == linked_uid.group(1)


def test_linked_view_embeds_the_structure_text(diff_result) -> None:
    diff, left, right = diff_result
    html = views.linked_view_html(
        diff, left, right, "L", "R", "ATOM      1  CA  ALA A   1\n", "pdb", "A"
    )
    assert "ATOM" in html
    assert "'pdb'" in html or '"pdb"' in html


# ---------------------------------------------------------------------------
# Determinism and the document shell
# ---------------------------------------------------------------------------


def test_element_ids_are_stable_across_calls(diff_result) -> None:
    """``hash()`` of a str is PYTHONHASHSEED-salted; md5 is not.

    Without this the ids changed on every interpreter restart, so none of these
    views could be compared between runs.
    """
    diff, left, right = diff_result
    assert views.contact_map_html(diff, "L", "R") == views.contact_map_html(diff, "L", "R")
    assert views.force_network_html(diff, left, right, "L", "R") == views.force_network_html(
        diff, left, right, "L", "R"
    )


def test_stable_uid_is_six_hex_characters() -> None:
    uid = views._stable_uid([1, 2, 3])
    assert re.fullmatch(r"[0-9a-f]{6}", uid)
    assert uid == views._stable_uid([1, 2, 3])
    assert uid != views._stable_uid([1, 2, 4])


def test_html_document_wraps_without_escaping() -> None:
    """Shiny escapes ``srcdoc`` itself; escaping here too would double-escape.

    This is the one behaviour that differs from the notebook's ``_voila_iframe``.
    """
    document = views.html_document("<div id='x'>&amp;</div><script>var a={};</script>")
    assert document.startswith("<!DOCTYPE html>")
    assert "<div id='x'>&amp;</div>" in document
    assert "<script>var a={};</script>" in document
    assert "&lt;" not in document


def test_voila_iframe_helper_was_not_ported() -> None:
    assert not hasattr(views, "_voila_iframe")


# ---------------------------------------------------------------------------
# Degenerate inputs
# ---------------------------------------------------------------------------


def test_contact_map_reports_no_overlap_instead_of_raising() -> None:
    left = nx.Graph()
    left.add_node("ALA_1", position=1, resname="ALA")
    right = nx.Graph()
    right.add_node("ALA_9", position=9, resname="ALA")
    diff = diff_rins(left, right)
    assert diff["matched_pos"] == []
    assert "No overlapping residues" in views.contact_map_html(diff, "L", "R")


def test_summary_renders_for_identical_networks(diff_result) -> None:
    _diff, left, _right = diff_result
    identical = diff_rins(left, left.copy())
    html = views.summary_html(identical, "L", "L")
    assert "1.000" in html


def test_alignment_result_is_not_fed_to_the_diff_views(diff_result) -> None:
    """align_rins returns a different shape; the app must branch on mode, and these
    keys are what it branches on."""
    _diff, left, right = diff_result
    alignment = align_rins(left, right)
    assert set(alignment) == {"mapping", "conserved", "only_G1", "only_G2", "jaccard"}
    assert "matched_pos" not in alignment
    assert "residue_impact" not in alignment
