from __future__ import annotations

import html as html_module
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple

from common.topology_bridge import build_view, load_structure
from topology_viewer.app import _chain_choices, preferred_chain

FIXTURES = Path(__file__).resolve().parents[2] / "notebooks" / "topology_viewer" / "fixtures"


def _inner_document(view_html: str) -> str:
    """The document build_view escaped into the iframe's ``srcdoc`` attribute.

    Assertions must run against this, not against the wrapper: every ``<`` in the
    view arrives as ``&lt;``, so a naive ``"<div class=..." not in html`` check
    passes no matter what the renderer emitted.
    """
    match = re.search(r'srcdoc="(.*?)"\s+style=', view_html, re.DOTALL)
    assert match, "build_view did not return an iframe with a srcdoc attribute"
    return html_module.unescape(match.group(1))


def _payload(view_html: str) -> dict:
    """The JSON blob the browser consumes.

    Assert annotation state here rather than by grepping for class names: both
    ``.topo-filters`` and ``.site-chip`` appear unconditionally in TOPOLOGY_CSS and
    inside TOPOLOGY_JS template strings, so a name match proves nothing.
    """
    document = _inner_document(view_html)
    match = re.search(
        r'<script type="application/json" data-role="payload">(.*?)</script>',
        document,
        re.DOTALL,
    )
    assert match, "no payload script found in the rendered view"
    return json.loads(match.group(1))


@dataclass
class _Ref:
    pdb_id: str
    chains: Dict[str, Optional[Tuple[int, int]]] = field(default_factory=dict)


class _Structure:
    def __init__(self, chains: dict[str, list[int]]) -> None:
        self.residues_by_chain = chains


def test_preferred_chain_picks_the_mapped_chain_not_the_biggest() -> None:
    """A complex often carries a larger unrelated chain; defaulting to it would
    draw a topology for the wrong protein."""
    structure = _Structure({"A": [0] * 900, "E": [0] * 300})
    refs = [_Ref("6Q2O", {"E": (1, 300), "F": (1, 300)})]
    assert preferred_chain(refs, "6Q2O", structure) == "E"


def test_preferred_chain_is_case_insensitive_on_the_entry_id() -> None:
    structure = _Structure({"E": [0]})
    assert preferred_chain([_Ref("6q2o", {"E": None})], "6Q2O", structure) == "E"


def test_preferred_chain_skips_chains_absent_from_the_file() -> None:
    structure = _Structure({"B": [0]})
    assert preferred_chain([_Ref("1ABC", {"E": None, "B": None})], "1ABC", structure) == "B"


def test_preferred_chain_returns_none_when_nothing_matches() -> None:
    structure = _Structure({"A": [0]})
    assert preferred_chain([_Ref("9XYZ", {"E": None})], "1ABC", structure) is None
    assert preferred_chain([], "1ABC", structure) is None
    assert preferred_chain(None, "1ABC", structure) is None


def test_chain_choices_inverts_chain_options_into_shiny_form() -> None:
    structure = load_structure(
        (FIXTURES / "annotated.pdb").read_text(), "annotated.pdb"
    )
    choices = _chain_choices(structure)
    # Structure.chain_options() yields (label, value); Shiny selects want value -> label.
    assert set(choices) == set(structure.residues_by_chain)
    for value, label in choices.items():
        assert label.startswith(value)
        assert "residues" in label


def test_build_view_renders_an_iframe_offline_for_every_fixture() -> None:
    """The whole 2D pipeline, exercised through the app's own import graph.

    No network: this is what keeps the suite offline while still proving the
    renderer, the layout engines and both secondary-structure routes work.
    """
    expected_source = {
        "annotated.pdb": "file (HELIX/SHEET)",
        "annotated.cif": "file (_struct_conf)",
        "bare.pdb": "built-in P-SEA",
    }
    for name, ss_source in expected_source.items():
        path = FIXTURES / name
        text = path.read_text()
        structure = load_structure(text, path.name)
        assert structure.ss_source == ss_source

        html = build_view(
            structure,
            structure.default_chain(),
            {"kind": "upload", "data": text, "format": structure.fmt},
            height=600,
        )
        assert "<iframe srcdoc=" in html
        document = _inner_document(html)
        assert 'class="topo-root"' in document
        assert '<div class="topo-filters">' not in document
        payload = _payload(html)
        assert payload["annotations"] is None
        assert payload["stats"]["residues"] > 0


def test_build_view_without_sites_draws_no_annotation_filter_row() -> None:
    """The file-mode guarantee: no sites in, no marks out."""
    path = FIXTURES / "bare.pdb"
    text = path.read_text()
    structure = load_structure(text, path.name)
    html = build_view(
        structure,
        structure.default_chain(),
        {"kind": "upload", "data": text, "format": structure.fmt},
        height=600,
        sites=None,
        numbering=None,
        accession="",
    )
    assert '<div class="topo-filters">' not in _inner_document(html)
    payload = _payload(html)
    assert payload["annotations"] is None
    assert all("sites" not in element for element in payload["elements"])


def test_build_view_with_sites_does_draw_the_filter_row() -> None:
    """Positive control for the two assertions above.

    Without this, ``'<div class="topo-filters">' not in html`` would still pass if
    the renderer stopped emitting the row altogether.
    """
    from common.topology_bridge import annotations_module as ann

    path = FIXTURES / "annotated.pdb"
    text = path.read_text()
    structure = load_structure(text, path.name)
    chain = structure.default_chain()
    positions = [residue.seq for residue in structure.residues_by_chain[chain][:3]]
    sites = [
        ann.Site(position=position, kind="ptm", residue="SER", source="Scop3P")
        for position in positions
    ]

    html = build_view(
        structure,
        chain,
        {"kind": "upload", "data": text, "format": structure.fmt},
        height=600,
        sites=sites,
        numbering=None,
        accession="P07949",
    )
    document = _inner_document(html)
    assert '<div class="topo-filters">' in document
    payload = _payload(html)
    assert payload["annotations"] is not None
    assert payload["annotations"]["accession"] == "P07949"
    assert payload["annotations"]["counts"]["ptm"] == len(positions)
