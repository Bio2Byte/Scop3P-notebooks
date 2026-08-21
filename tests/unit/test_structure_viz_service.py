from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pytest
import requests

from common.structure_viz import (
    B2B_METRIC_COLUMNS,
    B2B_THRESHOLDS,
    PTM_TABLE_COLUMNS,
    StructureOps,
    StructureVizService,
    b2b_interpretation,
    b2b_legend_html,
    b2b_value_range,
    build_ptm_table,
    chain_choices_for_pdb,
    merge_ptm_tables,
    numeric_b2b_columns,
    parse_pdb_xrefs,
    pdb_entry_choices,
    pseudocolor,
    residue_three_letter,
    uniprot_range_for_chain,
)


def _fill_colours(html: str) -> set[str]:
    """Every node fill colour the generated pyvis document mentions."""
    return {match.lower() for match in re.findall(r'"color":\s*"(#[0-9a-fA-F]{6})"', html)}


class _TextResponse:
    def __init__(self, text: str) -> None:
        self.text = text

    def raise_for_status(self) -> None:
        return None


class _JsonResponse:
    """Minimal stand-in for a Scop3P v1 JSON response."""

    def __init__(self, payload) -> None:  # noqa: ANN001
        self._payload = payload
        self.status_code = 200
        self.headers = {"content-type": "application/json"}

    def json(self):  # noqa: ANN201
        return self._payload

    def raise_for_status(self) -> None:
        return None


def test_fetch_ptms_normalizes_position_column(monkeypatch, tmp_path: Path) -> None:
    # v1 field names; Scop3PClient maps them to position/residue/name.
    monkeypatch.setattr(
        "common.services.requests.get",
        lambda *args, **kwargs: _JsonResponse(
            [
                {"modified_residue": "S", "uniprot_position": "10", "modification_name": "Phosphorylation"},
                {"modified_residue": "T", "uniprot_position": "bad", "modification_name": "Phosphorylation"},
            ]
        ),
    )

    dataframe = StructureVizService(tmp_path).fetch_ptms("P12345")
    assert dataframe["position"].tolist() == [10, pd.NA]
    assert str(dataframe["position"].dtype) == "Int64"


def test_fetch_ptms_returns_empty_dataframe(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        "common.services.requests.get",
        lambda *args, **kwargs: _JsonResponse([]),
    )

    dataframe = StructureVizService(tmp_path).fetch_ptms("P12345")
    assert dataframe.empty


def test_fetch_sequence_parses_fasta(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        "common.structure_viz.requests.get",
        lambda *args, **kwargs: _TextResponse(">sp|P12345|\nACD\nEFG\n"),
    )

    sequence = StructureVizService(tmp_path).fetch_sequence("P12345")
    assert sequence == "ACDEFG"


def test_normalize_b2b_prediction_uses_expected_columns_and_types(tmp_path: Path) -> None:
    service = StructureVizService(tmp_path)
    dataframe = service._normalize_b2b_prediction(
        {
            "seq": "AC",
            "backbone": ["0.1", "0.2"],
            "sidechain": ["0.3", "0.4"],
            "ppII": ["0.5", "0.6"],
            "coil": ["0.7", "0.8"],
            "sheet": ["0.9", "1.0"],
            "helix": ["1.1", "1.2"],
            "earlyFolding": ["1.3", "1.4"],
            "disoMine": ["1.5", "1.6"],
        }
    )

    assert list(dataframe.columns) == [
        "Position",
        "Amino acid",
        *B2B_METRIC_COLUMNS,
        *(f"{metric}_normalized" for metric in B2B_METRIC_COLUMNS),
    ]
    assert dataframe["Position"].tolist() == [1, 2]
    assert dataframe["Amino acid"].tolist() == ["AC", "AC"]
    assert dataframe.loc[0, "backbone"] == 0.1
    assert dataframe.loc[1, "disoMine"] == 1.6
    assert dataframe.loc[0, "backbone_normalized"] == 0.0
    assert dataframe.loc[1, "backbone_normalized"] == 1.0
    assert dataframe.loc[0, "disoMine_normalized"] == 0.0
    assert dataframe.loc[1, "disoMine_normalized"] == 1.0


def test_normalize_b2b_prediction_tolerates_missing_and_uneven_fields(tmp_path: Path) -> None:
    service = StructureVizService(tmp_path)
    dataframe = service._normalize_b2b_prediction(
        {
            "seq": "ACD",
            "backbone": [0.1, 0.2, 0.3],
            "sidechain": [0.4],
            "ppII": None,
            "coil": [0.5, 0.6, 0.7, 0.8],
            "earlyFolding": [0.9, 1.0],
        }
    )

    assert dataframe["Position"].tolist() == [1, 2, 3]
    assert dataframe["Amino acid"].tolist() == ["ACD", "ACD", "ACD"]
    assert dataframe.loc[0, "sidechain"] == 0.4
    assert pd.isna(dataframe.loc[1, "sidechain"])
    assert pd.isna(dataframe.loc[2, "sidechain"])
    assert dataframe["coil"].tolist() == [0.5, 0.6, 0.7]
    assert dataframe["ppII"].isna().all()
    assert dataframe["ppII_normalized"].isna().all()
    assert dataframe.loc[0, "coil_normalized"] == 0.0
    assert dataframe.loc[2, "coil_normalized"] == 1.0


def test_normalize_b2b_prediction_uses_zero_for_constant_metric_series(tmp_path: Path) -> None:
    service = StructureVizService(tmp_path)
    dataframe = service._normalize_b2b_prediction(
        {
            "seq": "AC",
            "backbone": [0.7, 0.7],
        }
    )

    assert dataframe["backbone"].tolist() == [0.7, 0.7]
    assert dataframe["backbone_normalized"].tolist() == [0.0, 0.0]


# ---------------------------------------------------------------------------
# PTM tables: residue normalisation, UniProt source, site-level merge
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value,expected", [
    ("Phosphoserine", "SER"),
    ("phosphotyrosine", "TYR"),
    ("Phosphothreonine", "THR"),
    ("N6-acetyllysine", "LYS"),
    ("SER", "SER"),
    ("ser", "SER"),
    ("S", "SER"),
    ("Y", "TYR"),
    ("", ""),
    ("something unmapped", ""),
])
def test_residue_three_letter(value: str, expected: str) -> None:
    assert residue_three_letter(value) == expected


def test_residue_three_letter_falls_back_to_the_sequence() -> None:
    assert residue_three_letter("", "MKSTY", 3) == "SER"
    assert residue_three_letter("unmapped", "MKSTY", 5) == "TYR"
    assert residue_three_letter("", "MKSTY", 99) == ""
    assert residue_three_letter("", "MKSTY", None) == ""


def test_build_ptm_table_normalises_residue_for_the_viewer_colour_map() -> None:
    """StructureViewerBuilder.ptm_html colours on SER/THR/TYR.

    Scop3P v1 reports the descriptive name ("Phosphoserine"), which matches no key in
    that map, so every PTM fell through to the default colour. The descriptive text is
    kept in `modification` rather than lost.
    """
    frame = build_ptm_table(
        pd.DataFrame([
            {"position": 105, "residue": "Phosphoserine", "name": "phosphorylation", "source": "PRIDE"},
            {"position": 687, "residue": "Phosphotyrosine", "name": "phosphorylation", "source": "UniProt"},
        ]),
        "P07949",
    )
    assert frame["residue"].tolist() == ["SER", "TYR"]
    assert frame["modification"].tolist() == ["Phosphoserine", "Phosphotyrosine"]
    assert frame["ACC_ID"].unique().tolist() == ["P07949"]
    assert list(frame.columns) == list(PTM_TABLE_COLUMNS)
    assert str(frame["position"].dtype) == "Int64"


def test_build_ptm_table_of_nothing_still_has_the_schema() -> None:
    frame = build_ptm_table(pd.DataFrame(), "P1")
    assert frame.empty
    assert list(frame.columns) == list(PTM_TABLE_COLUMNS)


def _ptm(position: int, residue: str, source: str, reference: str = "", evidence: str = ""):
    return {
        "position": position, "residue": residue, "modification": residue, "name": "mod",
        "evidence": evidence, "source": source, "reference": reference,
        "functionalScore": pd.NA,
        "feature_type": source, "ACC_ID": "P1",
    }


def test_merge_lists_a_shared_site_once() -> None:
    scop3p = pd.DataFrame([_ptm(10, "SER", "PRIDE"), _ptm(20, "THR", "PRIDE")])
    uniprot = pd.DataFrame([_ptm(10, "SER", "UniProt"), _ptm(30, "TYR", "UniProt")])
    merged = merge_ptm_tables(scop3p, uniprot)

    assert merged["position"].tolist() == [10, 20, 30]
    shared = merged[merged["position"] == 10].iloc[0]
    assert "PRIDE" in shared["source"] and "UniProt" in shared["source"]


def test_merge_folds_references_without_duplicating_a_citation() -> None:
    """Scop3P cites a bare PMID, UniProt cites PubMed:<pmid> -- the same paper."""
    scop3p = pd.DataFrame([_ptm(10, "SER", "PRIDE", reference="24560924")])
    uniprot = pd.DataFrame([_ptm(10, "SER", "UniProt", reference="PubMed:24560924")])
    merged = merge_ptm_tables(scop3p, uniprot)
    assert len(merged) == 1
    assert merged.iloc[0]["reference"] == "24560924"


def test_merge_borrows_uniprot_evidence_only_when_scop3p_has_none() -> None:
    scop3p = pd.DataFrame([_ptm(10, "SER", "PRIDE", evidence=""),
                           _ptm(20, "THR", "PRIDE", evidence="Experimental")])
    uniprot = pd.DataFrame([_ptm(10, "SER", "UniProt", evidence="ECO:0000269"),
                            _ptm(20, "THR", "UniProt", evidence="ECO:0000250")])
    merged = merge_ptm_tables(scop3p, uniprot).set_index("position")
    assert merged.at[10, "evidence"] == "ECO:0000269"
    assert merged.at[20, "evidence"] == "Experimental"


def test_merge_does_not_collapse_different_residues_at_one_position() -> None:
    """Site identity is accession + residue + position, so a disagreement about the
    residue is two rows, not one silently merged."""
    scop3p = pd.DataFrame([_ptm(10, "SER", "PRIDE")])
    uniprot = pd.DataFrame([_ptm(10, "THR", "UniProt")])
    assert len(merge_ptm_tables(scop3p, uniprot)) == 2


@pytest.mark.parametrize("scop3p_empty,uniprot_empty", [(True, False), (False, True), (True, True)])
def test_merge_handles_either_source_being_empty(scop3p_empty, uniprot_empty) -> None:
    empty = pd.DataFrame(columns=list(PTM_TABLE_COLUMNS))
    scop3p = empty if scop3p_empty else pd.DataFrame([_ptm(10, "SER", "PRIDE")])
    uniprot = empty if uniprot_empty else pd.DataFrame([_ptm(20, "THR", "UniProt")])
    merged = merge_ptm_tables(scop3p, uniprot)
    assert list(merged.columns) == list(PTM_TABLE_COLUMNS)
    assert len(merged) == (0 if scop3p_empty and uniprot_empty else 1)


def test_merge_treats_none_as_uniprot_disabled() -> None:
    scop3p = pd.DataFrame([_ptm(10, "SER", "PRIDE")])
    merged = merge_ptm_tables(scop3p, None)
    assert len(merged) == 1
    assert merged.iloc[0]["source"] == "PRIDE"


def test_fetch_uniprot_ptms_keeps_single_residue_features_only(monkeypatch, tmp_path: Path) -> None:
    payload = {
        "sequence": "MKSTYA",
        "features": [
            {"category": "PTM", "type": "MOD_RES", "begin": "3", "end": "3",
             "description": "Phosphoserine; by autocatalysis",
             "evidences": [{"code": "ECO:0000269", "source": {"name": "PubMed", "id": "123"}}]},
            {"category": "PTM", "type": "CARBOHYD", "begin": "4", "end": "6",
             "description": "N-linked glycan"},                      # a range: skipped
            {"category": "DOMAIN", "type": "DOMAIN", "begin": "1", "end": "1",
             "description": "Kinase"},                                # wrong category
        ],
    }
    monkeypatch.setattr("common.http_lookup.requests.get",
                        lambda *a, **kw: _JsonResponse(payload))

    frame = StructureVizService(tmp_path).fetch_uniprot_ptms("P1")
    assert len(frame) == 1
    row = frame.iloc[0]
    assert row["position"] == 3
    assert row["residue"] == "SER"
    assert row["name"] == "Phosphoserine"        # trimmed at the semicolon
    assert row["source"] == "UniProt"
    assert row["evidence"] == "ECO:0000269"
    assert row["reference"] == "PubMed:123"
    assert list(frame.columns) == list(PTM_TABLE_COLUMNS)


def test_fetch_uniprot_ptms_with_no_features(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr("common.http_lookup.requests.get",
                        lambda *a, **kw: _JsonResponse({"sequence": "MK", "features": []}))
    frame = StructureVizService(tmp_path).fetch_uniprot_ptms("P1")
    assert frame.empty
    assert list(frame.columns) == list(PTM_TABLE_COLUMNS)


def test_viewer_colour_map_covers_the_residues_the_sources_produce() -> None:
    """The 3D viewer colours on the three-letter code build_ptm_table now emits.

    This is the contract that silently broke when Scop3P v1 started reporting
    "Phosphoserine" instead of a residue code: colorMap[residue] missed every time and
    every PTM rendered in the fallback colour.
    """
    import re as _re

    from common.structure_viz import StructureViewerBuilder

    html = StructureViewerBuilder.ptm_html(
        "ATOM      1  CA  SER A 105\n", "P07949",
        [{"position": 105, "residue": "SER"}, {"position": 98, "residue": "ASN"}],
    )
    match = _re.search(r"const colorMap = (\{.*?\});", html, _re.DOTALL)
    assert match, "colorMap not found in the rendered viewer"
    keys = set(json.loads(match.group(1)))

    # Everything build_ptm_table can emit for the phospho and UniProt sources.
    for residue in ("SER", "THR", "TYR", "ASN", "LYS", "CYS"):
        assert residue in keys, f"{residue} would fall back to the default colour"
    assert residue_three_letter("Phosphoserine") in keys


# --------------------------------------------------------------------------------------
# PDB cross-references: turning a UniProt entry into the structure picker's choices
# --------------------------------------------------------------------------------------


def _xref(pdb_id: str, chains: str, *, method: str = "X-ray", resolution: str | None = "2.00 A"):
    """A UniProt PDB cross-reference as the REST payload spells it."""
    properties = [{"key": "Chains", "value": chains}, {"key": "Method", "value": method}]
    if resolution is not None:
        properties.append({"key": "Resolution", "value": resolution})
    return {"database": "PDB", "id": pdb_id, "properties": properties}


def test_parse_pdb_xrefs_reads_chains_method_and_resolution() -> None:
    refs = parse_pdb_xrefs({"uniProtKBCrossReferences": [_xref("2IVS", "A/B=705-1013")]})
    assert len(refs) == 1
    assert refs[0].pdb_id == "2IVS"
    assert refs[0].method == "X-ray"
    assert refs[0].resolution == "2.00 A"
    # "A/B=..." is one range shared by two chains, not a chain literally named "A/B".
    assert refs[0].chain_ranges == {"A": (705, 1013), "B": (705, 1013)}


def test_parse_pdb_xrefs_keeps_per_chain_ranges_apart() -> None:
    refs = parse_pdb_xrefs({"uniProtKBCrossReferences": [_xref("1ABC", "A=1-200, B=300-400")]})
    assert refs[0].chain_ranges == {"A": (1, 200), "B": (300, 400)}


def test_parse_pdb_xrefs_spans_a_discontinuous_chain() -> None:
    """A chain observed in two segments still has one outer span for slicing."""
    refs = parse_pdb_xrefs({"uniProtKBCrossReferences": [_xref("1ABC", "A=10-50; 70-120")]})
    assert refs[0].chain_ranges == {"A": (10, 120)}


def test_parse_pdb_xrefs_ignores_non_pdb_databases() -> None:
    payload = {
        "uniProtKBCrossReferences": [
            {"database": "AlphaFoldDB", "id": "P07949", "properties": []},
            _xref("2IVS", "A=705-1013"),
        ]
    }
    assert [ref.pdb_id for ref in parse_pdb_xrefs(payload)] == ["2IVS"]


def test_parse_pdb_xrefs_survives_a_reference_with_no_properties() -> None:
    payload = {"uniProtKBCrossReferences": [{"database": "PDB", "id": "9XYZ"}]}
    refs = parse_pdb_xrefs(payload)
    assert refs[0].pdb_id == "9XYZ"
    assert refs[0].chain_ranges == {}


def test_pdb_entry_choices_leads_with_a_placeholder_and_describes_each_entry() -> None:
    refs = parse_pdb_xrefs(
        {"uniProtKBCrossReferences": [_xref("2IVS", "A/B=705-1013"), _xref("2IVT", "A=705-1013")]}
    )
    choices = pdb_entry_choices(refs)
    assert list(choices)[0] == ""
    assert "2IVS" in choices and "2IVT" in choices
    assert "chains A, B" in choices["2IVS"]
    assert "2.00 A" in choices["2IVS"]


def test_pdb_entry_choices_of_nothing_explains_itself() -> None:
    """The picker must never render as an empty box."""
    choices = pdb_entry_choices([])
    assert list(choices) == [""]
    assert choices[""].strip() != ""


def test_chain_choices_for_pdb_labels_the_uniprot_range() -> None:
    refs = parse_pdb_xrefs({"uniProtKBCrossReferences": [_xref("2IVS", "A/B=705-1013")]})
    choices = chain_choices_for_pdb(refs, "2IVS")
    assert choices == {"A": "A (705-1013)", "B": "B (705-1013)"}


def test_chain_choices_for_an_unknown_entry_is_empty() -> None:
    assert chain_choices_for_pdb([], "2IVS") == {}


def test_uniprot_range_for_chain_round_trips() -> None:
    refs = parse_pdb_xrefs({"uniProtKBCrossReferences": [_xref("1ABC", "A=1-200, B=300-400")]})
    assert uniprot_range_for_chain(refs, "1ABC", "B") == (300, 400)
    assert uniprot_range_for_chain(refs, "1ABC", "Z") is None
    assert uniprot_range_for_chain(refs, "9ZZZ", "A") is None


def test_fetch_pdb_xrefs_uses_a_bounded_connect_timeout(monkeypatch, tmp_path: Path) -> None:
    """A slow UniProt must not hold the ASGI loop for the full read timeout."""
    seen: dict = {}

    def fake_get(url, headers=None, timeout=None):  # noqa: ANN001, ARG001
        seen["url"] = url
        seen["timeout"] = timeout
        return _JsonResponse({"uniProtKBCrossReferences": [_xref("2IVS", "A=705-1013")]})

    monkeypatch.setattr("common.http_lookup.requests.get", fake_get)
    service = StructureVizService(tmp_path)
    refs = service.fetch_pdb_xrefs("P07949")

    assert [ref.pdb_id for ref in refs] == ["2IVS"]
    connect, read = seen["timeout"]
    assert connect == StructureVizService.CONNECT_TIMEOUT
    assert read == StructureVizService.LOOKUP_READ_TIMEOUT
    # Both bounds must stay well under the service default, which is sized for file
    # downloads. A lookup inherits that only by mistake.
    assert connect < read < service.timeout


# --------------------------------------------------------------------------------------
# Bio2Byte colouring of the residue interaction network
# --------------------------------------------------------------------------------------


def test_pseudocolor_runs_green_to_red() -> None:
    low = pseudocolor(0.0, 1.0, 0.0)
    high = pseudocolor(0.0, 1.0, 1.0)
    assert low != high
    for colour in (low, high):
        assert colour.startswith("#") and len(colour) == 7
        int(colour[1:], 16)  # parses as hex


def test_pseudocolor_is_stable_when_the_range_is_flat() -> None:
    """A constant series must not divide by zero."""
    assert pseudocolor(0.5, 0.5, 0.5).startswith("#")


def test_pseudocolor_walks_the_hue_sweep_monotonically() -> None:
    """Successive values must not double back on the scale."""
    ramp = [pseudocolor(0.0, 1.0, step / 10) for step in range(11)]
    assert len(set(ramp)) > 5
    assert ramp[0] != ramp[-1]


def test_numeric_b2b_columns_skips_position_and_residue() -> None:
    frame = pd.DataFrame(
        {"Position": [1, 2], "Residue": ["A", "C"], "backbone": [0.1, 0.2], "note": ["x", "y"]}
    )
    columns = numeric_b2b_columns(frame)
    assert "backbone" in columns
    assert "Position" not in columns
    assert "Residue" not in columns
    assert "note" not in columns


def test_numeric_b2b_columns_of_nothing() -> None:
    assert numeric_b2b_columns(None) == []
    assert numeric_b2b_columns(pd.DataFrame()) == []


def test_b2b_value_range_ignores_unparseable_values() -> None:
    frame = pd.DataFrame({"backbone": [0.1, "n/a", 0.9]})
    assert b2b_value_range(frame, "backbone") == (0.1, 0.9)
    assert b2b_value_range(frame, "absent") is None


def test_b2b_legend_html_names_the_metric_and_its_bounds() -> None:
    html = b2b_legend_html("backbone", 0.25, 0.95)
    assert "backbone" in html
    assert "0.25" in html and "0.95" in html


@pytest.mark.parametrize("metric", sorted(B2B_THRESHOLDS))
def test_b2b_bands_change_where_mutation_effect_changes_its_label(metric: str) -> None:
    """The two apps must interpret the same predictor identically.

    Both surfaces describe Bio2Byte output to the same users, so if these drift one app
    calls a residue rigid while the other calls it flexible. mutation_effect keeps its
    boundaries as literals inside its label functions rather than as a table, so this
    compares the *boundaries* -- the values at which each side switches band -- instead
    of the label strings, which are worded differently on purpose ("context dependent"
    against "context-dependent").
    """
    from common.mutation_effect import MutationEffectInference

    label = dict(MutationEffectInference.LABEL_FUNCS)[metric][0]
    ramp = [step / 1000 for step in range(-200, 2201)]

    def switches(fn) -> list[int]:
        return [
            index
            for index in range(1, len(ramp))
            if fn(ramp[index]) != fn(ramp[index - 1])
        ]

    ours = switches(lambda value: b2b_interpretation(metric, value))
    theirs = switches(label)
    assert ours, f"{metric} reported one band across its whole range"
    assert ours == theirs


def test_b2b_interpretation_of_an_unknown_metric_is_not_a_crash() -> None:
    assert b2b_interpretation("not-a-metric", 0.5) is None


def test_b2b_interpretation_of_a_missing_value_is_not_a_crash() -> None:
    metric = next(iter(B2B_THRESHOLDS))
    assert b2b_interpretation(metric, None) is None
    assert b2b_interpretation(metric, "n/a") is None


def _triangle_graph():
    import networkx as nx

    graph = nx.Graph()
    for position, residue in ((10, "SER"), (11, "ALA"), (12, "TYR")):
        graph.add_node(position, resname=residue)
    graph.add_edges_from([(10, 11), (11, 12)])
    return graph


def test_rin_to_pyvis_html_colours_by_a_bio2byte_property(tmp_path: Path) -> None:
    """The overlay path must run against a real frame, not just a truthy stub."""
    frame = pd.DataFrame({"Position": [10, 11, 12], "backbone": [0.10, 0.55, 0.95]})
    out = StructureOps.rin_to_pyvis_html(
        _triangle_graph(),
        tmp_path / "rin.html",
        [10],
        [12],
        b2b_frame=frame,
        b2b_metric="backbone",
    )
    html = out.read_text(encoding="utf-8")
    # Three distinct property values must produce three distinct fills.
    fills = {pseudocolor(0.10, 0.95, value) for value in (0.10, 0.55, 0.95)}
    assert len(fills) == 3
    assert sum(fill.lower() in html.lower() for fill in fills) >= 2


def test_rin_to_pyvis_html_matches_values_on_position_not_row_order(tmp_path: Path) -> None:
    """A chain starting at residue 705 must not read row 0's value."""
    shifted = pd.DataFrame({"Position": [12, 11, 10], "backbone": [0.95, 0.55, 0.10]})
    aligned = pd.DataFrame({"Position": [10, 11, 12], "backbone": [0.10, 0.55, 0.95]})
    first = StructureOps.rin_to_pyvis_html(
        _triangle_graph(), tmp_path / "a.html", b2b_frame=shifted, b2b_metric="backbone"
    ).read_text(encoding="utf-8")
    second = StructureOps.rin_to_pyvis_html(
        _triangle_graph(), tmp_path / "b.html", b2b_frame=aligned, b2b_metric="backbone"
    ).read_text(encoding="utf-8")
    # Same position-to-value mapping, only the row order differs -> same colours.
    assert _fill_colours(first) == _fill_colours(second)


def test_rin_to_pyvis_html_falls_back_to_site_colours(tmp_path: Path) -> None:
    """No metric, an absent column or an empty frame must all still render."""
    frame = pd.DataFrame({"Position": [10], "backbone": [0.5]})
    for kwargs in (
        {},
        {"b2b_frame": frame, "b2b_metric": ""},
        {"b2b_frame": frame, "b2b_metric": "absent"},
        {"b2b_frame": pd.DataFrame(), "b2b_metric": "backbone"},
        {"b2b_frame": None, "b2b_metric": "backbone"},
    ):
        out = StructureOps.rin_to_pyvis_html(
            _triangle_graph(), tmp_path / "fallback.html", [10], [12], **kwargs
        )
        assert out.exists() and out.stat().st_size > 0


def test_a_flaky_connect_is_retried_once(monkeypatch, tmp_path: Path) -> None:
    """The observed UniProt failure is intermittent, not an outage.

    A connect times out while a direct request moments earlier succeeded in well under a
    second, so one retry recovers the run instead of leaving the structure pickers empty.
    """
    attempts = {"n": 0}

    def flaky_get(url, headers=None, timeout=None):  # noqa: ANN001, ARG001
        attempts["n"] += 1
        if attempts["n"] == 1:
            raise requests.exceptions.ConnectTimeout("connect timed out")
        return _JsonResponse({"uniProtKBCrossReferences": [_xref("2IVS", "A=705-1013")]})

    monkeypatch.setattr("common.http_lookup.requests.get", flaky_get)
    refs = StructureVizService(tmp_path).fetch_pdb_xrefs("P07949")
    assert attempts["n"] == 2
    assert [ref.pdb_id for ref in refs] == ["2IVS"]


def test_a_persistent_failure_still_raises(monkeypatch, tmp_path: Path) -> None:
    """Retrying must not turn a real outage into a silent empty result."""

    def always_fails(url, headers=None, timeout=None):  # noqa: ANN001, ARG001
        raise requests.exceptions.ConnectTimeout("connect timed out")

    monkeypatch.setattr("common.http_lookup.requests.get", always_fails)
    with pytest.raises(requests.exceptions.ConnectTimeout):
        StructureVizService(tmp_path).fetch_pdb_xrefs("P07949")


def test_it_does_not_retry_forever(monkeypatch, tmp_path: Path) -> None:
    attempts = {"n": 0}

    def always_fails(url, headers=None, timeout=None):  # noqa: ANN001, ARG001
        attempts["n"] += 1
        raise requests.exceptions.ConnectTimeout("connect timed out")

    monkeypatch.setattr("common.http_lookup.requests.get", always_fails)
    with pytest.raises(requests.exceptions.ConnectTimeout):
        StructureVizService(tmp_path).fetch_pdb_xrefs("P07949")
    assert attempts["n"] == 2, "exactly one retry, not an unbounded loop"


def test_a_stalled_host_is_bounded_by_the_read_timeout(monkeypatch, tmp_path: Path) -> None:
    """A host can accept the connection and then never answer.

    Observed against the EBI Proteins API: connect in 0.06s, then silence. Only the read
    timeout bounds that, and these lookups run inside a synchronous reactive effect, so
    the bound is how long one slow upstream freezes every connected session.
    """
    seen: dict = {}

    def stalling_get(url, headers=None, params=None, timeout=None):  # noqa: ANN001, ARG001
        seen["timeout"] = timeout
        raise requests.exceptions.ReadTimeout("host accepted then stalled")

    monkeypatch.setattr("common.http_lookup.requests.get", stalling_get)
    service = StructureVizService(tmp_path)
    with pytest.raises(requests.exceptions.ReadTimeout):
        service.fetch_uniprot_ptms("P07949")
    assert seen["timeout"][1] == StructureVizService.LOOKUP_READ_TIMEOUT
    assert seen["timeout"][1] < service.timeout


def test_file_downloads_keep_a_longer_timeout(monkeypatch, tmp_path: Path) -> None:
    """A large structure genuinely takes time; it must not inherit the lookup bound."""
    seen: dict = {}

    class _Bytes:
        content = b"ATOM\n"
        ok = True
        text = "ATOM\n"

        def raise_for_status(self) -> None:
            return None

    def capture(url, timeout=None, **kwargs):  # noqa: ANN001, ARG001
        seen["timeout"] = timeout
        return _Bytes()

    monkeypatch.setattr("common.http_lookup.requests.get", capture)
    service = StructureVizService(tmp_path)
    service.download_pdb("2IVT")
    read = seen["timeout"][1] if isinstance(seen["timeout"], tuple) else seen["timeout"]
    assert read > StructureVizService.LOOKUP_READ_TIMEOUT
