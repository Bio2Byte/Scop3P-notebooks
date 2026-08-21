from __future__ import annotations

from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pytest
import requests
from Bio.PDB import PDBParser

from common.rinalign import (
    MAX_ALIGNMENT_NODES,
    RESIDUE_GROUPS,
    THREE_TO_ONE,
    RINAlignService,
    StructureEntry,
    align_rins,
    build_rin,
    diff_rins,
    get_cb,
    restype_sim,
    rin_html,
    wl_sigs,
)

FIXTURES = Path(__file__).resolve().parents[2] / "notebooks" / "topology_viewer" / "fixtures"


# ---------------------------------------------------------------------------
# HTTP stubs
# ---------------------------------------------------------------------------


class _Response:
    def __init__(self, payload: Any = None, status_code: int = 200, text: str = "",
                 content_type: str = "application/json") -> None:
        self._payload = payload
        self.status_code = status_code
        self.text = text
        self.headers = {"content-type": content_type}

    def json(self) -> Any:
        if self._payload is None:
            raise ValueError("no json body")
        return self._payload

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")


def _patch_get(monkeypatch: pytest.MonkeyPatch, handler) -> list[dict[str, Any]]:
    """Replace requests.get and record every call, so timeouts can be asserted."""
    calls: list[dict[str, Any]] = []

    def fake_get(url, **kwargs):  # noqa: ANN001
        calls.append({"url": url, **kwargs})
        return handler(url, **kwargs)

    monkeypatch.setattr("common.rinalign.requests.get", fake_get)
    return calls


@pytest.fixture()
def service(tmp_path: Path) -> RINAlignService:
    return RINAlignService(workdir=tmp_path / "work", timeout=7)


# ---------------------------------------------------------------------------
# fetch_uniprot_info
# ---------------------------------------------------------------------------


def test_fetch_uniprot_info_extracts_every_field(service, monkeypatch) -> None:
    payload = {
        "primaryAccession": "P04637",
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Cellular tumor antigen p53"}}},
        "genes": [{"geneName": {"value": "TP53"}}],
        "organism": {"scientificName": "Homo sapiens"},
        "sequence": {"length": 393},
    }
    _patch_get(monkeypatch, lambda url, **kw: _Response(payload))
    info = service.fetch_uniprot_info("P04637")
    assert info["accession"] == "P04637"
    assert info["protein_name"] == "Cellular tumor antigen p53"
    assert info["gene_name"] == "TP53"
    assert info["organism"] == "Homo sapiens"
    assert info["length"] == 393
    assert info["_raw"] is payload


def test_fetch_uniprot_info_falls_back_to_submitted_name(service, monkeypatch) -> None:
    """``recommendedName`` is absent for unreviewed entries."""
    payload = {"proteinDescription": {"submittedName": [{"fullName": {"value": "Uncharacterized"}}]}}
    _patch_get(monkeypatch, lambda url, **kw: _Response(payload))
    assert service.fetch_uniprot_info("X99999")["protein_name"] == "Uncharacterized"


def test_fetch_uniprot_info_accepts_a_bare_string_full_name(service, monkeypatch) -> None:
    payload = {"proteinDescription": {"recommendedName": {"fullName": "Plain string"}}}
    _patch_get(monkeypatch, lambda url, **kw: _Response(payload))
    assert service.fetch_uniprot_info("X1")["protein_name"] == "Plain string"


def test_fetch_uniprot_info_raises_on_a_missing_accession(service, monkeypatch) -> None:
    _patch_get(monkeypatch, lambda url, **kw: _Response(status_code=404))
    with pytest.raises(ValueError, match="404"):
        service.fetch_uniprot_info("NOPE")


def test_metadata_lookups_use_the_shared_bounded_policy(service, monkeypatch) -> None:
    """A hung upstream must not pin the worker forever; the notebook had no timeouts.

    These lookups now share one policy with every other protocol -- a short connect bound
    and a read bound well under what a file download needs -- so a protocol cannot quietly
    opt out of it. Passing ``timeout=self.timeout`` at a call site would silently reinstate
    the long read timeout, which is exactly what this pins against.
    """
    from common import http_lookup as policy

    calls = _patch_get(monkeypatch, lambda url, **kw: _Response({"sequence": {"length": 1}}))
    service.fetch_uniprot_info("P04637")
    assert calls[0]["timeout"] == (policy.CONNECT_TIMEOUT, policy.READ_TIMEOUT)
    assert policy.CONNECT_TIMEOUT < policy.READ_TIMEOUT < service.timeout * 10


def test_a_flaky_lookup_is_retried_once(service, monkeypatch) -> None:
    """Observed for real: UniProt dropping a TLS handshake mid-run.

    A direct request moments earlier succeeded, so one retry recovers the run instead of
    leaving the structure dropdown empty -- which reads to the user as "this protein has
    no structures".
    """
    attempts = {"n": 0}

    def flaky(url, **kwargs):  # noqa: ANN001, ARG001
        attempts["n"] += 1
        if attempts["n"] == 1:
            raise requests.exceptions.SSLError("UNEXPECTED_EOF_WHILE_READING")
        return _Response({"sequence": {"length": 1}})

    _patch_get(monkeypatch, flaky)
    service.fetch_uniprot_info("P04637")
    assert attempts["n"] == 2


def test_a_persistent_lookup_failure_still_raises(service, monkeypatch) -> None:
    """Retrying must never turn an outage into a silently empty dropdown."""

    def always_fails(url, **kwargs):  # noqa: ANN001, ARG001
        raise requests.exceptions.SSLError("handshake failed")

    _patch_get(monkeypatch, always_fails)
    with pytest.raises(requests.exceptions.SSLError):
        service.fetch_uniprot_info("P04637")


# ---------------------------------------------------------------------------
# fetch_pdb_structures -- strategy 1, the intricate parser
# ---------------------------------------------------------------------------


def _xref(pdb_id: str, chains: str, method: str = "X-ray", resolution: str = "2.20 A") -> dict:
    return {
        "database": "PDB",
        "id": pdb_id,
        "properties": [
            {"key": "Method", "value": method},
            {"key": "Resolution", "value": resolution},
            {"key": "Chains", "value": chains},
        ],
    }


def test_uniprot_xrefs_expand_slash_and_comma_chain_groups(service) -> None:
    data = {"uniProtKBCrossReferences": [_xref("1abc", "A/B=1-100, C=10-200")]}
    entries = service.fetch_pdb_structures("P1", uniprot_data=data)
    assert [(e.pdb_id, e.chain_id, e.unp_start, e.unp_end) for e in entries] == [
        ("1ABC", "A", 1, 100),
        ("1ABC", "B", 1, 100),
        ("1ABC", "C", 10, 200),
    ]
    assert all(e.resolution == 2.2 for e in entries)
    assert all(e.source == "PDB" for e in entries)


def test_uniprot_xrefs_deduplicate_repeated_entry_chain_pairs(service) -> None:
    data = {"uniProtKBCrossReferences": [_xref("1ABC", "A=1-50"), _xref("1abc", "A=1-50")]}
    assert len(service.fetch_pdb_structures("P1", uniprot_data=data)) == 1


def test_unparseable_resolution_becomes_none_not_an_error(service) -> None:
    data = {"uniProtKBCrossReferences": [_xref("1ABC", "A=1-50", resolution="-")]}
    entry = service.fetch_pdb_structures("P1", uniprot_data=data)[0]
    assert entry.resolution is None
    assert "X-ray" in entry.label  # method is used when resolution is unavailable


def test_missing_residue_range_leaves_the_bounds_unknown(service) -> None:
    data = {"uniProtKBCrossReferences": [_xref("1ABC", "A")]}
    entry = service.fetch_pdb_structures("P1", uniprot_data=data)[0]
    assert (entry.unp_start, entry.unp_end) == ("?", "?")


def test_entries_sort_by_resolution_with_unknowns_last(service) -> None:
    data = {
        "uniProtKBCrossReferences": [
            _xref("3AAA", "A=1-9", resolution="3.00 A"),
            _xref("1AAA", "A=1-9", resolution="-", method="NMR"),
            _xref("2AAA", "A=1-9", resolution="1.50 A"),
        ]
    }
    entries = service.fetch_pdb_structures("P1", uniprot_data=data)
    assert [e.pdb_id for e in entries] == ["2AAA", "3AAA", "1AAA"]


def test_non_pdb_cross_references_are_ignored(service) -> None:
    data = {"uniProtKBCrossReferences": [{"database": "AlphaFoldDB", "id": "P1"}]}
    assert service.fetch_pdb_structures("P1", uniprot_data=data) == []


# ---------------------------------------------------------------------------
# fetch_pdb_structures -- fallback strategies
# ---------------------------------------------------------------------------


def test_best_structures_fallback_handles_case_variant_keys(service, monkeypatch) -> None:
    payload = {"p04637": [{"pdb_id": "1tup", "chain_id": "A", "resolution": 2.2,
                           "unp_start": 94, "unp_end": 312}]}

    def handler(url, **kw):
        return _Response(payload) if "best_structures" in url else _Response(status_code=404)

    _patch_get(monkeypatch, handler)
    entries = service.fetch_pdb_structures("P04637", uniprot_data=None)
    assert len(entries) == 1
    assert (entries[0].pdb_id, entries[0].chain_id, entries[0].resolution) == ("1TUP", "A", 2.2)


def test_sifts_fallback_reads_nested_residue_numbers(service, monkeypatch) -> None:
    payload = {
        "P1": {
            "PDB": {
                "1abc": [
                    {"chain_id": "A", "start": {"residue_number": 5},
                     "end": {"residue_number": 55}, "resolution": 1.9}
                ]
            }
        }
    }

    def handler(url, **kw):
        if "best_structures" in url:
            return _Response(status_code=404)
        return _Response(payload)

    _patch_get(monkeypatch, handler)
    entries = service.fetch_pdb_structures("P1", uniprot_data=None)
    assert len(entries) == 1
    assert (entries[0].pdb_id, entries[0].unp_start, entries[0].unp_end) == ("1ABC", 5, 55)


def test_a_network_outage_is_logged_and_yields_no_structures(service, monkeypatch, caplog) -> None:
    """The notebook's bare ``except: pass`` made an outage look like 'no PDB entries'."""
    def handler(url, **kw):
        raise ConnectionError("dns failure")

    _patch_get(monkeypatch, handler)
    with caplog.at_level("WARNING"):
        assert service.fetch_pdb_structures("P1", uniprot_data=None) == []
    assert any("failed" in record.message for record in caplog.records)


# ---------------------------------------------------------------------------
# check_alphafold / download_structure
# ---------------------------------------------------------------------------


def test_check_alphafold_builds_an_entry(service, monkeypatch) -> None:
    payload = [{"pdbUrl": "https://af/AF-P1-F1.pdb", "cifUrl": "https://af/AF-P1-F1.cif"}]
    _patch_get(monkeypatch, lambda url, **kw: _Response(payload))
    entry = service.check_alphafold("P1")
    assert entry is not None
    assert (entry.pdb_id, entry.chain_id, entry.source) == ("AF-P1", "A", "AlphaFold")


@pytest.mark.parametrize("payload,status", [(None, 404), ([], 200)])
def test_check_alphafold_returns_none_when_there_is_no_model(service, monkeypatch, payload, status) -> None:
    _patch_get(monkeypatch, lambda url, **kw: _Response(payload, status_code=status))
    assert service.check_alphafold("P1") is None


def test_download_falls_back_from_cif_to_pdb(service, monkeypatch) -> None:
    def handler(url, **kw):
        if url.endswith(".cif"):
            return _Response(status_code=404)
        return _Response(text="ATOM      1  CA  ALA A   1\n")

    calls = _patch_get(monkeypatch, handler)
    path, text = service.download_structure(
        StructureEntry(pdb_id="1ABC", chain_id="A", label="l", source="PDB")
    )
    assert path.suffix == ".pdb"
    assert text.startswith("ATOM")
    assert [c["url"].rsplit(".", 1)[-1] for c in calls] == ["cif", "pdb"]


def test_download_is_cached_so_repeat_clicks_do_not_refetch(service, monkeypatch) -> None:
    """The notebook mkdtemp'd per call and leaked a directory on every Generate."""
    calls = _patch_get(monkeypatch, lambda url, **kw: _Response(text="ATOM\n"))
    entry = StructureEntry(pdb_id="1ABC", chain_id="A", label="l", source="PDB")
    first, _ = service.download_structure(entry)
    second, _ = service.download_structure(entry)
    assert first == second
    assert len(calls) == 1


def test_download_prefers_the_alphafold_pdb_url(service, monkeypatch) -> None:
    _patch_get(monkeypatch, lambda url, **kw: _Response(text="ATOM\n"))
    entry = StructureEntry(
        pdb_id="AF-P1", chain_id="A", label="l", source="AlphaFold",
        pdb_url="https://af/AF-P1-F1.pdb", cif_url="https://af/AF-P1-F1.cif",
    )
    path, _ = service.download_structure(entry)
    assert path.suffix == ".pdb"


def test_download_reports_the_failure_reason(service, monkeypatch) -> None:
    _patch_get(monkeypatch, lambda url, **kw: _Response(status_code=503))
    with pytest.raises(ValueError, match="503"):
        service.download_structure(
            StructureEntry(pdb_id="1ABC", chain_id="A", label="l", source="PDB")
        )


# ---------------------------------------------------------------------------
# PTM / variant positions
# ---------------------------------------------------------------------------


def test_scop3p_positions_are_ints(service, monkeypatch) -> None:
    """PTM positions now come through Scop3PClient, so the v1 field name applies and
    the request goes to common.services rather than being built here."""
    payload = [
        {"uniprot_position": 15},
        {"uniprot_position": "27"},
        {"uniprot_position": None},
        {},
    ]
    monkeypatch.setattr("common.services.requests.get", lambda *a, **kw: _Response(payload))
    assert service.fetch_scop3p_ptm_positions("P1") == {15, 27}


def test_scop3p_positions_are_empty_for_an_uncovered_protein(service, monkeypatch) -> None:
    monkeypatch.setattr("common.services.requests.get", lambda *a, **kw: _Response([]))
    assert service.fetch_scop3p_ptm_positions("P0DTD1") == set()


def test_a_broken_scop3p_endpoint_is_reported_not_silently_empty(service, monkeypatch) -> None:
    """The bug this replaced: the pre-v1 URL returned the app's HTML with a 200, the
    broad except swallowed the JSON error, and every protein reported zero PTMs. A
    transport or endpoint failure must now surface."""
    from common.services import Scop3PApiError

    monkeypatch.setattr(
        "common.services.requests.get",
        lambda *a, **kw: _Response(None, content_type="text/html"),
    )
    with pytest.raises(Scop3PApiError):
        service.fetch_scop3p_ptm_positions("P07949")


def test_uniprot_ptm_positions_keep_only_single_residue_features(service, monkeypatch) -> None:
    payload = {
        "features": [
            {"category": "PTM", "begin": "10", "end": "10"},
            {"category": "PTM", "begin": "20", "end": "25"},   # a range: no single site
            {"category": "DOMAIN", "begin": "30", "end": "30"},  # wrong category
        ]
    }
    _patch_get(monkeypatch, lambda url, **kw: _Response(payload))
    assert service.fetch_uniprot_ptm_positions("P1") == {10}


def test_variant_positions_filter_on_disease_association(service, monkeypatch) -> None:
    payload = {
        "features": [
            {"type": "VARIANT", "begin": "10", "association": [{"disease": True}]},
            {"type": "VARIANT", "begin": "20", "association": [{"disease": False}]},
            {"type": "VARIANT", "begin": "30"},
            {"type": "MUTAGEN", "begin": "40", "association": [{"disease": True}]},
        ]
    }
    _patch_get(monkeypatch, lambda url, **kw: _Response(payload))
    assert service.fetch_uniprot_variant_positions("P1", disease_only=True) == {10}
    assert service.fetch_uniprot_variant_positions("P1", disease_only=False) == {10, 20, 30}


@pytest.mark.parametrize(
    "method",
    ["fetch_uniprot_ptm_positions", "fetch_uniprot_variant_positions"],
)
def test_uniprot_annotation_fetchers_degrade_to_an_empty_set(service, monkeypatch, method) -> None:
    """A UniProt annotation source being down must not block the comparison.

    Scop3P is deliberately not in this list: it is the primary PTM source, so a broken
    endpoint there is reported rather than quietly reduced to "no sites". See
    test_a_broken_scop3p_endpoint_is_reported_not_silently_empty.
    """
    _patch_get(monkeypatch, lambda url, **kw: (_ for _ in ()).throw(ConnectionError("down")))
    assert getattr(service, method)("P1") == set()


# ---------------------------------------------------------------------------
# Residue tables
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("residue", ["MSE", "PTR", "SEP", "TPO"])
def test_modified_residues_become_network_nodes(residue: str) -> None:
    """Without these, phospho-residues are dropped as heteroatoms and the network
    loses exactly the positions the whole tool is about."""
    assert residue in THREE_TO_ONE


def test_phospho_residues_share_a_group_with_their_parents() -> None:
    assert RESIDUE_GROUPS["SEP"] == RESIDUE_GROUPS["SER"]
    assert RESIDUE_GROUPS["TPO"] == RESIDUE_GROUPS["THR"]
    assert RESIDUE_GROUPS["PTR"] == RESIDUE_GROUPS["TYR"]


def test_selenomethionine_has_no_similarity_group_yet() -> None:
    """Documents a known gap carried over from the notebook, deliberately unchanged.

    MSE is in THREE_TO_ONE (mapped to 'M', so it becomes a network node) but was
    never added to RESIDUE_GROUPS alongside PTR/SEP/TPO. The consequence is confined
    to cross-protein alignment: restype_sim gives MSE group '?', so a
    selenomethionine scores 0.0 against a methionine instead of 1.0/0.5.
    Selenomethionine is common in crystal structures, so this is worth a decision --
    but changing it changes published alignment scores, so parity is kept here and
    the finding is recorded in docs/use-cases/rinalign.md instead.
    """
    assert "MSE" not in RESIDUE_GROUPS
    assert restype_sim("MSE", "MET") == 0.0


# ---------------------------------------------------------------------------
# build_rin
# ---------------------------------------------------------------------------


def _write_pdb(path: Path, rows: list[tuple[int, str, str, float, float, float]]) -> Path:
    lines = []
    for serial, (seq, resname, atom, x, y, z) in enumerate(rows, start=1):
        lines.append(
            f"ATOM  {serial:>5}  {atom:<3} {resname:>3} A{seq:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 50.00           C"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


def test_build_rin_node_attributes_and_ids(tmp_path: Path) -> None:
    path = _write_pdb(tmp_path / "t.pdb", [
        (1, "ALA", "CA", 0.0, 0.0, 0.0), (1, "ALA", "CB", 0.5, 0.0, 0.0),
        (5, "GLY", "CA", 7.0, 0.0, 0.0),
        (9, "MSE", "CA", 20.0, 0.0, 0.0), (9, "MSE", "CB", 20.5, 0.0, 0.0),
    ])
    graph, residues = build_rin(path, cutoff=8.0)
    assert set(graph.nodes) == {"ALA_1", "GLY_5", "MSE_9"}
    assert graph.nodes["ALA_1"]["position"] == 1
    assert graph.nodes["ALA_1"]["resname"] == "ALA"
    assert graph.nodes["MSE_9"]["one_letter"] == "M"
    assert graph.nodes["GLY_5"]["chain"] == "A"
    assert set(residues) == set(graph.nodes)


def test_glycine_falls_back_to_ca_when_cb_is_absent(tmp_path: Path) -> None:
    path = _write_pdb(tmp_path / "g.pdb", [(1, "GLY", "CA", 0.0, 0.0, 0.0)])
    graph, residues = build_rin(path, cutoff=8.0)
    assert "GLY_1" in graph
    assert np.allclose(residues["GLY_1"]["cb"], [0.0, 0.0, 0.0])


def test_sequence_adjacent_pairs_in_one_chain_are_excluded(tmp_path: Path) -> None:
    """Backbone neighbours are always in contact, so they carry no information."""
    path = _write_pdb(tmp_path / "adj.pdb", [
        (1, "ALA", "CB", 0.0, 0.0, 0.0),
        (2, "ALA", "CB", 1.0, 0.0, 0.0),   # adjacent -> excluded
        (3, "ALA", "CB", 2.0, 0.0, 0.0),   # |1-3| == 2 -> kept
    ])
    graph, _ = build_rin(path, cutoff=8.0)
    edges = {tuple(sorted(e)) for e in graph.edges}
    assert ("ALA_1", "ALA_2") not in edges
    assert ("ALA_2", "ALA_3") not in edges
    assert ("ALA_1", "ALA_3") in edges


def test_cutoff_is_inclusive_and_distances_are_rounded(tmp_path: Path) -> None:
    path = _write_pdb(tmp_path / "c.pdb", [
        (1, "ALA", "CB", 0.0, 0.0, 0.0),
        (4, "ALA", "CB", 7.0, 0.0, 0.0),
    ])
    assert build_rin(path, cutoff=8.0)[0].number_of_edges() == 1
    assert build_rin(path, cutoff=6.0)[0].number_of_edges() == 0
    graph, _ = build_rin(path, cutoff=7.0)  # exactly at the cutoff
    assert graph.number_of_edges() == 1
    assert graph.edges[("ALA_1", "ALA_4")]["distance"] == 7.0


def test_chain_filter_keeps_only_the_requested_chain(tmp_path: Path) -> None:
    path = tmp_path / "two.pdb"
    path.write_text(
        "ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 50.00           C\n"
        "ATOM      2  CB  ALA B   1       1.000   0.000   0.000  1.00 50.00           C\n"
        "END\n"
    )
    graph, _ = build_rin(path, chain_id="A", cutoff=8.0)
    assert list(graph.nodes) == ["ALA_1"]
    assert graph.nodes["ALA_1"]["chain"] == "A"


def test_unknown_residues_are_skipped(tmp_path: Path) -> None:
    path = _write_pdb(tmp_path / "u.pdb", [
        (1, "ALA", "CB", 0.0, 0.0, 0.0),
        (2, "XYZ", "CB", 1.0, 0.0, 0.0),
    ])
    assert list(build_rin(path, cutoff=8.0)[0].nodes) == ["ALA_1"]


def _naive_edges(path: Path, cutoff: float) -> set[tuple[str, str]]:
    """The notebook's O(n^2) contact loop, kept as the reference implementation."""
    structure = PDBParser(QUIET=True).get_structure("p", str(path))
    residues = []
    for chain in structure[0]:
        for residue in chain:
            hetero, seq, _ = residue.get_id()
            name = residue.get_resname().strip()
            if name not in THREE_TO_ONE or (hetero != " " and name not in THREE_TO_ONE):
                continue
            coordinates = get_cb(residue)
            if coordinates is None:
                continue
            residues.append((f"{name}_{seq}", seq, chain.get_id(), coordinates))
    edges = set()
    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            if residues[i][2] == residues[j][2] and abs(residues[i][1] - residues[j][1]) <= 1:
                continue
            if np.linalg.norm(residues[i][3] - residues[j][3]) <= cutoff:
                edges.add(tuple(sorted((residues[i][0], residues[j][0]))))
    return edges


@pytest.mark.parametrize("name", ["annotated.pdb", "bare.pdb", "big.pdb"])
@pytest.mark.parametrize("cutoff", [6.0, 8.0, 11.5])
def test_kdtree_contacts_match_the_naive_loop_exactly(name: str, cutoff: float) -> None:
    """The KD-tree rewrite must not change the science.

    ``KDTree.query_pairs(r)`` selects pairs at ``distance <= r``, the same comparison
    the notebook made, so the edge sets are identical. ``annotated.cif`` is excluded:
    the synthetic fixture omits ``_atom_site.label_alt_id``, which BioPython's
    MMCIFParser requires and every real PDBe or RCSB file provides.
    """
    path = FIXTURES / name
    graph, _ = build_rin(path, cutoff=cutoff)
    assert {tuple(sorted(e)) for e in graph.edges} == _naive_edges(path, cutoff)


def test_rin_html_reports_the_graph_size() -> None:
    graph = nx.Graph([("A_1", "A_3"), ("A_3", "A_5")])
    html = rin_html(graph, "1ABC A")
    assert "1ABC A" in html and ">3<" in html and ">2<" in html


# ---------------------------------------------------------------------------
# diff_rins
# ---------------------------------------------------------------------------


def _positional_graph(spec: dict[int, str], edges: list[tuple[int, int]]) -> nx.Graph:
    graph = nx.Graph()
    names = {position: f"{resname}_{position}" for position, resname in spec.items()}
    for position, resname in spec.items():
        graph.add_node(names[position], position=position, resname=resname)
    for left, right in edges:
        graph.add_edge(names[left], names[right])
    return graph


def test_diff_splits_edges_into_conserved_lost_and_gained() -> None:
    left = _positional_graph({1: "ALA", 2: "SER", 3: "LEU", 4: "GLY"},
                             [(1, 3), (1, 4), (2, 4)])
    right = _positional_graph({1: "ALA", 2: "SER", 3: "LEU", 4: "GLY"},
                              [(1, 3), (2, 3)])
    result = diff_rins(left, right)
    assert result["conserved"] == [(1, 3)]
    assert result["lost"] == [(1, 4), (2, 4)]
    assert result["gained"] == [(2, 3)]
    assert result["jaccard"] == pytest.approx(1 / 4)


def test_diff_reports_residue_identity_changes_as_mutations() -> None:
    left = _positional_graph({1: "ALA", 2: "ARG"}, [])
    right = _positional_graph({1: "ALA", 2: "TRP"}, [])
    assert diff_rins(left, right)["mutations"] == [
        {"position": 2, "left": "ARG", "right": "TRP"}
    ]


def test_diff_compares_only_positions_present_in_both_structures() -> None:
    """A contact is only 'lost' if the residue exists in both models; otherwise it is
    merely unresolved, and reporting it as lost would be wrong."""
    left = _positional_graph({1: "ALA", 5: "SER", 9: "LEU"}, [(1, 5), (1, 9)])
    right = _positional_graph({1: "ALA", 5: "SER"}, [(1, 5)])
    result = diff_rins(left, right)
    assert result["matched_pos"] == [1, 5]
    assert result["only_left_pos"] == [9]
    assert result["only_right_pos"] == []
    assert result["lost"] == []          # (1, 9) is unmatched, not lost
    assert result["conserved"] == [(1, 5)]
    assert (1, 9) in result["onlyA_edges"]


def test_residue_impact_is_sorted_by_absolute_net_change() -> None:
    left = _positional_graph({1: "ALA", 2: "SER", 3: "LEU", 4: "GLY", 5: "VAL"},
                             [(1, 3), (1, 4), (1, 5)])
    right = _positional_graph({1: "ALA", 2: "SER", 3: "LEU", 4: "GLY", 5: "VAL"},
                              [(2, 4)])
    impact = diff_rins(left, right)["residue_impact"]
    magnitudes = [abs(row["net_change"]) for row in impact]
    assert magnitudes == sorted(magnitudes, reverse=True)
    by_position = {row["position"]: row for row in impact}
    assert by_position[1]["lost"] == 3 and by_position[1]["net_change"] == -3
    assert by_position[2]["gained"] == 1 and by_position[2]["net_change"] == 1


def test_impact_rows_flag_mutated_positions() -> None:
    left = _positional_graph({1: "ALA", 2: "ARG"}, [])
    right = _positional_graph({1: "ALA", 2: "TRP"}, [])
    flags = {row["position"]: row["is_mutation"] for row in diff_rins(left, right)["residue_impact"]}
    assert flags == {1: False, 2: True}


def test_identical_networks_score_a_perfect_jaccard() -> None:
    graph = _positional_graph({1: "ALA", 2: "SER", 3: "LEU"}, [(1, 3)])
    result = diff_rins(graph, graph.copy())
    assert result["jaccard"] == 1.0
    assert result["lost"] == [] and result["gained"] == []


def test_two_edgeless_networks_are_treated_as_identical() -> None:
    graph = _positional_graph({1: "ALA"}, [])
    assert diff_rins(graph, graph.copy())["jaccard"] == 1.0


# ---------------------------------------------------------------------------
# restype_sim / wl_sigs / align_rins
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("first,second,expected", [
    ("ALA", "ALA", 1.0),
    ("ALA", "LEU", 0.5),      # both hydrophobic
    ("ASP", "GLU", 0.5),      # both negative
    ("ALA", "ASP", 0.0),
    ("XXX", "YYY", 0.0),      # both unknown: must not score as a match
    ("XXX", "XXX", 1.0),      # identical strings still score 1.0
])
def test_restype_sim(first: str, second: str, expected: float) -> None:
    assert restype_sim(first, second) == expected


def test_wl_sigs_returns_depth_plus_one_labels_per_node() -> None:
    graph = _positional_graph({1: "ALA", 2: "SER", 3: "LEU"}, [(1, 3), (2, 3)])
    signatures = wl_sigs(graph, k=3)
    assert set(signatures) == set(graph.nodes)
    assert all(len(labels) == 4 for labels in signatures.values())


def test_wl_sigs_agree_on_isomorphic_graphs() -> None:
    first = _positional_graph({1: "ALA", 3: "ALA", 5: "ALA", 7: "ALA"},
                              [(1, 3), (3, 5), (5, 7)])
    second = _positional_graph({2: "ALA", 4: "ALA", 6: "ALA", 8: "ALA"},
                               [(2, 4), (4, 6), (6, 8)])
    assert sorted(map(tuple, wl_sigs(first, 3).values())) == sorted(
        map(tuple, wl_sigs(second, 3).values())
    )


def test_align_isomorphic_graphs_reaches_a_perfect_score() -> None:
    first = _positional_graph({1: "ALA", 3: "SER", 5: "LEU", 7: "GLY"},
                              [(1, 3), (3, 5), (5, 7)])
    second = _positional_graph({2: "ALA", 4: "SER", 6: "LEU", 8: "GLY"},
                               [(2, 4), (4, 6), (6, 8)])
    result = align_rins(first, second)
    assert len(result["mapping"]) == 4
    assert result["jaccard"] == 1.0
    assert result["only_G1"] == [] and result["only_G2"] == []


def test_align_mapping_is_one_to_one() -> None:
    first = _positional_graph({1: "ALA", 3: "SER", 5: "LEU"}, [(1, 3)])
    second = _positional_graph({2: "ALA", 4: "SER", 6: "TRP"}, [(2, 4)])
    mapping = align_rins(first, second)["mapping"]
    assert len({a for a, _, _ in mapping}) == len(mapping)
    assert len({b for _, b, _ in mapping}) == len(mapping)


def test_align_scores_are_bounded_by_one() -> None:
    first = _positional_graph({1: "ALA", 3: "SER"}, [(1, 3)])
    second = _positional_graph({2: "ALA", 4: "SER"}, [(2, 4)])
    assert all(0.0 <= score <= 1.0 for _, _, score in align_rins(first, second)["mapping"])


def test_align_refuses_a_network_that_would_exhaust_the_worker() -> None:
    """A dense n1 x n2 matrix filled from a Python loop blocks the event loop for
    every session, so an oversized request is rejected with an actionable message."""
    big = nx.Graph()
    for index in range(MAX_ALIGNMENT_NODES + 1):
        big.add_node(f"ALA_{index}", position=index, resname="ALA")
    small = _positional_graph({1: "ALA"}, [])
    with pytest.raises(ValueError, match="limited to"):
        align_rins(big, small)
    with pytest.raises(ValueError, match="limited to"):
        align_rins(small, big)
