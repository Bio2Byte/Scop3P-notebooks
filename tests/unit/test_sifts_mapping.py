"""SIFTS residue numbering: UniProt positions onto a PDB entry's author numbering.

The failure this guards against is silent. A PTM at UniProt 918 drawn on an entry whose
author numbering is offset by one lands on residue 917, renders perfectly, and is wrong.
Nothing raises, so it has to be pinned by test.

Payload shapes here are copied from live PDBe responses (2IVT, 1A3N), including the
awkward one: ``author_residue_number`` is null at one end of most real segments.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from common.structure_viz import (
    MAPPING_SOURCE_LABELS,
    SIFTS_SOURCES,
    PositionMapping,
    StructureVizService,
    identity_mapping,
    offset_mapping,
    parse_pdbe_uniprot_mappings,
    parse_sifts_map_from_cif,
    remap_positions,
    remap_site_rows,
)


def _segment(chain, unp_start, unp_end, author_start, author_end, label_start=1, label_end=1):
    return {
        "chain_id": chain,
        "struct_asym_id": chain,
        "unp_start": unp_start,
        "unp_end": unp_end,
        "start": {"author_residue_number": author_start, "residue_number": label_start},
        "end": {"author_residue_number": author_end, "residue_number": label_end},
    }


def _payload(pdb_id, accession, *segments):
    return {pdb_id.lower(): {"UniProt": {accession: {"mappings": list(segments)}}}}


# --------------------------------------------------------------------------------------
# Tier 1: the PDBe SIFTS API
# --------------------------------------------------------------------------------------


def test_api_mapping_expands_a_segment_residue_by_residue() -> None:
    """1A3N chain A: author 1-141 against UniProt 2-142, the cleaved initiator Met."""
    payload = _payload("1a3n", "P69905", _segment("A", 2, 142, 1, 141))
    mapping = parse_pdbe_uniprot_mappings(payload, "1A3N", "P69905", "A")
    assert len(mapping) == 141
    assert mapping[1] == 2
    assert mapping[141] == 142
    assert {uniprot - author for author, uniprot in mapping.items()} == {1}


def test_api_mapping_reconstructs_a_null_author_end() -> None:
    """PDBe reports a null author number at one end of most real segments.

    2IVT: author start 703, author end null, UniProt 703-1013. The end is recoverable
    from the start plus the UniProt span, because a SIFTS segment is colinear.
    """
    payload = _payload("2ivt", "P07949", _segment("A", 703, 1013, 703, None, 4, 314))
    mapping = parse_pdbe_uniprot_mappings(payload, "2IVT", "P07949", "A")
    assert len(mapping) == 311
    assert mapping[703] == 703
    assert mapping[1013] == 1013


def test_api_mapping_reconstructs_a_null_author_start() -> None:
    """1A3N chain B: author start null, end 146, UniProt 2-147."""
    payload = _payload("1a3n", "P68871", _segment("B", 2, 147, None, 146))
    mapping = parse_pdbe_uniprot_mappings(payload, "1A3N", "P68871", "B")
    assert mapping[1] == 2
    assert mapping[146] == 147


def test_api_mapping_never_substitutes_label_numbering_for_author_numbering() -> None:
    """The regression that would corrupt every position by ~700.

    2IVT's label numbering runs 4-314 while its author numbering runs 703-1013. If a
    missing author number fell back to ``residue_number``, the mapping would be built in
    the wrong coordinate system and every mark would be ~700 residues out.
    """
    payload = _payload("2ivt", "P07949", _segment("A", 703, 1013, 703, None, 4, 314))
    mapping = parse_pdbe_uniprot_mappings(payload, "2IVT", "P07949", "A")
    assert 4 not in mapping, "label numbering leaked into the author-numbered map"
    assert 314 not in mapping
    assert min(mapping) == 703


def test_api_mapping_skips_a_segment_with_no_author_numbers_at_all() -> None:
    payload = _payload("9xyz", "P00001", _segment("A", 1, 50, None, None))
    assert parse_pdbe_uniprot_mappings(payload, "9XYZ", "P00001", "A") == {}


def test_api_mapping_filters_by_chain() -> None:
    payload = _payload(
        "1a3n", "P69905", _segment("A", 2, 4, 1, 3), _segment("C", 2, 4, 501, 503)
    )
    assert set(parse_pdbe_uniprot_mappings(payload, "1A3N", "P69905", "A")) == {1, 2, 3}
    assert set(parse_pdbe_uniprot_mappings(payload, "1A3N", "P69905", "C")) == {501, 502, 503}


def test_api_mapping_accepts_an_isoform_keyed_block() -> None:
    payload = _payload("1abc", "P07949-2", _segment("A", 1, 3, 1, 3))
    assert parse_pdbe_uniprot_mappings(payload, "1ABC", "P07949", "A")


def test_api_mapping_handles_a_descending_segment() -> None:
    """Segments can run backwards, so the step must be signed."""
    payload = _payload("1abc", "P00001", _segment("A", 10, 12, 30, 28))
    mapping = parse_pdbe_uniprot_mappings(payload, "1ABC", "P00001", "A")
    assert mapping == {30: 10, 29: 11, 28: 12}


def test_api_mapping_of_an_unrelated_entry_is_empty() -> None:
    payload = _payload("1abc", "P00001", _segment("A", 1, 3, 1, 3))
    assert parse_pdbe_uniprot_mappings(payload, "9ZZZ", "P00001", "A") == {}


# --------------------------------------------------------------------------------------
# Tier 2: the SIFTS-enriched mmCIF
# --------------------------------------------------------------------------------------

_CIF_WITH_SIFTS = """\
data_2IVT
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_sifts_xref_db_acc
_atom_site.pdbx_sifts_xref_db_num
ATOM CA A 4 A 703 P07949 703
ATOM CB A 4 A 703 P07949 703
ATOM CA A 5 A 704 P07949 704
ATOM CA B 5 B 904 P07949 704
ATOM CA A 6 A 705 ? ?
#
"""

_CIF_WITHOUT_SIFTS = """\
data_1ABC
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
ATOM CA A 1
ATOM CA A 2
#
"""


def test_cif_parser_reads_the_sifts_columns() -> None:
    mapping = parse_sifts_map_from_cif(_CIF_WITH_SIFTS, "P07949", "A")
    assert mapping == {703: 703, 704: 704}


def test_cif_parser_filters_by_chain() -> None:
    assert parse_sifts_map_from_cif(_CIF_WITH_SIFTS, "P07949", "B") == {904: 704}


def test_cif_parser_ignores_unmapped_residues() -> None:
    """A residue with "?" for its UniProt number is absent, not residue zero."""
    assert 705 not in parse_sifts_map_from_cif(_CIF_WITH_SIFTS, "P07949", "A")


def test_cif_parser_returns_nothing_when_sifts_columns_are_absent() -> None:
    """A plain RCSB mmCIF has atom sites but no SIFTS, and must not be mistaken for one."""
    assert parse_sifts_map_from_cif(_CIF_WITHOUT_SIFTS, "P07949", "A") == {}


def test_cif_parser_filters_by_accession() -> None:
    assert parse_sifts_map_from_cif(_CIF_WITH_SIFTS, "P99999", "A") == {}


def test_cif_parser_tolerates_junk() -> None:
    for text in ("", "not a cif at all", "loop_\n_atom_site.foo\n"):
        assert parse_sifts_map_from_cif(text, "P07949", "A") == {}


# --------------------------------------------------------------------------------------
# Tier 3: the offset fallback
# --------------------------------------------------------------------------------------


def test_offset_mapping_lines_up_two_ranges() -> None:
    assert offset_mapping((705, 707), (1, 3)) == {1: 705, 2: 706, 3: 707}


def test_offset_mapping_truncates_to_the_shorter_range() -> None:
    mapping = offset_mapping((1, 100), (1, 10))
    assert len(mapping) == 10


@pytest.mark.parametrize(
    "uniprot_range, pdb_range", [(None, (1, 3)), ((1, 3), None), (None, None)]
)
def test_offset_mapping_needs_both_ranges(uniprot_range, pdb_range) -> None:
    assert offset_mapping(uniprot_range, pdb_range) == {}


# --------------------------------------------------------------------------------------
# Applying a mapping
# --------------------------------------------------------------------------------------


def _mapping(pdb_to_uniprot, source="sifts-api"):
    inverted: dict[int, list[int]] = {}
    for pdb_resi, uniprot in pdb_to_uniprot.items():
        inverted.setdefault(uniprot, []).append(pdb_resi)
    return PositionMapping(pdb_to_uniprot, inverted, source, "1ABC", "A")


def test_remap_positions_translates_into_author_numbering() -> None:
    mapping = _mapping({1: 2, 2: 3, 3: 4})
    assert remap_positions([2, 4], mapping) == [1, 3]


def test_remap_positions_drops_what_the_structure_does_not_contain() -> None:
    """Better a missing mark than one drawn on the wrong residue."""
    assert remap_positions([2, 999], _mapping({1: 2})) == [1]


def test_remap_positions_marks_every_chain_holding_the_residue() -> None:
    """A homodimer has one sequence position in two places, and both should be marked."""
    mapping = _mapping({1: 5, 501: 5})
    assert remap_positions([5], mapping) == [1, 501]


def test_remap_positions_is_a_no_op_for_identity() -> None:
    assert remap_positions([10, 20], identity_mapping()) == [10, 20]


def test_remap_site_rows_keeps_the_row_and_records_the_original() -> None:
    rows = [{"position": 2, "residue": "SER"}]
    out = remap_site_rows(rows, _mapping({1: 2}))
    assert out == [{"position": 1, "residue": "SER", "uniprot_position": 2}]


def test_remap_site_rows_does_not_mutate_its_input() -> None:
    rows = [{"position": 2, "residue": "SER"}]
    remap_site_rows(rows, _mapping({1: 2}))
    assert rows == [{"position": 2, "residue": "SER"}]


def test_remap_site_rows_skips_an_unparseable_position() -> None:
    assert remap_site_rows([{"position": "n/a"}], _mapping({1: 2})) == []


def test_remap_site_rows_passes_identity_through_unchanged() -> None:
    rows = [{"position": 918, "residue": "TYR"}]
    assert remap_site_rows(rows, identity_mapping()) == rows


# --------------------------------------------------------------------------------------
# Provenance, which is what the user is shown
# --------------------------------------------------------------------------------------


def test_every_source_has_a_label() -> None:
    assert SIFTS_SOURCES <= set(MAPPING_SOURCE_LABELS)
    assert set(MAPPING_SOURCE_LABELS) == {
        "sifts-api",
        "sifts-mmcif",
        "chain-range-offset",
        "direct",
    }


@pytest.mark.parametrize("source", sorted(MAPPING_SOURCE_LABELS))
def test_describe_names_the_method_it_used(source: str) -> None:
    mapping = _mapping({1: 2}, source=source)
    described = mapping.describe()
    assert described.endswith(".")
    assert MAPPING_SOURCE_LABELS[source].split()[0] in described


def test_only_a_real_sifts_source_claims_to_be_one() -> None:
    """"Mapped via SIFTS" asserts authority, so a guess must never say it.

    Checked on the claim form rather than the bare word: the offset fallback's wording
    mentions SIFTS to say it was *unavailable*, which is the opposite of a claim.
    """
    for source in MAPPING_SOURCE_LABELS:
        described = _mapping({1: 2}, source=source).describe()
        claims = "mapped via SIFTS" in described
        assert claims == (source in SIFTS_SOURCES), f"{source}: {described}"


def test_identity_is_not_reported_as_sifts() -> None:
    mapping = identity_mapping("1ABC", "A")
    assert mapping.is_identity
    assert not mapping.is_sifts


def test_describe_reports_how_many_residues_were_aligned() -> None:
    assert "3 residues" in _mapping({1: 2, 2: 3, 3: 4}).describe()


# --------------------------------------------------------------------------------------
# Tier selection
# --------------------------------------------------------------------------------------


class _Response:
    def __init__(self, payload=None, text="", ok=True) -> None:
        self._payload = payload
        self.text = text
        self.ok = ok

    def json(self):
        if self._payload is None:
            raise ValueError("no json")
        return self._payload

    def raise_for_status(self) -> None:
        if not self.ok:
            raise RuntimeError("http error")


def test_tier_1_wins_when_the_api_answers(monkeypatch, tmp_path: Path) -> None:
    payload = _payload("1abc", "P00001", _segment("A", 2, 4, 1, 3))
    monkeypatch.setattr(
        "common.http_lookup.requests.get",
        lambda url, **kwargs: _Response(payload=payload),
    )
    mapping = StructureVizService(tmp_path).build_position_mapping("1ABC", "P00001", "A")
    assert mapping.source == "sifts-api"
    assert mapping.is_sifts
    assert mapping.to_pdb(2) == [1]


def test_tier_2_takes_over_when_the_api_fails(monkeypatch, tmp_path: Path) -> None:
    def fake_get(url, **kwargs):  # noqa: ANN001, ARG001
        if "api/mappings" in url:
            raise RuntimeError("PDBe API down")
        return _Response(text=_CIF_WITH_SIFTS + "x" * 1200)

    monkeypatch.setattr("common.http_lookup.requests.get", fake_get)
    mapping = StructureVizService(tmp_path).build_position_mapping("2IVT", "P07949", "A")
    assert mapping.source == "sifts-mmcif"
    assert mapping.is_sifts


def test_tier_3_takes_over_when_sifts_is_unavailable(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        "common.http_lookup.requests.get",
        lambda url, **kwargs: (_ for _ in ()).throw(RuntimeError("offline")),
    )
    pdb = tmp_path / "x.pdb"
    pdb.write_text(
        "ATOM      1  CA  MET A   1      0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  ALA A   3      0.000   0.000   0.000  1.00  0.00           C\n"
    )
    mapping = StructureVizService(tmp_path).build_position_mapping(
        "1ABC", "P00001", "A", uniprot_range=(705, 707), pdb_path=pdb
    )
    assert mapping.source == "chain-range-offset"
    assert not mapping.is_sifts, "an inferred offset must not claim SIFTS authority"
    assert mapping.to_pdb(705) == [1]


def test_tier_4_is_direct_when_nothing_else_is_available(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        "common.http_lookup.requests.get",
        lambda url, **kwargs: (_ for _ in ()).throw(RuntimeError("offline")),
    )
    mapping = StructureVizService(tmp_path).build_position_mapping("1ABC", "P00001", "A")
    assert mapping.source == "direct"
    assert mapping.is_identity
    assert mapping.to_pdb(918) == [918]


def test_no_pdb_id_means_identity_without_any_request(tmp_path: Path) -> None:
    """AlphaFold models are UniProt-numbered; asking PDBe about them is meaningless."""

    def explode(*args, **kwargs):  # noqa: ANN002, ANN003, ARG001
        raise AssertionError("no HTTP call should be made without a PDB entry")

    service = StructureVizService(tmp_path)
    service.fetch_pdbe_sifts_map = explode  # type: ignore[method-assign]
    mapping = service.build_position_mapping("", "P00001", "A")
    assert mapping.is_identity


def test_the_mapping_is_cached_per_entry_chain(monkeypatch, tmp_path: Path) -> None:
    calls = {"n": 0}
    payload = _payload("1abc", "P00001", _segment("A", 2, 4, 1, 3))

    def counting_get(url, **kwargs):  # noqa: ANN001, ARG001
        calls["n"] += 1
        return _Response(payload=payload)

    monkeypatch.setattr("common.http_lookup.requests.get", counting_get)
    service = StructureVizService(tmp_path)
    first = service.build_position_mapping("1ABC", "P00001", "A")
    second = service.build_position_mapping("1ABC", "P00001", "A")
    assert first is second
    assert calls["n"] == 1
    before = calls["n"]
    other_chain = service.build_position_mapping("1ABC", "P00001", "B")
    assert calls["n"] > before, "a different chain must not reuse another chain's mapping"
    assert other_chain is not first


def test_build_position_mapping_uses_a_bounded_connect_timeout(monkeypatch, tmp_path: Path) -> None:
    seen: dict = {}

    def fake_get(url, headers=None, timeout=None):  # noqa: ANN001, ARG001
        seen["timeout"] = timeout
        return _Response(payload=_payload("1abc", "P00001", _segment("A", 2, 4, 1, 3)))

    monkeypatch.setattr("common.http_lookup.requests.get", fake_get)
    service = StructureVizService(tmp_path)
    service.build_position_mapping("1ABC", "P00001", "A")
    assert seen["timeout"] == (
        StructureVizService.CONNECT_TIMEOUT,
        StructureVizService.LOOKUP_READ_TIMEOUT,
    )
