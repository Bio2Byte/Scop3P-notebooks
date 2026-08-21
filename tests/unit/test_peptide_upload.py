from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from common.models import PeptideSelectionMode
from common.peptide_mapper import (
    PeptideColumnMapping,
    PeptideMapperService,
    build_upload_mapping,
    detect_peptide_columns,
    map_selection,
    mapped_residue_rows,
    normalize_protein_id,
    peptides_for_protein,
    protein_choices,
    read_peptide_table,
)

NOTEBOOK_HEADERS = ["ACC_ID", "Pep_seq", "pep_start", "pep_end", "UP_POS"]
MAXQUANT_HEADERS = ["Proteins", "Sequence", "Start position", "End position", "Position"]


# ---------------------------------------------------------------------------
# Column detection
# ---------------------------------------------------------------------------


def test_detects_the_notebooks_own_headers() -> None:
    mapping = detect_peptide_columns(NOTEBOOK_HEADERS)
    assert mapping.is_complete()
    assert (mapping.protein, mapping.sequence) == ("ACC_ID", "Pep_seq")
    assert (mapping.start, mapping.end, mapping.position) == ("pep_start", "pep_end", "UP_POS")


def test_detects_maxquant_style_headers() -> None:
    """"Start position" must match the "start" candidate despite the space."""
    mapping = detect_peptide_columns(MAXQUANT_HEADERS)
    assert mapping.is_complete()
    assert mapping.protein == "Proteins"
    assert mapping.sequence == "Sequence"
    assert mapping.start == "Start position"
    assert mapping.end == "End position"
    assert mapping.position == "Position"


def test_detects_repo_canonical_headers() -> None:
    mapping = detect_peptide_columns(
        ["accession", "peptideSequence", "peptideStart", "peptideEnd", "uniprotPosition"]
    )
    assert mapping.is_complete()


def test_detection_is_case_insensitive() -> None:
    mapping = detect_peptide_columns(["acc_id", "PEP_SEQ", "Pep_Start", "PEP_end", "up_pos"])
    assert mapping.is_complete()


def test_a_missing_field_is_reported_not_guessed() -> None:
    """The notebook fell back to columns[0], which silently built a wrong table."""
    mapping = detect_peptide_columns(["ACC_ID", "Pep_seq", "pep_start", "pep_end"])
    assert not mapping.is_complete()
    assert mapping.missing == ("position",)
    assert mapping.position is None
    assert mapping.protein == "ACC_ID"


def test_completely_unrecognised_headers_report_every_field() -> None:
    mapping = detect_peptide_columns(["a", "b", "c"])
    assert set(mapping.missing) == {"protein", "sequence", "start", "end", "position"}


# ---------------------------------------------------------------------------
# Protein identifier normalisation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value,expected", [
    ("P07949", "P07949"),
    ("  P07949 ", "P07949"),
    ("sp|P07949|RET_HUMAN", "P07949"),
    ("tr|A0A123|A0A123_HUMAN", "A0A123"),
    ("P07949;Q12345", "P07949"),
    ("P07949,Q12345", "P07949"),
    ("sp|P07949|RET_HUMAN;sp|Q12345|X", "P07949"),
    ("P07949-2", "P07949-2"),   # isoform suffix is meaningful; keep it
    ("", ""),
])
def test_normalize_protein_id(value: str, expected: str) -> None:
    assert normalize_protein_id(value) == expected


# ---------------------------------------------------------------------------
# read_peptide_table
# ---------------------------------------------------------------------------


def test_reads_tab_separated(tmp_path: Path) -> None:
    path = tmp_path / "p.tsv"
    path.write_text("ACC_ID\tPep_seq\nP07949\tSSFGDVLLSK\n")
    frame, delimiter = read_peptide_table(path)
    assert delimiter == "tab"
    assert list(frame.columns) == ["ACC_ID", "Pep_seq"]


def test_falls_back_to_comma_separated(tmp_path: Path) -> None:
    path = tmp_path / "p.csv"
    path.write_text("ACC_ID,Pep_seq\nP07949,SSFGDVLLSK\n")
    frame, delimiter = read_peptide_table(path)
    assert delimiter == "comma"
    assert frame.shape == (1, 2)


def test_a_single_column_result_is_rejected(tmp_path: Path) -> None:
    """Otherwise column detection fails later with a far less useful message."""
    path = tmp_path / "p.txt"
    path.write_text("ACC_ID|Pep_seq\nP07949|SSFGDVLLSK\n")
    with pytest.raises(ValueError, match="delimiter"):
        read_peptide_table(path)


# ---------------------------------------------------------------------------
# build_upload_mapping
# ---------------------------------------------------------------------------


def _raw(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=NOTEBOOK_HEADERS)


def _mapping() -> PeptideColumnMapping:
    return PeptideColumnMapping(
        protein="ACC_ID", sequence="Pep_seq", start="pep_start", end="pep_end", position="UP_POS"
    )


def test_mapping_emits_the_canonical_schema() -> None:
    frame = build_upload_mapping(
        _raw([{"ACC_ID": "P07949", "Pep_seq": "SSFGDVLLSK",
               "pep_start": 100, "pep_end": 109, "UP_POS": 104}]),
        _mapping(),
    )
    assert set(frame.columns) == {
        "accession", "peptideSequence", "peptideStart", "peptideEnd",
        "uniprotPosition", "label",
    }
    row = frame.iloc[0]
    assert row["accession"] == "P07949"
    assert (row["peptideStart"], row["peptideEnd"], row["uniprotPosition"]) == (100, 109, 104)
    assert row["label"] == "SSFGDVLLSK (100-109) @UP_POS=104"


def test_mapping_normalises_fasta_style_protein_ids() -> None:
    frame = build_upload_mapping(
        _raw([{"ACC_ID": "sp|P07949|RET_HUMAN", "Pep_seq": "SSFGDVLLSK",
               "pep_start": 100, "pep_end": 109, "UP_POS": 104}]),
        _mapping(),
    )
    assert frame.iloc[0]["accession"] == "P07949"


def test_mapping_drops_unusable_rows() -> None:
    frame = build_upload_mapping(
        _raw([
            {"ACC_ID": "P07949", "Pep_seq": "AAA", "pep_start": 10, "pep_end": 20, "UP_POS": 15},
            {"ACC_ID": "P07949", "Pep_seq": "BBB", "pep_start": "x", "pep_end": 20, "UP_POS": 15},
            {"ACC_ID": "P07949", "Pep_seq": "CCC", "pep_start": 30, "pep_end": 20, "UP_POS": 35},
            {"ACC_ID": "P07949", "Pep_seq": "DDD", "pep_start": 0, "pep_end": 5, "UP_POS": 3},
            {"ACC_ID": "", "Pep_seq": "EEE", "pep_start": 1, "pep_end": 9, "UP_POS": 4},
            {"ACC_ID": "P07949", "Pep_seq": "FFF", "pep_start": 40, "pep_end": 50, "UP_POS": None},
        ]),
        _mapping(),
    )
    assert list(frame["peptideSequence"]) == ["AAA"]


def test_mapping_index_is_contiguous_after_dropping_rows() -> None:
    """map_selection addresses ALL_ROWS selections with .loc against iterrows keys."""
    frame = build_upload_mapping(
        _raw([
            {"ACC_ID": "P1", "Pep_seq": "BAD", "pep_start": "x", "pep_end": 2, "UP_POS": 1},
            {"ACC_ID": "P1", "Pep_seq": "AAA", "pep_start": 10, "pep_end": 20, "UP_POS": 15},
            {"ACC_ID": "P1", "Pep_seq": "BBB", "pep_start": 30, "pep_end": 40, "UP_POS": 35},
        ]),
        _mapping(),
    )
    assert list(frame.index) == [0, 1]


def test_mapping_requires_every_field() -> None:
    with pytest.raises(ValueError, match="position"):
        build_upload_mapping(_raw([]), PeptideColumnMapping(
            protein="ACC_ID", sequence="Pep_seq", start="pep_start", end="pep_end"
        ))


def test_mapping_of_an_empty_table_is_empty() -> None:
    assert build_upload_mapping(_raw([]), _mapping()).empty


# ---------------------------------------------------------------------------
# The reuse claim, executable
# ---------------------------------------------------------------------------


UPLOAD_ROWS = [
    {"ACC_ID": "P07949", "Pep_seq": "SSFGDVLLSK", "pep_start": 100, "pep_end": 109, "UP_POS": 104},
    {"ACC_ID": "P07949", "Pep_seq": "SSFGDVLLSK", "pep_start": 100, "pep_end": 109, "UP_POS": 106},
    {"ACC_ID": "P07949", "Pep_seq": "AGQKPIYIVM", "pep_start": 105, "pep_end": 114, "UP_POS": 110},
    {"ACC_ID": "Q12345", "Pep_seq": "MMMKKK", "pep_start": 1, "pep_end": 6, "UP_POS": 3},
]


@pytest.mark.parametrize("mode", list(PeptideSelectionMode))
def test_uploaded_tables_flow_through_the_scop3p_pipeline_unchanged(mode) -> None:
    """The point of normalising the upload: nothing downstream needs to know.

    build_options and map_selection are the Scop3P path's own helpers, used here on
    an uploaded frame with no adaptation.
    """
    frame = build_upload_mapping(_raw(UPLOAD_ROWS), _mapping())
    subset = peptides_for_protein(frame, "P07949")
    assert len(subset) == 3

    options = PeptideMapperService.build_options(subset, mode)
    assert options
    keys = [value for _label, value in options]

    union_ranges, intersection, modifications = map_selection(
        dataframe_all=subset,
        dataframe_filtered=subset,
        selected_keys=keys,
        mode=mode,
        mods_scope="All peptides",
    )
    assert union_ranges == [(100, 114)]
    # The two spans overlap at 105-109, which is what the red intersection shows.
    assert intersection == list(range(105, 110))
    assert modifications == [104, 106, 110]


def test_upload_labels_omit_maxscore_when_there_is_no_score_column() -> None:
    frame = build_upload_mapping(_raw(UPLOAD_ROWS), _mapping())
    options = PeptideMapperService.build_options(frame, PeptideSelectionMode.UNIQUE_SPANS)
    labels = [label for label, _ in options]
    assert all("maxScore" not in label for label in labels)
    assert any("modSites=" in label for label in labels)


def test_scop3p_labels_keep_maxscore() -> None:
    """The score-tolerance change must not alter the Scop3P path's labels."""
    scop3p = pd.DataFrame([
        {"peptideSequence": "AAA", "peptideStart": 1, "peptideEnd": 3,
         "uniprotPosition": 2, "score": 0.9, "label": "AAA"},
    ])
    label = PeptideMapperService.build_options(scop3p, PeptideSelectionMode.UNIQUE_SPANS)[0][0]
    assert label == "AAA (1-3) | modSites=1 maxScore=0.9"


def test_uploaded_frames_are_searchable_with_the_existing_filter_syntax() -> None:
    frame = build_upload_mapping(_raw(UPLOAD_ROWS), _mapping())
    assert len(PeptideMapperService.filter_peptides(frame, "SSFGD")) == 2
    assert len(PeptideMapperService.filter_peptides(frame, ">=110")) == 1
    assert len(PeptideMapperService.filter_peptides(frame, "100-109")) == 3
    assert len(PeptideMapperService.filter_peptides(frame, "")) == 4


# ---------------------------------------------------------------------------
# Protein picker and export helpers
# ---------------------------------------------------------------------------


def test_protein_choices_counts_peptides_and_sites() -> None:
    frame = build_upload_mapping(_raw(UPLOAD_ROWS), _mapping())
    choices = protein_choices(frame)
    assert list(choices) == ["P07949", "Q12345"]          # sorted
    assert choices["P07949"] == "P07949 (2 peptides, 3 sites)"
    assert choices["Q12345"] == "Q12345 (1 peptides, 1 sites)"


def test_protein_choices_of_an_empty_frame() -> None:
    assert protein_choices(pd.DataFrame()) == {}


def test_peptides_for_protein_reindexes() -> None:
    frame = build_upload_mapping(_raw(UPLOAD_ROWS), _mapping())
    subset = peptides_for_protein(frame, "Q12345")
    assert len(subset) == 1
    assert list(subset.index) == [0]


def test_peptides_for_an_unknown_protein_is_empty() -> None:
    frame = build_upload_mapping(_raw(UPLOAD_ROWS), _mapping())
    assert peptides_for_protein(frame, "NOPE").empty
    assert peptides_for_protein(frame, "").empty


def test_mapped_residue_rows_shape() -> None:
    table = mapped_residue_rows([(10, 12), (20, 21)], [11], [10, 21])
    assert list(table.columns) == ["type", "start", "end", "position"]
    assert list(table["type"]) == [
        "peptide_span", "peptide_span", "intersection", "modification", "modification"
    ]
    assert table.iloc[0]["start"] == 10 and table.iloc[0]["end"] == 12
    assert table.iloc[2]["position"] == 11


def test_mapped_residue_rows_of_nothing_is_an_empty_table() -> None:
    table = mapped_residue_rows([], [], [])
    assert table.empty
    assert list(table.columns) == ["type", "start", "end", "position"]


# ---------------------------------------------------------------------------
# The app's source-tab wiring
# ---------------------------------------------------------------------------


def test_upload_tab_constant_matches_the_rendered_panel_title() -> None:
    """`active_accession()` compares input.source_tabs() to this string.

    A navset with an id reports its ACTIVE PANEL TITLE, so if the panel is renamed
    and the constant is not, every upload-mode render silently falls back to reading
    the Scop3P accession field instead. Nothing errors; the wrong protein is fetched.
    """
    from peptide_mapper.app import UPLOAD_TAB, app_ui

    rendered = str(app_ui.tagify())
    assert f">{UPLOAD_TAB}<" in rendered, (
        f"no nav panel titled {UPLOAD_TAB!r}; the constant and the UI have drifted"
    )


def test_the_scop3p_source_tab_is_still_first() -> None:
    """It must stay the initially-active panel: the Scop3P flow is the default, and
    the upload controls are the opt-in."""
    from peptide_mapper.app import UPLOAD_TAB, app_ui

    rendered = str(app_ui.tagify())
    assert rendered.index(">Scop3P peptides<") < rendered.index(f">{UPLOAD_TAB}<")
