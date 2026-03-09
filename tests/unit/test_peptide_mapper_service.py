from __future__ import annotations

import pandas as pd

from apps.common.models import PeptideSelectionMode
from apps.common.peptide_mapper import PeptideMapperService, map_selection, positions_to_ranges


def _sample_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "peptideSequence": "AAAA",
                "peptideStart": 10,
                "peptideEnd": 20,
                "uniprotPosition": 12,
                "score": 0.9,
                "label": "AAAA (10-20) @S12 score=0.9",
            },
            {
                "peptideSequence": "BBBB",
                "peptideStart": 18,
                "peptideEnd": 30,
                "uniprotPosition": 22,
                "score": 0.8,
                "label": "BBBB (18-30) @T22 score=0.8",
            },
            {
                "peptideSequence": "AAAA",
                "peptideStart": 10,
                "peptideEnd": 20,
                "uniprotPosition": 14,
                "score": 0.7,
                "label": "AAAA (10-20) @Y14 score=0.7",
            },
        ]
    )


def test_positions_to_ranges_merges_consecutive_positions() -> None:
    assert positions_to_ranges([1, 2, 3, 7, 8, 10]) == [(1, 3), (7, 8), (10, 10)]


def test_filter_peptides_range_query() -> None:
    dataframe = _sample_df()
    filtered = PeptideMapperService.filter_peptides(dataframe, "15-23")
    assert len(filtered) == 3


def test_filter_peptides_sequence_query() -> None:
    dataframe = _sample_df()
    filtered = PeptideMapperService.filter_peptides(dataframe, "bbbb")
    assert len(filtered) == 1
    assert filtered.iloc[0]["peptideSequence"] == "BBBB"


def test_build_options_unique_spans() -> None:
    dataframe = _sample_df()
    options = PeptideMapperService.build_options(dataframe, PeptideSelectionMode.UNIQUE_SPANS)
    assert len(options) == 2
    assert options[0][1].startswith("span|")


def test_map_selection_unique_spans_selected_mods() -> None:
    dataframe = _sample_df()
    filtered = dataframe.copy()

    union_ranges, intersection, mods = map_selection(
        dataframe_all=dataframe,
        dataframe_filtered=filtered,
        selected_keys=["span|AAAA|10|20", "span|BBBB|18|30"],
        mode=PeptideSelectionMode.UNIQUE_SPANS,
        mods_scope="Selected peptides only",
    )

    assert union_ranges == [(10, 30)]
    assert intersection == [18, 19, 20]
    assert mods == [12, 14, 22]


def test_map_selection_all_rows_all_mods() -> None:
    dataframe = _sample_df()
    filtered = dataframe.iloc[[0, 1]].copy()

    union_ranges, intersection, mods = map_selection(
        dataframe_all=dataframe,
        dataframe_filtered=filtered,
        selected_keys=["row|0", "row|1"],
        mode=PeptideSelectionMode.ALL_ROWS,
        mods_scope="All protein mods",
    )

    assert union_ranges == [(10, 30)]
    assert intersection == [18, 19, 20]
    assert mods == [12, 14, 22]
