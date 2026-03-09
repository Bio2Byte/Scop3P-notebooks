from __future__ import annotations

from dataclasses import dataclass
import re

import pandas as pd

from .models import PeptideSelectionMode


@dataclass(slots=True, frozen=True)
class ParsedSearch:
    query: str


def positions_to_ranges(positions: list[int]) -> list[tuple[int, int]]:
    if not positions:
        return []

    ordered = sorted(set(int(position) for position in positions))
    ranges: list[tuple[int, int]] = []
    start = ordered[0]
    previous = ordered[0]

    for value in ordered[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append((start, previous))
        start = value
        previous = value

    ranges.append((start, previous))
    return ranges


class PeptideMapperService:
    """Pure mapping/filtering logic for the Peptide Mapper app."""

    range_pattern = re.compile(r"^(\d+)\s*-\s*(\d+)$")
    gte_pattern = re.compile(r"^>=\s*(\d+)$")
    lte_pattern = re.compile(r"^<=\s*(\d+)$")

    @classmethod
    def filter_peptides(cls, dataframe: pd.DataFrame, query: str) -> pd.DataFrame:
        if dataframe is None or dataframe.empty:
            return dataframe

        if not query:
            return dataframe

        clean_query = query.strip()

        range_match = cls.range_pattern.match(clean_query)
        if range_match:
            start, end = int(range_match.group(1)), int(range_match.group(2))
            return dataframe[(dataframe["peptideStart"] <= end) & (dataframe["peptideEnd"] >= start)]

        gte_match = cls.gte_pattern.match(clean_query)
        if gte_match:
            start = int(gte_match.group(1))
            return dataframe[dataframe["peptideEnd"] >= start]

        lte_match = cls.lte_pattern.match(clean_query)
        if lte_match:
            end = int(lte_match.group(1))
            return dataframe[dataframe["peptideStart"] <= end]

        if clean_query.isdigit():
            position = int(clean_query)
            return dataframe[(dataframe["peptideStart"] <= position) & (dataframe["peptideEnd"] >= position)]

        return dataframe[
            dataframe["peptideSequence"].astype(str).str.contains(clean_query, case=False, na=False)
        ]

    @staticmethod
    def build_options(dataframe: pd.DataFrame, mode: PeptideSelectionMode) -> list[tuple[str, str]]:
        if dataframe is None or dataframe.empty:
            return []

        if mode == PeptideSelectionMode.UNIQUE_SPANS:
            grouped = (
                dataframe.groupby(["peptideSequence", "peptideStart", "peptideEnd"], as_index=False)
                .agg(n_mod_sites=("uniprotPosition", "nunique"), max_score=("score", "max"))
            )
            options: list[tuple[str, str]] = []
            for _, row in grouped.iterrows():
                peptide_sequence = str(row["peptideSequence"])
                start = int(row["peptideStart"])
                end = int(row["peptideEnd"])
                max_score = row["max_score"]
                key = f"span|{peptide_sequence}|{start}|{end}"
                label = (
                    f"{peptide_sequence} ({start}-{end}) | "
                    f"modSites={int(row['n_mod_sites'])} maxScore={max_score}"
                )
                options.append((label, key))
            return options

        options = []
        for index, row in dataframe.iterrows():
            options.append((str(row["label"]), f"row|{int(index)}"))
        return options


def map_selection(
    *,
    dataframe_all: pd.DataFrame,
    dataframe_filtered: pd.DataFrame,
    selected_keys: list[str],
    mode: PeptideSelectionMode,
    mods_scope: str,
) -> tuple[list[tuple[int, int]], list[int], list[int]]:
    spans: list[tuple[int, int]] = []
    mod_positions: list[int] = []

    if mode == PeptideSelectionMode.UNIQUE_SPANS:
        parsed = []
        for key in selected_keys:
            _, peptide_sequence, start, end = key.split("|", maxsplit=3)
            parsed.append((peptide_sequence, int(start), int(end)))
            spans.append((int(start), int(end)))

        if mods_scope == "Selected peptides only":
            for peptide_sequence, start, end in parsed:
                subset = dataframe_all[
                    (dataframe_all["peptideSequence"] == peptide_sequence)
                    & (dataframe_all["peptideStart"] == start)
                    & (dataframe_all["peptideEnd"] == end)
                ]
                mod_positions.extend(subset["uniprotPosition"].dropna().astype(int).tolist())
        else:
            mod_positions = dataframe_all["uniprotPosition"].dropna().astype(int).tolist()
    else:
        indices = [int(value.split("|", maxsplit=1)[1]) for value in selected_keys]
        subset = dataframe_filtered.loc[indices].copy()
        spans = [
            (int(row["peptideStart"]), int(row["peptideEnd"]))
            for _, row in subset.iterrows()
        ]
        if mods_scope == "Selected peptides only":
            mod_positions = subset["uniprotPosition"].dropna().astype(int).tolist()
        else:
            mod_positions = dataframe_all["uniprotPosition"].dropna().astype(int).tolist()

    position_lists = [list(range(start, end + 1)) for start, end in spans]
    union_positions = sorted(set(position for positions in position_lists for position in positions))
    intersection_positions = (
        sorted(set(position_lists[0]).intersection(*map(set, position_lists[1:])))
        if len(position_lists) > 1
        else []
    )

    return positions_to_ranges(union_positions), intersection_positions, sorted(set(mod_positions))
