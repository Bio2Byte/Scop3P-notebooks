from __future__ import annotations

from dataclasses import dataclass
import re

import pandas as pd

from .models import PeptideSelectionMode


@dataclass(slots=True, frozen=True)
class ParsedSearch:
    query: str


def _format_score(value: object) -> str:
    """A peptide score short enough to read in a dropdown chip.

    Four significant figures rather than fixed decimals. The upstream aggregate carries
    full float precision (``278.273526557391``), which is noise in a label and helped push
    the chip's remove control out of view. Significant figures are the right choice here
    because the column holds both intensities in the hundreds and probabilities near zero:
    a fixed ``.2f`` would render a score of 1e-5 as "0.00" and lose it entirely.

    A value that already reads short is unchanged: ``0.9`` stays ``0.9``.
    """
    try:
        return f"{float(value):.4g}"
    except (TypeError, ValueError):
        return str(value)


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
            # `score` only exists on the Scop3P path. An uploaded peptide table has
            # no score column, so aggregate it only when present and drop the
            # maxScore fragment from the label rather than showing "maxScore=<NA>".
            # Scop3P labels are unchanged.
            has_score = "score" in dataframe.columns
            aggregations = {"n_mod_sites": ("uniprotPosition", "nunique")}
            if has_score:
                aggregations["max_score"] = ("score", "max")
            grouped = (
                dataframe.groupby(["peptideSequence", "peptideStart", "peptideEnd"], as_index=False)
                .agg(**aggregations)
            )
            options: list[tuple[str, str]] = []
            for _, row in grouped.iterrows():
                peptide_sequence = str(row["peptideSequence"])
                start = int(row["peptideStart"])
                end = int(row["peptideEnd"])
                key = f"span|{peptide_sequence}|{start}|{end}"
                label = f"{peptide_sequence} ({start}-{end}) | modSites={int(row['n_mod_sites'])}"
                if has_score:
                    # Rounded: the raw aggregate arrives as 278.273526557391, which
                    # added 12 meaningless digits to every chip in the selector.
                    label = f"{label} maxScore={_format_score(row['max_score'])}"
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


# ---------------------------------------------------------------------------
# Uploaded peptide tables
#
# Ported from notebooks/Peptide_mapper_fileupload_voila.ipynb. The design decision
# that makes this cheap: an uploaded table is normalised into exactly the column
# schema the Scop3P path already produces, so filter_peptides, build_options and
# map_selection above are reused unchanged and everything downstream of the peptide
# list -- search, selection, AlphaFold download, NGL rendering, exports -- is
# source-agnostic.
# ---------------------------------------------------------------------------

#: Header names to look for, per canonical field. The notebook's five ``guess_column``
#: candidate lists, extended with the headers MaxQuant, FragPipe and DIA-NN actually
#: emit so a real search-engine export is recognised without hand-mapping.
PEPTIDE_COLUMN_CANDIDATES: dict[str, tuple[str, ...]] = {
    "protein": (
        "ACC_ID", "Protein", "protein", "Uniprot", "UniProt", "accession",
        "Protein ID", "Proteins", "Leading razor protein", "Protein.Group",
        "Protein.Ids", "Master Protein Accessions",
    ),
    "sequence": (
        "Pep_seq", "peptideSequence", "peptide", "Peptide", "Sequence",
        "Modified sequence", "Stripped.Sequence", "Annotated Sequence",
    ),
    "start": (
        "pep_start", "peptideStart", "PeptideStart", "start", "Start position",
        "Start", "pep_start_prot",
    ),
    "end": (
        "pep_end", "peptideEnd", "PeptideEnd", "end", "End position", "End",
        "pep_end_prot",
    ),
    "position": (
        "UP_POS", "uniprotPosition", "UniprotPosition", "modpos_prot", "Position",
        "Positions within proteins",
    ),
}

_REQUIRED_FIELDS = ("protein", "sequence", "start", "end", "position")


@dataclass(slots=True, frozen=True)
class PeptideColumnMapping:
    """Which uploaded column supplies each canonical field."""

    protein: str | None = None
    sequence: str | None = None
    start: str | None = None
    end: str | None = None
    position: str | None = None

    @property
    def missing(self) -> tuple[str, ...]:
        return tuple(field for field in _REQUIRED_FIELDS if getattr(self, field) is None)

    def is_complete(self) -> bool:
        return not self.missing


def _normalise_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def detect_peptide_columns(columns: list[str] | tuple[str, ...]) -> PeptideColumnMapping:
    """Guess which column supplies each field, in three passes of decreasing strictness.

    Unlike the notebook, an unmatched field stays ``None`` instead of falling back to
    ``columns[0]``. That fallback produced a mapping that looked valid, built a
    nonsense peptide table, and then failed somewhere much later; a missing field
    named up front is far easier to act on.
    """
    available = list(columns)
    exact = {str(column).lower(): column for column in available}
    normalised = {_normalise_header(column): column for column in available}

    resolved: dict[str, str | None] = {}
    for field, candidates in PEPTIDE_COLUMN_CANDIDATES.items():
        match: str | None = None

        # 1. Case-insensitive exact match, which is what the notebook did.
        for candidate in candidates:
            if candidate.lower() in exact:
                match = exact[candidate.lower()]
                break

        # 2. Ignore punctuation and spacing, so "Start position" == "start_position".
        if match is None:
            for candidate in candidates:
                key = _normalise_header(candidate)
                if key in normalised:
                    match = normalised[key]
                    break

        # 3. Containment, for headers carrying extra qualifiers -- "Start position"
        #    or "Positions within proteins". Only the candidate-inside-column
        #    direction: the reverse lets a one-character column name like "a" match
        #    "accession" and mis-map the whole table. Both sides are length-guarded
        #    for the same reason.
        if match is None:
            for candidate in candidates:
                key = _normalise_header(candidate)
                if len(key) < 4:
                    continue
                for column_key, column in normalised.items():
                    if len(column_key) >= 4 and key in column_key:
                        match = column
                        break
                if match is not None:
                    break

        resolved[field] = match

    return PeptideColumnMapping(**resolved)


def normalize_protein_id(value: object) -> str:
    """Reduce a protein identifier to a bare UniProt accession.

    Search engines emit ``sp|P07949|RET_HUMAN``, semicolon-separated protein groups,
    and padded values; the AlphaFold URL needs the accession alone. The notebook did
    none of this, so any FASTA-style identifier produced a 404 on download.
    """
    text = str(value).strip()
    if not text:
        return ""
    text = text.split(";")[0].strip()   # protein group -> leading entry
    text = text.split(",")[0].strip()
    if "|" in text:
        parts = [part for part in text.split("|") if part]
        if len(parts) >= 2 and parts[0].lower() in {"sp", "tr"}:
            return parts[1].strip()
        return parts[1].strip() if len(parts) >= 2 else parts[0].strip()
    return text


def read_peptide_table(path, sep: str | None = None) -> tuple[pd.DataFrame, str]:
    """Read an uploaded peptide table, tab first then comma.

    Returns the frame and the delimiter used. A single-column result means the
    delimiter guess was wrong, which is worth saying explicitly rather than letting
    column detection fail mysteriously afterwards.
    """
    if sep is not None:
        return pd.read_csv(path, sep=sep), sep

    for candidate, name in (("\t", "tab"), (",", "comma")):
        try:
            frame = pd.read_csv(path, sep=candidate)
        except Exception:
            continue
        if frame.shape[1] > 1:
            return frame, name
    raise ValueError(
        "Could not read the table as tab- or comma-separated with more than one "
        "column. Check the file's delimiter."
    )


def build_upload_mapping(
    dataframe: pd.DataFrame, mapping: PeptideColumnMapping
) -> pd.DataFrame:
    """Normalise an uploaded table into the canonical peptide schema.

    Emits ``accession``, ``peptideSequence``, ``peptideStart``, ``peptideEnd`` and
    ``uniprotPosition``, plus the ``label`` that ``build_options`` reads in ALL_ROWS
    mode. The index is reset because ``map_selection`` addresses rows by ``.loc``
    against keys minted from ``iterrows()``, so it has to be unique and contiguous.
    """
    if not mapping.is_complete():
        raise ValueError(
            "Select a column for each of: " + ", ".join(mapping.missing) + "."
        )
    if dataframe is None or dataframe.empty:
        return pd.DataFrame()

    columns = {
        "accession": mapping.protein,
        "peptideSequence": mapping.sequence,
        "peptideStart": mapping.start,
        "peptideEnd": mapping.end,
        "uniprotPosition": mapping.position,
    }
    frame = dataframe[list(columns.values())].copy()
    frame.columns = list(columns)

    frame["accession"] = frame["accession"].map(normalize_protein_id)
    frame["peptideSequence"] = frame["peptideSequence"].astype(str).str.strip()
    for numeric in ("peptideStart", "peptideEnd", "uniprotPosition"):
        frame[numeric] = pd.to_numeric(frame[numeric], errors="coerce")

    frame = frame.replace({"accession": "", "peptideSequence": ""}, pd.NA)
    frame = frame.dropna(
        subset=["accession", "peptideSequence", "peptideStart", "peptideEnd", "uniprotPosition"]
    )
    for numeric in ("peptideStart", "peptideEnd", "uniprotPosition"):
        frame[numeric] = frame[numeric].astype(int)

    # Notebook parity: 1-indexed starts, and spans that do not run backwards.
    frame = frame[(frame["peptideStart"] >= 1) & (frame["peptideEnd"] >= frame["peptideStart"])]
    frame = frame.reset_index(drop=True)

    frame["label"] = [
        f"{row.peptideSequence} ({row.peptideStart}-{row.peptideEnd}) "
        f"@UP_POS={row.uniprotPosition}"
        for row in frame.itertuples()
    ]
    return frame


def protein_choices(dataframe: pd.DataFrame) -> dict[str, str]:
    """``{accession: "accession (n peptides, m sites)"}``, ready for update_select."""
    if dataframe is None or dataframe.empty:
        return {}
    grouped = dataframe.groupby("accession").agg(
        peptides=("peptideSequence", "nunique"),
        sites=("uniprotPosition", "nunique"),
    )
    return {
        str(accession): f"{accession} ({row.peptides} peptides, {row.sites} sites)"
        for accession, row in grouped.sort_index().iterrows()
    }


def peptides_for_protein(dataframe: pd.DataFrame, accession: str) -> pd.DataFrame:
    """Rows for one protein, re-indexed so map_selection's .loc keys stay valid."""
    if dataframe is None or dataframe.empty or not accession:
        return pd.DataFrame(columns=getattr(dataframe, "columns", None))
    return dataframe[dataframe["accession"] == accession].reset_index(drop=True)


def mapped_residue_rows(
    union_ranges: list[tuple[int, int]],
    intersection_positions: list[int],
    modification_positions: list[int],
) -> pd.DataFrame:
    """The exported residue table: one row per mapped span, residue or site."""
    rows: list[dict[str, object]] = []
    for start, end in union_ranges or []:
        rows.append({"type": "peptide_span", "start": start, "end": end, "position": ""})
    for position in intersection_positions or []:
        rows.append({"type": "intersection", "start": "", "end": "", "position": position})
    for position in modification_positions or []:
        rows.append({"type": "modification", "start": "", "end": "", "position": position})
    return pd.DataFrame(rows, columns=["type", "start", "end", "position"])
