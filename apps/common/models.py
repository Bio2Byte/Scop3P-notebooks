from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class PeptideSelectionMode(StrEnum):
    UNIQUE_SPANS = "Unique peptide spans"
    ALL_ROWS = "All rows"


@dataclass(slots=True, frozen=True)
class PeptideRow:
    peptide_sequence: str
    peptide_start: int
    peptide_end: int
    uniprot_position: int | None
    modified_residue: str
    score: float | None
