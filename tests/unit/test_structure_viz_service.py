from __future__ import annotations

from pathlib import Path

import pandas as pd

from common.structure_viz import StructureVizService


class _TextResponse:
    def __init__(self, text: str) -> None:
        self.text = text

    def raise_for_status(self) -> None:
        return None


def test_fetch_ptms_normalizes_position_column(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        "apps.common.services.Scop3pRestApi.fetch_modifications",
        lambda self, accession: {
            "modifications": [
                {"residue": "S", "position": "10", "name": "Phosphorylation"},
                {"residue": "T", "position": "bad", "name": "Phosphorylation"},
            ]
        },
    )

    dataframe = StructureVizService(tmp_path).fetch_ptms("P12345")
    assert dataframe["position"].tolist() == [10, pd.NA]
    assert str(dataframe["position"].dtype) == "Int64"


def test_fetch_ptms_returns_empty_dataframe(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        "apps.common.services.Scop3pRestApi.fetch_modifications",
        lambda self, accession: {"modifications": []},
    )

    dataframe = StructureVizService(tmp_path).fetch_ptms("P12345")
    assert dataframe.empty


def test_fetch_sequence_parses_fasta(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        "apps.common.structure_viz.requests.get",
        lambda *args, **kwargs: _TextResponse(">sp|P12345|\nACD\nEFG\n"),
    )

    sequence = StructureVizService(tmp_path).fetch_sequence("P12345")
    assert sequence == "ACDEFG"
