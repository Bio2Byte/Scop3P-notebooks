from __future__ import annotations

from pathlib import Path

import pandas as pd

from common.structure_viz import B2B_METRIC_COLUMNS, StructureVizService


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
    assert dataframe["Amino acid"].tolist() == ["A", "C"]
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
    assert dataframe["Amino acid"].tolist() == ["A", "C", "D"]
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
