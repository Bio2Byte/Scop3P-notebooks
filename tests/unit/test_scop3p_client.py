from __future__ import annotations

import pandas as pd

from common.services import Scop3PClient


def test_fetch_peptides_modifications_normalizes_payload(monkeypatch) -> None:
    payload = {
        "peptides": [
            {
                "peptideSequence": "AAAA",
                "peptideStart": "10",
                "peptideEnd": "20",
                "peptideModificationPosition": "3",
                "uniprotPosition": "12",
                "modifiedResidue": "S",
                "score": "0.95",
            }
        ]
    }

    monkeypatch.setattr(
        "apps.common.services.Scop3pRestApi.fetch_peptides",
        lambda self, accession: payload,
    )

    client = Scop3PClient()
    dataframe = client.fetch_peptides_modifications("O00571")

    assert isinstance(dataframe, pd.DataFrame)
    assert dataframe.iloc[0]["peptideStart"] == 10
    assert dataframe.iloc[0]["uniprotPosition"] == 12
    assert "label" in dataframe.columns


def test_fetch_peptides_modifications_returns_empty_dataframe(monkeypatch) -> None:
    monkeypatch.setattr(
        "apps.common.services.Scop3pRestApi.fetch_peptides",
        lambda self, accession: {"peptides": []},
    )

    dataframe = Scop3PClient().fetch_peptides_modifications("O00571")
    assert isinstance(dataframe, pd.DataFrame)
    assert dataframe.empty


def test_fetch_modifications_normalizes_position_column(monkeypatch) -> None:
    monkeypatch.setattr(
        "apps.common.services.Scop3pRestApi.fetch_modifications",
        lambda self, accession: {
            "modifications": [
                {"position": "10", "residue": "S", "name": "Phosphorylation"},
                {"position": "bad", "residue": "T", "name": "Phosphorylation"},
            ]
        },
    )

    dataframe = Scop3PClient().fetch_modifications("O00571")
    assert dataframe["position"].tolist() == [10, pd.NA]
    assert str(dataframe["position"].dtype) == "Int64"


def test_format_label_uses_placeholders_for_missing_numeric_values() -> None:
    row = pd.Series(
        {
            "peptideSequence": "AAAA",
            "peptideStart": pd.NA,
            "peptideEnd": pd.NA,
            "modifiedResidue": "S",
            "uniprotPosition": pd.NA,
            "score": "",
        }
    )

    label = Scop3PClient._format_label(row)
    assert label == "AAAA (?-?) @S? score="
