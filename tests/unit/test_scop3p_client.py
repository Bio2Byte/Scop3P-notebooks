from __future__ import annotations

import pandas as pd

from apps.common.services import Scop3PClient


class _DummyResponse:
    def __init__(self, payload: dict) -> None:
        self._payload = payload

    def raise_for_status(self) -> None:
        return None

    def json(self) -> dict:
        return self._payload


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

    def _mock_get(*args, **kwargs):
        return _DummyResponse(payload)

    monkeypatch.setattr("apps.common.services.requests.get", _mock_get)

    client = Scop3PClient()
    dataframe = client.fetch_peptides_modifications("O00571")

    assert isinstance(dataframe, pd.DataFrame)
    assert dataframe.iloc[0]["peptideStart"] == 10
    assert dataframe.iloc[0]["uniprotPosition"] == 12
    assert "label" in dataframe.columns
