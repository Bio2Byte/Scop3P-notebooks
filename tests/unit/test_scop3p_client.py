"""Scop3P v1 client.

The Scop3P API moved the accession from a query parameter into the path and renamed
every field to snake_case. It also serves its single-page app from a
``GET /scop3p/{catchall}`` route, so a request to a retired endpoint returns 200 OK
with HTML rather than a 404 -- which is what turned that migration into an opaque
``JSONDecodeError`` deep inside requests. Both the field mapping and the content-type
guard are pinned here.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from common.services import Scop3PApiError, Scop3PClient


class _Response:
    def __init__(self, payload: Any = None, status_code: int = 200,
                 content_type: str = "application/json", text: str = "") -> None:
        self._payload = payload
        self.status_code = status_code
        self.headers = {"content-type": content_type}
        self.text = text

    def json(self) -> Any:
        if self._payload is None:
            raise ValueError("Expecting value: line 1 column 1 (char 0)")
        return self._payload


def _patch(monkeypatch, handler) -> list[str]:
    urls: list[str] = []

    def fake_get(url, **kwargs):  # noqa: ANN001
        urls.append(url)
        return handler(url, **kwargs)

    # The request is issued by common.http_lookup now, which is where the shared
    # timeout/retry policy lives, so that is where it has to be patched.
    monkeypatch.setattr("common.http_lookup.requests.get", fake_get)
    return urls


# ---------------------------------------------------------------------------
# URLs
# ---------------------------------------------------------------------------


def test_client_targets_the_v1_path_form(monkeypatch) -> None:
    """v1 takes the accession in the path; the pre-v1 query form is retired."""
    urls = _patch(monkeypatch, lambda url, **kw: _Response([]))
    client = Scop3PClient()
    client.fetch_modifications("P07949")
    client.fetch_peptides_modifications("P07949")

    assert urls == [
        "https://iomics.ugent.be/scop3p/api/v1/proteins/P07949/modifications",
        "https://iomics.ugent.be/scop3p/api/v1/proteins/P07949/peptides",
    ]
    assert all("?accession=" not in url for url in urls)


def test_base_url_is_overridable(monkeypatch) -> None:
    urls = _patch(monkeypatch, lambda url, **kw: _Response([]))
    Scop3PClient(base_url="https://example.test/api/v1/").fetch_modifications("P1")
    assert urls == ["https://example.test/api/v1/proteins/P1/modifications"]


def test_requests_carry_the_configured_timeout(monkeypatch) -> None:
    seen: dict = {}

    def handler(url, **kwargs):  # noqa: ANN001
        seen.update(kwargs)
        return _Response([])

    _patch(monkeypatch, handler)
    Scop3PClient(timeout=7).fetch_modifications("P1")
    # Scop3P is reached through the shared lookup policy, so the bound is the policy's --
    # a short connect and a read timeout well under what a file download would want --
    # rather than whatever timeout this client was constructed with.
    from common import http_lookup as policy

    assert seen["timeout"] == (policy.CONNECT_TIMEOUT, policy.READ_TIMEOUT)
    assert seen["headers"]["Accept"] == "application/json"


# ---------------------------------------------------------------------------
# The catch-all: 200 OK + HTML
# ---------------------------------------------------------------------------


def test_html_from_the_spa_catchall_raises_something_actionable(monkeypatch) -> None:
    """The regression this client exists to prevent.

    A retired endpoint answers 200 with the app shell, so raise_for_status() passes
    and .json() fails with "Expecting value: line 1 column 1". The message must name
    the URL and point at the OpenAPI document instead.
    """
    _patch(monkeypatch, lambda url, **kw: _Response(
        payload=None, content_type="text/html; charset=utf-8", text="<!DOCTYPE html>"))

    with pytest.raises(Scop3PApiError) as excinfo:
        Scop3PClient().fetch_modifications("P07949")

    message = str(excinfo.value)
    assert "text/html" in message
    assert "proteins/P07949/modifications" in message
    assert "openapi.json" in message


def test_http_error_is_reported_with_its_status(monkeypatch) -> None:
    _patch(monkeypatch, lambda url, **kw: _Response(status_code=503))
    with pytest.raises(Scop3PApiError, match="503"):
        Scop3PClient().fetch_modifications("P1")


def test_transport_failure_is_wrapped(monkeypatch) -> None:
    import requests

    def handler(url, **kwargs):  # noqa: ANN001
        raise requests.ConnectionError("name resolution failed")

    _patch(monkeypatch, handler)
    with pytest.raises(Scop3PApiError, match="Could not reach Scop3P"):
        Scop3PClient().fetch_modifications("P1")


def test_malformed_json_is_wrapped(monkeypatch) -> None:
    _patch(monkeypatch, lambda url, **kw: _Response(payload=None))
    with pytest.raises(Scop3PApiError, match="malformed JSON"):
        Scop3PClient().fetch_modifications("P1")


# ---------------------------------------------------------------------------
# modifications
# ---------------------------------------------------------------------------


V1_MODIFICATIONS = [
    {
        "uniprot_position": 687,
        "modification_name": "phosphorylation",
        "modified_residue": "Phosphotyrosine",
        "best_probability": 98.6,
        "evidence_terms": "Experimental",
        "pubmed": "24560924",
        "source": "UniProt, PRIDE",
    },
    {
        "uniprot_position": "bad",
        "modification_name": "phosphorylation",
        "modified_residue": "Phosphoserine",
        "source": "PRIDE",
    },
]


def test_v1_modification_fields_map_to_the_toolkit_column_names(monkeypatch) -> None:
    """Downstream code reads position/residue/name; v1 renamed all three."""
    _patch(monkeypatch, lambda url, **kw: _Response(V1_MODIFICATIONS))
    dataframe = Scop3PClient().fetch_modifications("P07949")

    assert list(dataframe.columns) == [
        "position", "residue", "name", "source", "evidence", "reference", "functionalScore",
    ]
    first = dataframe.iloc[0]
    assert first["position"] == 687
    assert first["residue"] == "Phosphotyrosine"
    assert first["name"] == "phosphorylation"
    assert first["functionalScore"] == 98.6
    assert first["reference"] == "24560924"


def test_unparseable_position_becomes_na_not_an_error(monkeypatch) -> None:
    _patch(monkeypatch, lambda url, **kw: _Response(V1_MODIFICATIONS))
    dataframe = Scop3PClient().fetch_modifications("P07949")
    assert dataframe["position"].tolist() == [687, pd.NA]
    assert str(dataframe["position"].dtype) == "Int64"


def test_an_uncovered_protein_is_an_empty_frame_not_an_error(monkeypatch) -> None:
    """v1 answers 200 with [] for both an uncovered protein and an unknown accession.
    Neither is a failure; Scop3P mainly covers human phosphoproteins."""
    _patch(monkeypatch, lambda url, **kw: _Response([]))
    dataframe = Scop3PClient().fetch_modifications("P0DTD1")
    assert isinstance(dataframe, pd.DataFrame)
    assert dataframe.empty


def test_the_pre_v1_wrapped_shape_is_still_accepted(monkeypatch) -> None:
    """Tolerated so a pinned or proxied older deployment keeps working."""
    _patch(monkeypatch, lambda url, **kw: _Response(
        {"modifications": [{"position": 10, "residue": "S", "name": "Phosphorylation"}]}))
    dataframe = Scop3PClient().fetch_modifications("O00571")
    assert dataframe["position"].tolist() == [10]


# ---------------------------------------------------------------------------
# peptides
# ---------------------------------------------------------------------------


V1_PEPTIDES = [
    {
        "peptide_sequence": "RPSLDSMENQVSVDAFK",
        "peptide_start": 694,
        "peptide_end": 710,
        "peptide_modification_position": 3,
        "uniprot_position": 696,
        "score": 113.118789309606,
        "singly_phosphorylated": 1,
        "modifications_str": "3|[21]Phospho[S]",
    }
]


def test_v1_peptide_fields_map_to_the_toolkit_column_names(monkeypatch) -> None:
    _patch(monkeypatch, lambda url, **kw: _Response(V1_PEPTIDES))
    dataframe = Scop3PClient().fetch_peptides_modifications("P07949")

    for column in ("peptideSequence", "peptideStart", "peptideEnd",
                   "peptideModificationPosition", "uniprotPosition", "score", "label"):
        assert column in dataframe.columns, column
    row = dataframe.iloc[0]
    assert row["peptideStart"] == 694
    assert row["peptideEnd"] == 710
    assert row["uniprotPosition"] == 696
    assert str(dataframe["peptideStart"].dtype) == "Int64"


def test_modified_residue_is_recovered_from_the_sequence(monkeypatch) -> None:
    """v1 dropped modifiedResidue, but the modification position is 1-based within the
    peptide, so the residue is that offset into the sequence. Cross-checked against the
    API's own uniprot_position: 694 + 3 - 1 == 696."""
    _patch(monkeypatch, lambda url, **kw: _Response(V1_PEPTIDES))
    row = Scop3PClient().fetch_peptides_modifications("P07949").iloc[0]
    assert row["modifiedResidue"] == "S"
    assert row["peptideStart"] + row["peptideModificationPosition"] - 1 == row["uniprotPosition"]
    assert row["label"] == "RPSLDSMENQVSVDAFK (694-710) @S696 score=113.118789309606"


def test_modified_residue_survives_a_missing_or_out_of_range_position(monkeypatch) -> None:
    _patch(monkeypatch, lambda url, **kw: _Response([
        {"peptide_sequence": "AAAA", "peptide_start": 1, "peptide_end": 4,
         "peptide_modification_position": None, "uniprot_position": 2},
        {"peptide_sequence": "AAAA", "peptide_start": 1, "peptide_end": 4,
         "peptide_modification_position": 99, "uniprot_position": 2},
    ]))
    dataframe = Scop3PClient().fetch_peptides_modifications("P1")
    assert dataframe["modifiedResidue"].tolist() == ["", ""]


def test_an_explicit_modified_residue_is_not_overwritten(monkeypatch) -> None:
    _patch(monkeypatch, lambda url, **kw: _Response([
        {"peptide_sequence": "AAAA", "peptide_start": 1, "peptide_end": 4,
         "peptide_modification_position": 2, "uniprot_position": 2, "modifiedResidue": "T"},
    ]))
    assert Scop3PClient().fetch_peptides_modifications("P1").iloc[0]["modifiedResidue"] == "T"


def test_no_peptides_is_an_empty_frame(monkeypatch) -> None:
    _patch(monkeypatch, lambda url, **kw: _Response([]))
    dataframe = Scop3PClient().fetch_peptides_modifications("O00571")
    assert isinstance(dataframe, pd.DataFrame)
    assert dataframe.empty


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
    assert Scop3PClient._format_label(row) == "AAAA (?-?) @S? score="
