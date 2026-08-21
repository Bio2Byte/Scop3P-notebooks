from __future__ import annotations

from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
import urllib.request

import pandas as pd
import requests

from common.cache import memoize
from common.http_lookup import lookup as http_lookup
from common.logging_utils import get_logger


LOGGER = get_logger("scop3p.common.services")


class Scop3PApiError(RuntimeError):
    """Scop3P returned something that is not the JSON we asked for."""


class Scop3PClient:
    """Client for the Scop3P v1 REST API.

    Talks to the API directly rather than through the ``scop3p`` PyPI package. That
    package (1.1.0, the latest release) still targets the pre-v1 query-string
    endpoints -- ``/scop3p/api/modifications?accession=...`` and friends -- which the
    current deployment no longer serves.

    Two things about that failure are worth recording, because they are what made it
    hard to diagnose:

    * The API has a ``GET /scop3p/{catchall}`` route that serves the single-page-app
      HTML. A request to a retired endpoint therefore comes back **200 OK** with
      ``content-type: text/html``, so ``raise_for_status()`` is happy and the error
      only surfaces as ``JSONDecodeError: Expecting value: line 1 column 1`` from deep
      inside ``requests``. :meth:`_get_json` checks the content type for exactly this
      reason and raises something actionable instead.
    * v1 moved the accession from a query parameter into the path and renamed every
      field from camelCase to snake_case. The frames returned here keep the original
      column names, so every consumer in the toolkit is unaffected by that rename.
    """

    #: v1 base. Endpoints are ``/proteins/{accession}/<resource>``.
    DEFAULT_BASE_URL = "https://iomics.ugent.be/scop3p/api/v1"

    #: v1 snake_case -> the column names this toolkit has always used.
    _MODIFICATION_FIELDS = {
        "uniprot_position": "position",
        "modified_residue": "residue",
        "modification_name": "name",
        "source": "source",
        "evidence_terms": "evidence",
        "pubmed": "reference",
        "best_probability": "functionalScore",
    }
    _PEPTIDE_FIELDS = {
        "peptide_sequence": "peptideSequence",
        "peptide_start": "peptideStart",
        "peptide_end": "peptideEnd",
        "peptide_modification_position": "peptideModificationPosition",
        "uniprot_position": "uniprotPosition",
        "score": "score",
    }

    def __init__(self, base_url: str | None = None, timeout: int = 30) -> None:
        self.base_url = (base_url or self.DEFAULT_BASE_URL).rstrip("/")
        self.timeout = timeout

    # -- transport ----------------------------------------------------------

    def _get_json(self, path: str) -> Any:
        url = f"{self.base_url}/{path.lstrip('/')}"
        try:
            response = http_lookup(
                url, headers={"Accept": "application/json"}, logger=LOGGER
            )
        except requests.RequestException as error:
            raise Scop3PApiError(f"Could not reach Scop3P at {url}: {error}") from error

        if response.status_code >= 400:
            raise Scop3PApiError(f"Scop3P returned HTTP {response.status_code} for {url}")

        content_type = (response.headers.get("content-type") or "").split(";")[0].strip()
        if content_type.lower() != "application/json":
            # The SPA catch-all answered, which means this endpoint does not exist.
            raise Scop3PApiError(
                f"Scop3P served {content_type or 'an unknown content type'} rather than "
                f"JSON for {url}. That is the single-page-app catch-all responding, so "
                f"the endpoint has most likely moved; check "
                f"{self.base_url}/openapi.json for the current paths."
            )

        try:
            return response.json()
        except ValueError as error:
            raise Scop3PApiError(f"Scop3P returned malformed JSON for {url}: {error}") from error

    @staticmethod
    def _records(payload: Any, legacy_key: str) -> list[dict]:
        """v1 returns a bare list; the pre-v1 API wrapped it in a dict."""
        if isinstance(payload, list):
            return [item for item in payload if isinstance(item, dict)]
        if isinstance(payload, dict):
            nested = payload.get(legacy_key)
            if isinstance(nested, list):
                return [item for item in nested if isinstance(item, dict)]
        return []

    @classmethod
    def _rename(cls, dataframe: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
        present = {source: target for source, target in mapping.items() if source in dataframe.columns}
        return dataframe.rename(columns=present)

    # -- endpoints ----------------------------------------------------------

    def fetch_peptides_modifications(self, accession: str) -> pd.DataFrame:
        payload = self._get_json(f"proteins/{accession}/peptides")
        dataframe = pd.DataFrame(self._records(payload, "peptides"))
        if dataframe.empty:
            return dataframe

        dataframe = self._rename(dataframe, self._PEPTIDE_FIELDS)

        numeric_columns = [
            "peptideStart",
            "peptideEnd",
            "peptideModificationPosition",
            "uniprotPosition",
        ]
        for column in numeric_columns:
            if column in dataframe.columns:
                dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce").astype("Int64")

        if "score" in dataframe.columns:
            dataframe["score"] = pd.to_numeric(dataframe["score"], errors="coerce")

        # v1 dropped the modifiedResidue field, but it is recoverable: the modification
        # position is 1-based within the peptide, so peptide_start + position - 1 is the
        # UniProt position, and the residue is just that offset into the sequence.
        if "modifiedResidue" not in dataframe.columns:
            dataframe["modifiedResidue"] = [
                self._modified_residue(row) for _, row in dataframe.iterrows()
            ]

        dataframe["label"] = dataframe.apply(self._format_label, axis=1)
        return dataframe

    @memoize(name="scop3p.modifications")
    def fetch_modifications(self, accession: str) -> pd.DataFrame:
        payload = self._get_json(f"proteins/{accession}/modifications")
        dataframe = pd.DataFrame(self._records(payload, "modifications"))
        if dataframe.empty:
            return dataframe

        dataframe = self._rename(dataframe, self._MODIFICATION_FIELDS)

        keep = [
            column
            for column in [
                "position",
                "residue",
                "name",
                "source",
                "evidence",
                "reference",
                "functionalScore",
                "specificSinglyPhosphorylated",
            ]
            if column in dataframe.columns
        ]
        dataframe = dataframe[keep].copy()
        if "position" in dataframe.columns:
            dataframe["position"] = pd.to_numeric(dataframe["position"], errors="coerce").astype("Int64")
        return dataframe

    # -- helpers ------------------------------------------------------------

    @staticmethod
    def _modified_residue(row: pd.Series) -> str:
        sequence = str(row.get("peptideSequence") or "")
        position = row.get("peptideModificationPosition")
        if not sequence or pd.isna(position):
            return ""
        index = int(position) - 1
        return sequence[index] if 0 <= index < len(sequence) else ""

    @staticmethod
    def _format_label(row: pd.Series) -> str:
        peptide_sequence = str(row.get("peptideSequence", ""))
        start = row.get("peptideStart")
        end = row.get("peptideEnd")
        modified_residue = str(row.get("modifiedResidue", "") or "")
        uniprot_position = row.get("uniprotPosition")
        score = row.get("score", "")

        start_txt = str(int(start)) if pd.notna(start) else "?"
        end_txt = str(int(end)) if pd.notna(end) else "?"
        position_txt = str(int(uniprot_position)) if pd.notna(uniprot_position) else "?"

        return (
            f"{peptide_sequence} ({start_txt}-{end_txt}) "
            f"@{modified_residue}{position_txt} score={score}"
        )


class AlphaFoldService:
    """AlphaFold download helper with fallback model versions."""

    def __init__(self, cache_dir: Path, versions: tuple[str, ...] = ("v6", "v4")) -> None:
        self.cache_dir = cache_dir
        self.versions = versions
        self.base_url = "https://alphafold.ebi.ac.uk/files"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def download_pdb(self, accession: str) -> Path:
        out_path = self.cache_dir / f"{accession}.pdb"
        last_error: Exception | None = None

        for version in self.versions:
            filename = f"AF-{accession}-F1-model_{version}.pdb"
            url = f"{self.base_url}/{filename}"
            try:
                urllib.request.urlretrieve(url, out_path)
                if out_path.stat().st_size < 1000:
                    raise RuntimeError(f"Downloaded file too small from {url}")
                return out_path
            except (HTTPError, URLError, RuntimeError) as error:
                last_error = error

        raise RuntimeError(
            f"Could not download AlphaFold PDB for {accession}. "
            f"Last error: {last_error}"
        )
