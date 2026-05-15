from __future__ import annotations

from pathlib import Path
from urllib.error import HTTPError, URLError
import urllib.request

import pandas as pd
from scop3p_api_client.api import Scop3pRestApi


class Scop3PClient:
    """API client for Scop3P endpoints used by converted apps."""

    def __init__(self, base_url: str = "https://iomics.ugent.be/scop3p/api", timeout: int = 30) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.api = Scop3pRestApi(default_timeout=timeout)

    def fetch_peptides_modifications(self, accession: str) -> pd.DataFrame:
        payload = self.api.fetch_peptides(accession)
        df = pd.DataFrame(payload.get("peptides", []))
        if df.empty:
            return df

        numeric_columns = ["peptideStart", "peptideEnd", "peptideModificationPosition", "uniprotPosition"]
        for column in numeric_columns:
            if column in df.columns:
                df[column] = pd.to_numeric(df[column], errors="coerce").astype("Int64")

        if "score" in df.columns:
            df["score"] = pd.to_numeric(df["score"], errors="coerce")

        df["label"] = df.apply(self._format_label, axis=1)
        return df

    def fetch_modifications(self, accession: str) -> pd.DataFrame:
        payload = self.api.fetch_modifications(accession)
        if isinstance(payload, list):
            payload = payload[0] if payload else {}

        dataframe = pd.DataFrame(payload.get("modifications", []))
        if dataframe.empty:
            return dataframe

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

    @staticmethod
    def _format_label(row: pd.Series) -> str:
        peptide_sequence = str(row.get("peptideSequence", ""))
        start = row.get("peptideStart")
        end = row.get("peptideEnd")
        modified_residue = str(row.get("modifiedResidue", ""))
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
