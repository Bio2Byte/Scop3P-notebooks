from __future__ import annotations

from pathlib import Path
from urllib.error import HTTPError, URLError
import urllib.request

import pandas as pd
import requests


class Scop3PClient:
    """API client for Scop3P endpoints used by converted apps."""

    def __init__(self, base_url: str = "https://iomics.ugent.be/scop3p/api", timeout: int = 30) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def fetch_peptides_modifications(self, accession: str) -> pd.DataFrame:
        url = f"{self.base_url}/get-peptides-modifications"
        response = requests.get(url, params={"accession": accession}, timeout=self.timeout)
        response.raise_for_status()
        payload = response.json()

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
