from __future__ import annotations

import colorsys
import html
import json
import os
import subprocess
import tempfile
import uuid
from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import re

import pandas as pd
import pyvis.network
import requests
from Bio.PDB import PDBIO, PDBParser, Select
from b2bTools import SingleSeq, constants
from scipy.spatial import KDTree

from common import http_lookup as http_lookup_module
from common.cache import (
    PREDICTION_MAX_ENTRIES,
    PREDICTION_TTL_SECONDS,
    memoize,
    shared_structure_dir,
    structure_file_lock,
)
from common.http_lookup import lookup as http_lookup
from common.logging_utils import get_logger, quiet_third_party
from common.structure_labels import (
    CHOOSE_ENTRY_PLACEHOLDER,
    LOOKUP_FAILED_PLACEHOLDER,
    NO_STRUCTURES_PLACEHOLDER,
    chain_label,
    structure_option_label,
)

from .services import Scop3PClient

LOGGER = get_logger("scop3p.common.structure_viz")


B2B_METRIC_COLUMNS = (
    "backbone",
    "sidechain",
    "ppII",
    "coil",
    "sheet",
    "helix",
    "earlyFolding",
    "disoMine",
)
B2B_NORMALIZED_SUFFIX = "_normalized"


class StructureVizService:
    #: Annotation-lookup timeouts, from the policy shared with every other protocol.
    #: File downloads deliberately keep the longer ``self.timeout``.
    CONNECT_TIMEOUT = http_lookup_module.CONNECT_TIMEOUT
    LOOKUP_READ_TIMEOUT = http_lookup_module.READ_TIMEOUT

    def __init__(self, workdir: Path, timeout: int = 60) -> None:
        self.workdir = workdir
        self.timeout = timeout
        self.scop3p_client = Scop3PClient(timeout=timeout)
        self.workdir.mkdir(parents=True, exist_ok=True)
        # Resolving a mapping costs a PDBe round trip, and every tab that draws a site
        # wants the same one. Keyed on (entry, accession, chain).
        self._mapping_cache: dict[tuple[str, str, str], PositionMapping] = {}

    @memoize(name="scop3p.modifications.table")
    def fetch_ptms(self, accession: str) -> pd.DataFrame:
        dataframe = self.scop3p_client.fetch_modifications(accession)
        if dataframe.empty:
            return dataframe
        return dataframe

    @memoize(name="uniprot.ptm.features")
    def fetch_uniprot_ptms(self, accession: str) -> pd.DataFrame:
        """Single-residue PTM features from the EBI Proteins API.

        ``categories=PTM`` without a type filter, so phosphorylation, acetylation,
        methylation, glycosylation, lipidation and the rest all come through when
        annotated. Only single-residue features are kept (``begin == end``) because the
        rest of the app maps PTMs onto residue positions. Descriptions are trimmed at
        the first semicolon: "Phosphotyrosine; by autocatalysis" becomes
        "Phosphotyrosine".
        """
        response = self._lookup(
            f"{PROTEINS_FEATURES_URL}/{accession}",
            headers={"Accept": "application/json"},
            params=[("categories", "PTM")],
        )
        response.raise_for_status()
        payload = response.json()

        sequence = payload.get("sequence") or ""
        rows: list[dict[str, object]] = []
        for feature in payload.get("features", []) or []:
            if feature.get("category") != "PTM":
                continue
            begin, end = feature.get("begin"), feature.get("end")
            if begin is None or end is None or str(begin) != str(end):
                continue
            try:
                position = int(begin)
            except (TypeError, ValueError):
                continue

            description = feature.get("description") or feature.get("type") or "PTM"
            name = str(description).split(";", 1)[0].strip()
            evidence, reference = _format_uniprot_evidence(feature.get("evidences"))

            rows.append(
                {
                    "position": position,
                    "residue": residue_three_letter(name, sequence, position),
                    "modification": name,
                    "name": name,
                    "evidence": evidence,
                    "source": "UniProt",
                    "reference": reference,
                    "functionalScore": pd.NA,
                    "feature_type": feature.get("type"),
                    "ACC_ID": accession,
                }
            )

        frame = pd.DataFrame(rows, columns=list(PTM_TABLE_COLUMNS))
        if not frame.empty:
            frame["position"] = pd.to_numeric(frame["position"], errors="coerce").astype("Int64")
        return frame

    def fetch_pdbe_sifts_map(
        self, pdb_id: str, accession: str, chain: str | None = None
    ) -> dict[int, int]:
        """Tier 1: the PDBe SIFTS API's residue-level mapping for one chain."""
        response = self._lookup(
            f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}",
            headers={"Accept": "application/json"},
        )
        response.raise_for_status()
        return parse_pdbe_uniprot_mappings(response.json(), pdb_id, accession, chain)

    def download_pdbe_updated_cif(self, pdb_id: str) -> Path | None:
        """The SIFTS-enriched mmCIF, falling back to RCSB's plain one."""
        pdb_id = (pdb_id or "").strip().lower()
        if not pdb_id:
            return None
        target = shared_structure_dir() / f"{pdb_id}_updated.cif"
        if target.exists() and target.stat().st_size > 1000:
            return target
        for url in (
            f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_updated.cif",
            f"https://files.rcsb.org/download/{pdb_id.upper()}.cif",
        ):
            try:
                response = requests.get(url, timeout=(self.CONNECT_TIMEOUT, self.timeout))
            except requests.RequestException:
                continue
            if response.ok and response.text and len(response.text) > 1000:
                target.write_text(response.text, encoding="utf-8")
                return target
        return None

    def build_position_mapping(
        self,
        pdb_id: str,
        accession: str,
        chain: str | None = None,
        uniprot_range: tuple[int, int] | None = None,
        pdb_path: Path | None = None,
    ) -> PositionMapping:
        """Resolve PDB author numbering against UniProt numbering, best source first.

        Never raises: a mapping is always returned, and its ``source`` says how much to
        trust it. Each tier's failure is logged rather than swallowed, because "SIFTS was
        unreachable" and "SIFTS has no mapping for this chain" lead to the same degraded
        result but are different problems.
        """
        chain_key = (chain or "").strip().upper()[:1]
        cache_key = (pdb_id.upper(), (accession or "").upper(), chain_key)
        cached = self._mapping_cache.get(cache_key)
        if cached is not None:
            return cached

        mapping: dict[int, int] = {}
        source = "direct"

        if pdb_id and accession:
            try:
                mapping = self.fetch_pdbe_sifts_map(pdb_id, accession, chain_key or None)
                if mapping:
                    source = "sifts-api"
            except Exception as error:  # noqa: BLE001
                LOGGER.warning(
                    "sifts api lookup failed pdb=%s chain=%s error=%s",
                    pdb_id, chain_key or "-", error,
                    extra={"event": "position_mapping"},
                )

            if not mapping:
                try:
                    cif_path = self.download_pdbe_updated_cif(pdb_id)
                    if cif_path is not None:
                        mapping = parse_sifts_map_from_cif(
                            cif_path.read_text(encoding="utf-8", errors="ignore"),
                            accession,
                            chain_key or None,
                        )
                        if mapping:
                            source = "sifts-mmcif"
                except Exception as error:  # noqa: BLE001
                    LOGGER.warning(
                        "sifts mmcif fallback failed pdb=%s error=%s", pdb_id, error,
                        extra={"event": "position_mapping"},
                    )

        if not mapping and uniprot_range and pdb_path is not None:
            try:
                observed = StructureOps.chain_range_from_pdb(Path(pdb_path), chain_key or "A")
                mapping = offset_mapping(uniprot_range, observed)
                if mapping:
                    source = "chain-range-offset"
            except Exception as error:  # noqa: BLE001
                LOGGER.warning(
                    "chain-range offset failed pdb=%s error=%s", pdb_id, error,
                    extra={"event": "position_mapping"},
                )

        resolved = PositionMapping(
            pdb_to_uniprot=mapping,
            uniprot_to_pdb=_invert(mapping),
            source=source,
            pdb_id=(pdb_id or "").upper(),
            chain=chain_key,
        )
        LOGGER.info(
            "position mapping resolved pdb=%s chain=%s source=%s residues=%s",
            resolved.pdb_id or "-", resolved.chain or "-", resolved.source, len(mapping),
            extra={"event": "position_mapping"},
        )
        self._mapping_cache[cache_key] = resolved
        return resolved

    @memoize(name="uniprot.pdb.xrefs")
    def fetch_pdb_xrefs(self, accession: str) -> list[PdbXref]:
        """PDB entries cross-referenced from this accession's UniProt entry.

        Uses a short connect timeout rather than the service default. This call sits on
        the critical path of "Set protein", and a synchronous effect blocks the ASGI
        loop for every connected session while it waits -- a slow UniProt froze the app
        for the full 60 seconds before this was bounded.
        """
        response = self._lookup(
            f"https://rest.uniprot.org/uniprotkb/{accession}.json",
            headers={"Accept": "application/json"},
        )
        response.raise_for_status()
        return parse_pdb_xrefs(response.json())

    @memoize(name="uniprot.disease.variants")
    def fetch_variants(self, accession: str) -> pd.DataFrame:
        url = f"https://www.ebi.ac.uk/proteins/api/variation/{accession}"
        response = self._lookup(url, headers={"Accept": "application/json"})
        response.raise_for_status()
        payload = response.json()
        rows: list[dict[str, object]] = []
        for feature in payload.get("features", []):
            if feature.get("type") != "VARIANT":
                continue
            for association in feature.get("association", []):
                if association.get("disease") is not True:
                    continue
                try:
                    position = int(feature.get("begin"))
                except Exception:
                    position = None
                rows.append(
                    {
                        "ACC_ID": accession,
                        "position": position,
                        "WT": feature.get("wildType"),
                        "MT": feature.get("mutatedType"),
                        "consequence": feature.get("consequenceType"),
                        "disease_name": association.get("name"),
                    }
                )
        return pd.DataFrame(rows)

    @memoize(name="uniprot.sequence.fasta")
    def fetch_sequence(self, accession: str) -> str:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        response = self._lookup(url)
        response.raise_for_status()
        lines = response.text.splitlines()
        return "".join(line.strip() for line in lines if line and not line.startswith(">"))

    @memoize(
        name="b2b.prediction.table",
        ttl_seconds=PREDICTION_TTL_SECONDS,
        max_entries=PREDICTION_MAX_ENTRIES,
    )
    def predict_b2b(self, accession: str, sequence: str) -> pd.DataFrame:
        """Run Bio2Byte and normalise the result.

        Cached on **(accession, sequence)**, and the sequence part is not optional: the
        Mutation Effect protocol predicts a mutated sequence under the same accession, so
        an accession-only key would hand back the wild-type prediction for the mutant and
        silently destroy the comparison that protocol exists to make.
        """
        with tempfile.NamedTemporaryFile(prefix="seq_", suffix=".fasta", mode="w") as fasta_file:
            fasta_file.write(f">{accession}\n{sequence}\n")
            fasta_file.flush()
            tools = []
            for name in ["TOOL_DYNAMINE", "TOOL_DISOMINE", "TOOL_EFOLDMINE"]:
                if hasattr(constants, name):
                    tools.append(getattr(constants, name))
            # b2bTools reports progress with print() and its dependencies warn on every
            # prediction; both are captured to DEBUG so they do not bury the trail.
            with quiet_third_party(LOGGER, event="b2b_predict"):
                predictor = SingleSeq(fasta_file.name)
                prediction = (
                    predictor.predict(tools=tools).get_all_predictions()
                    if tools
                    else predictor.predict().get_all_predictions()
                )
        
        protein = prediction.get("proteins", {}).get(accession, {})
        return self._normalize_b2b_prediction(protein)

    @staticmethod
    def _normalize_b2b_prediction(protein: dict[str, object]) -> pd.DataFrame:
        sequence = "".join(protein.get("seq", ""))
        size = (
            len(sequence)
            or len(protein.get("backbone", []))
            or len(protein.get("sidechain", []))
            or len(protein.get("ppII", []))
            or len(protein.get("coil", []))
            or len(protein.get("sheet", []))
            or len(protein.get("helix", []))
            or len(protein.get("earlyFolding", []))
            or len(protein.get("disoMine", []))
        )

        def _coerce_series(value: object) -> list[object]:
            if isinstance(value, (list, tuple)):
                values = list(value)
            elif hasattr(value, "tolist"):
                values = list(value.tolist())  # type: ignore[call-arg]
            else:
                values = []

            if len(values) < size:
                values.extend([None] * (size - len(values)))
            elif len(values) > size:
                values = values[:size]
            return values

        dataframe = pd.DataFrame(
            {
                "Position": list(range(1, size + 1)),
                "Amino acid": protein.get("seq"),
                "backbone": _coerce_series(protein.get("backbone")),
                "sidechain": _coerce_series(protein.get("sidechain")),
                "ppII": _coerce_series(protein.get("ppII")),
                "coil": _coerce_series(protein.get("coil")),
                "sheet": _coerce_series(protein.get("sheet")),
                "helix": _coerce_series(protein.get("helix")),
                "earlyFolding": _coerce_series(protein.get("earlyFolding")),
                "disoMine": _coerce_series(protein.get("disoMine")),
            }
        )
        for column in B2B_METRIC_COLUMNS:
            dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce")
            dataframe[StructureVizService.b2b_metric_column(column, normalized=True)] = (
                StructureVizService._min_max_normalize_series(dataframe[column])
            )
        return dataframe

    @staticmethod
    def b2b_metric_column(metric: str, *, normalized: bool = False) -> str:
        return f"{metric}{B2B_NORMALIZED_SUFFIX}" if normalized else metric

    @staticmethod
    def _min_max_normalize_series(series: pd.Series) -> pd.Series:
        numeric = pd.to_numeric(series, errors="coerce")
        non_null = numeric.dropna()
        if non_null.empty:
            return pd.Series([pd.NA] * len(numeric), index=numeric.index, dtype="Float64")

        minimum = float(non_null.min())
        maximum = float(non_null.max())
        if minimum == maximum:
            return pd.Series(
                [0.0 if pd.notna(value) else pd.NA for value in numeric],
                index=numeric.index,
                dtype="Float64",
            )

        normalized = (numeric - minimum) / (maximum - minimum)
        return normalized.astype("Float64")

    def download_alphafold_pdb(self, accession: str) -> Path:
        out_path = self.workdir / f"AF-{accession}-F1-model_v6.pdb"
        url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.pdb"
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        out_path.write_bytes(response.content)
        return out_path

    def _lookup(self, url: str, **kwargs) -> requests.Response:
        """An annotation lookup under the shared HTTP policy.

        See common.http_lookup for why the timeouts are what they are. Kept as a method
        so call sites read naturally and so the logger is attached automatically.
        """
        return http_lookup(url, logger=LOGGER, **kwargs)

    def pdb_path_for(self, pdb_id: str) -> Path | None:
        """Where download_pdb would put this entry, so a caller can tell whether the
        structure in hand came from that entry or from an upload.

        In the shared structure directory, not this session's workdir: the file is an
        immutable upstream artefact and every session wants the same one.
        """
        pdb_key = (pdb_id or "").strip().upper()
        return shared_structure_dir() / f"{pdb_key}.pdb" if pdb_key else None

    def download_pdb(self, pdb_id: str) -> Path:
        """Fetch a PDB entry once per process, not once per session."""
        out_path = self.pdb_path_for(pdb_id)
        if out_path is None:
            raise ValueError("a PDB ID is required")
        pdb_key = pdb_id.strip().upper()
        if out_path.exists() and out_path.stat().st_size > 0:
            return out_path

        # One lock per file: two tabs asking for the same entry at the same moment would
        # otherwise both download it, and a reader could see a half-written file.
        with structure_file_lock(out_path.name):
            if out_path.exists() and out_path.stat().st_size > 0:
                return out_path
            url = f"https://files.rcsb.org/download/{pdb_key}.pdb"
            response = requests.get(url, timeout=self.timeout)
            response.raise_for_status()
            out_path.write_bytes(response.content)
            return out_path

    def resolve_uploaded_or_remote_pdb(
        self,
        upload: list[dict[str, object]] | None,
        pdb_id: str | None = None,
        *,
        target_name: str | None = None,
    ) -> Path | None:
        if upload:
            upload_row = upload[0]
            source = Path(str(upload_row["datapath"]))
            filename = target_name or Path(str(upload_row["name"])).name
            target = self.workdir / filename
            target.write_bytes(source.read_bytes())
            return target

        if pdb_id and pdb_id.strip():
            downloaded = self.download_pdb(pdb_id)
            if target_name is None or downloaded.name == target_name:
                return downloaded
            target = self.workdir / target_name
            target.write_bytes(downloaded.read_bytes())
            return target

        return None


PROTEINS_FEATURES_URL = "https://www.ebi.ac.uk/proteins/api/features"


def _format_uniprot_evidence(evidences: object) -> tuple[str, str]:
    """Condense UniProt feature evidence into (codes, literature references)."""
    codes: list[str] = []
    references: list[str] = []
    for evidence in evidences or []:
        if not isinstance(evidence, dict):
            continue
        code = evidence.get("code")
        if code:
            codes.append(str(code))
        source = evidence.get("source") or {}
        name, identifier = source.get("name"), source.get("id")
        if name and identifier:
            references.append(f"{name}:{identifier}")
        elif identifier:
            references.append(str(identifier))
    return "; ".join(dict.fromkeys(codes)), "; ".join(dict.fromkeys(references))


# ---------------------------------------------------------------------------
# PTM tables: Scop3P and UniProt, on one schema
#
# Ported from notebooks/Scop3P_PTM_structure_viz_voila_app.ipynb, which grew the
# UniProt PTM source and the site-level merge after the first Shiny conversion.
# ---------------------------------------------------------------------------

#: Columns every PTM table carries, whatever its source.
PTM_TABLE_COLUMNS = (
    "position",
    "residue",
    "modification",
    "name",
    "evidence",
    "source",
    "reference",
    "functionalScore",
    # specificSinglyPhosphorylated is intentionally absent: the pre-v1 Scop3P payload
    # carried it, the v1 payload does not, and UniProt never did, so it could only ever
    # be an empty column.
    "feature_type",
    "ACC_ID",
)

_AA1_TO_AA3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL", "U": "SEC", "O": "PYL",
}

#: Scop3P names a modified residue descriptively; UniProt features and the 3D viewer's
#: colour map both key on the three-letter code. Without this mapping every Scop3P PTM
#: fell through to the viewer's default colour, and the site-level merge below could
#: never match a Scop3P row against a UniProt one.
_MODIFIED_RESIDUE_TO_CODE = {
    "phosphoserine": "SER",
    "phosphothreonine": "THR",
    "phosphotyrosine": "TYR",
    "phosphohistidine": "HIS",
    "phosphoaspartate": "ASP",
    "phosphoglutamate": "GLU",
    "n6-acetyllysine": "LYS",
    "omega-n-methylarginine": "ARG",
    "symmetric dimethylarginine": "ARG",
    "asymmetric dimethylarginine": "ARG",
    "n6,n6-dimethyllysine": "LYS",
    "n6-methyllysine": "LYS",
    "s-nitrosocysteine": "CYS",
}


def residue_three_letter(value: object, sequence: str = "", position: object = None) -> str:
    """Best three-letter residue code for a PTM record.

    Tries, in order: a known descriptive modified-residue name ("Phosphoserine" ->
    "SER"), an already-valid three-letter code, a one-letter code, then the residue at
    ``position`` in ``sequence``. Returns "" when none of those apply, which the callers
    treat as "unknown" rather than guessing.
    """
    text = str(value or "").strip()
    lowered = text.lower()

    if lowered in _MODIFIED_RESIDUE_TO_CODE:
        return _MODIFIED_RESIDUE_TO_CODE[lowered]

    upper = text.upper()
    if len(upper) == 3 and upper in set(_AA1_TO_AA3.values()):
        return upper
    if len(upper) == 1 and upper in _AA1_TO_AA3:
        return _AA1_TO_AA3[upper]

    try:
        index = int(position)
    except (TypeError, ValueError):
        return ""
    if sequence and 1 <= index <= len(sequence):
        letter = sequence[index - 1].upper()
        return _AA1_TO_AA3.get(letter, letter)
    return ""


def build_ptm_table(
    dataframe: pd.DataFrame,
    accession: str,
    *,
    source: str = "Scop3P",
    sequence: str = "",
) -> pd.DataFrame:
    """Put a Scop3P modification frame on :data:`PTM_TABLE_COLUMNS`.

    ``residue`` becomes the three-letter code so the 3D viewer can colour it and so it
    can be matched against UniProt; the descriptive Scop3P text is preserved in
    ``modification`` rather than discarded.
    """
    if dataframe is None or dataframe.empty:
        return pd.DataFrame(columns=list(PTM_TABLE_COLUMNS))

    out = dataframe.copy()
    out["modification"] = out["residue"] if "residue" in out.columns else pd.NA
    out["residue"] = [
        residue_three_letter(row.get("residue"), sequence, row.get("position"))
        for _, row in out.iterrows()
    ]
    out["ACC_ID"] = accession
    out["feature_type"] = source
    if "source" not in out.columns:
        out["source"] = source
    else:
        out["source"] = out["source"].fillna(source)
        out.loc[out["source"].astype(str).str.strip().eq(""), "source"] = source

    for column in PTM_TABLE_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    out["position"] = pd.to_numeric(out["position"], errors="coerce").astype("Int64")
    return out[list(PTM_TABLE_COLUMNS)].reset_index(drop=True)


def merge_ptm_tables(scop3p: pd.DataFrame, uniprot: pd.DataFrame) -> pd.DataFrame:
    """Merge two PTM tables without listing the same residue site twice.

    Site identity is accession + three-letter residue + position. Where both sources
    describe the same site the Scop3P row wins -- its naming and evidence are the more
    specific -- and the UniProt reference is folded in, so nothing is lost. UniProt-only
    sites are appended.
    """
    columns = list(PTM_TABLE_COLUMNS)
    empty = pd.DataFrame(columns=columns)
    left = empty if scop3p is None or scop3p.empty else scop3p.copy()
    right = empty if uniprot is None or uniprot.empty else uniprot.copy()

    if left.empty and right.empty:
        return empty
    if left.empty:
        return right[columns].sort_values("position", na_position="last").reset_index(drop=True)
    if right.empty:
        return left[columns].sort_values("position", na_position="last").reset_index(drop=True)

    def _key(frame: pd.DataFrame) -> list[tuple]:
        return [
            (str(row.get("ACC_ID")), str(row.get("residue")), str(row.get("position")))
            for _, row in frame.iterrows()
        ]

    left = left.reset_index(drop=True)
    right = right.reset_index(drop=True)
    left_keys = _key(left)
    right_keys = _key(right)
    index_by_key = {key: position for position, key in enumerate(left_keys)}

    for offset, key in enumerate(right_keys):
        target = index_by_key.get(key)
        if target is None:
            continue
        # Scop3P cites a bare PMID, UniProt cites "PubMed:<pmid>". Compare on the bare
        # id so the same paper is not listed twice in two notations.
        seen: dict[str, str] = {}
        for value in (left.at[target, "reference"], right.at[offset, "reference"]):
            if pd.isna(value):
                continue
            for part in re.split(r"[;,]", str(value)):
                citation = part.strip()
                if not citation:
                    continue
                seen.setdefault(citation.split(":")[-1].strip(), citation)
        merged_reference = "; ".join(seen.values())
        if merged_reference:
            left.at[target, "reference"] = merged_reference
        if pd.isna(left.at[target, "evidence"]) or not str(left.at[target, "evidence"]).strip():
            left.at[target, "evidence"] = right.at[offset, "evidence"]
        # Record that both sources saw this site. Split on both separators: Scop3P
        # already comma-joins its own sources ("UniProt, PRIDE"), so splitting on ";"
        # alone produced "UniProt, PRIDE; UniProt".
        left.at[target, "source"] = ", ".join(
            dict.fromkeys(
                part.strip()
                for value in (left.at[target, "source"], right.at[offset, "source"])
                for part in re.split(r"[;,]", str(value))
                if part.strip()
            )
        )

    unmatched = right[[key not in index_by_key for key in right_keys]]
    out = pd.concat([left, unmatched], ignore_index=True)[columns]
    return out.sort_values(["position", "source"], na_position="last").reset_index(drop=True)


# ---------------------------------------------------------------------------
# PDB entries for an accession
#
# Ported from notebooks/Scop3P_PTM_structure_viz_voila_app.ipynb. Choosing a PDB entry
# from a dropdown of the accession's own cross-references (rather than typing an ID)
# brings two things the text field could not: the chains that actually contain this
# protein, and the UniProt residue range each of those chains covers. The range is what
# lets a viewer or a RIN be built on the mapped region instead of the whole asymmetric
# unit.
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class PdbXref:
    """One PDB entry cross-referenced from a UniProt entry."""

    pdb_id: str
    chain_ranges: dict[str, tuple[int, int] | None]
    method: str = ""
    resolution: str = ""

    def label(self) -> str:
        """Delegated to common.structure_labels so all protocols read alike."""
        return structure_option_label(
            self.pdb_id,
            method=self.method,
            resolution=self.resolution,
            chains=self.chain_ranges,
        )


def parse_pdb_xrefs(payload: dict) -> list[PdbXref]:
    """Pull PDB cross-references, with per-chain UniProt ranges, out of a UniProt entry.

    UniProt writes the Chains property in several shapes -- ``A=1-200``,
    ``A=1-200, B=5-150``, ``A/B=1-120``, and ``A=10-50; 70-120`` for a chain covering
    two discontinuous stretches. The last of those is why the range is taken as
    min(starts) to max(ends) rather than the first pair.
    """
    refs: list[PdbXref] = []
    for xref in payload.get("uniProtKBCrossReferences", []) or []:
        if xref.get("database") != "PDB":
            continue
        pdb_id = (xref.get("id") or "").upper().strip()
        if not pdb_id:
            continue

        properties = {
            item.get("key"): item.get("value")
            for item in (xref.get("properties") or [])
            if isinstance(item, dict)
        }
        chains_raw = properties.get("Chains") or properties.get("Chain") or ""

        chain_ranges: dict[str, tuple[int, int] | None] = {}
        # Split on a comma only where the next token starts a new "<chains>=" group, so
        # "A=10-50; 70-120" stays with its chain.
        for part in [p.strip() for p in re.split(r",\s*(?=[A-Za-z0-9/]+=)", chains_raw) if p.strip()]:
            if "=" not in part:
                continue
            chains_text, ranges_text = part.split("=", 1)
            starts: list[int] = []
            ends: list[int] = []
            for match in re.finditer(r"(\d+)\s*-\s*(\d+)", ranges_text):
                starts.append(int(match.group(1)))
                ends.append(int(match.group(2)))
            chains = [c.strip().upper()[:1] for c in re.split(r"[/\s]+", chains_text.strip()) if c.strip()]
            for chain in chains:
                if starts and ends:
                    chain_ranges[chain] = (min(starts), max(ends))
                else:
                    chain_ranges.setdefault(chain, None)

        refs.append(
            PdbXref(
                pdb_id=pdb_id,
                chain_ranges=chain_ranges,
                method=str(properties.get("Method") or "").strip(),
                resolution=str(properties.get("Resolution") or "").strip(),
            )
        )
    refs.sort(key=lambda ref: ref.pdb_id)
    return refs


def pdb_entry_choices(
    refs: list[PdbXref], *, lookup_failed: bool = False
) -> dict[str, str]:
    """``{pdb_id: label}`` for a select, with a blank first entry.

    ``lookup_failed`` is separate from "no refs" on purpose. An empty picker after a
    failed request looked identical to a protein with no structures, so the user was told
    a falsehood about their protein and had no reason to press the button again.
    """
    if lookup_failed:
        return {"": LOOKUP_FAILED_PLACEHOLDER}
    if not refs:
        return {"": NO_STRUCTURES_PLACEHOLDER}
    choices = {"": CHOOSE_ENTRY_PLACEHOLDER}
    for ref in refs:
        choices[ref.pdb_id] = ref.label()
    return choices


def chain_choices_for_pdb(refs: list[PdbXref], pdb_id: str) -> dict[str, str]:
    """``{chain: "A (1-200)"}`` for the chains of one entry that hold this protein."""
    wanted = (pdb_id or "").upper().strip()
    ranges: dict[str, tuple[int, int] | None] = {}
    for ref in refs:
        if ref.pdb_id.upper() == wanted:
            ranges.update(ref.chain_ranges)
    return {chain: chain_label(chain, ranges[chain]) for chain in sorted(ranges)}


def uniprot_range_for_chain(
    refs: list[PdbXref], pdb_id: str, chain: str
) -> tuple[int, int] | None:
    """The UniProt residue range one chain of one entry covers, if UniProt states it."""
    wanted_pdb = (pdb_id or "").upper().strip()
    wanted_chain = (chain or "").upper().strip()[:1]
    if not wanted_pdb or not wanted_chain:
        return None
    for ref in refs:
        if ref.pdb_id.upper() != wanted_pdb:
            continue
        span = ref.chain_ranges.get(wanted_chain)
        if span:
            return int(span[0]), int(span[1])
    return None


# ---------------------------------------------------------------------------
# SIFTS residue-level numbering
# ---------------------------------------------------------------------------
# A PDB entry numbers its residues however the depositors chose (author numbering).
# UniProt numbers the canonical sequence from 1. The two frequently disagree -- 2IVT
# numbers the RET kinase domain 705-1013 to match UniProt, but plenty of entries start
# at 1, renumber a construct, or carry expression tags -- so painting a UniProt-numbered
# PTM straight onto a PDB structure can silently mark the wrong residue.
#
# SIFTS (Structure Integration with Function, Taxonomy and Sequence) is the EBI's
# residue-level correspondence between the two, and is the authority here. This ports the
# notebook's four-tier resolution, best first:
#
#   1. PDBe SIFTS API      -- residue-level mapping segments, the authoritative answer.
#   2. SIFTS-enriched mmCIF -- PDBe's "updated" mmCIF carries the same correspondence per
#                              atom in _atom_site.pdbx_sifts_xref_db_num.
#   3. Chain-range offset   -- no SIFTS available: line up UniProt's stated range for the
#                              chain against the residue numbers actually present.
#   4. Direct numbering     -- assume they agree. Correct for AlphaFold models, which are
#                              built on the UniProt sequence, and a guess anywhere else.
#
# Which tier answered is carried on the result and shown to the user, because "the marks
# are where SIFTS says" and "the marks assume the numbering agrees" are very different
# claims to make about a figure.

#: How a mapping was obtained, worst to best, with the wording shown in the UI.
MAPPING_SOURCE_LABELS: dict[str, str] = {
    "sifts-api": "SIFTS residue-level mapping (PDBe API)",
    "sifts-mmcif": "SIFTS numbering from PDBe's enriched mmCIF",
    "chain-range-offset": "offset inferred from UniProt's chain range (no SIFTS available)",
    "direct": "direct numbering, UniProt positions used as-is",
}

#: Tiers that came from SIFTS proper, as opposed to being inferred or assumed.
SIFTS_SOURCES = frozenset({"sifts-api", "sifts-mmcif"})


@dataclass(frozen=True, slots=True)
class PositionMapping:
    """A residue-number correspondence between one PDB chain and a UniProt sequence.

    ``uniprot_to_pdb`` maps one UniProt position to a *list* of author numbers: a
    homodimer with the protein in two chains legitimately has two structural residues for
    one sequence position, and both should be marked.
    """

    pdb_to_uniprot: dict[int, int]
    uniprot_to_pdb: dict[int, list[int]]
    source: str
    pdb_id: str = ""
    chain: str = ""

    @property
    def is_sifts(self) -> bool:
        return self.source in SIFTS_SOURCES

    @property
    def is_identity(self) -> bool:
        """True when the mapping asserts the two numberings agree."""
        return self.source == "direct"

    def to_pdb(self, uniprot_position: int) -> list[int]:
        """The author residue number(s) for a UniProt position."""
        if self.is_identity:
            return [int(uniprot_position)]
        return list(self.uniprot_to_pdb.get(int(uniprot_position), []))

    def describe(self) -> str:
        """One sentence naming the method, for the status line and the viewer panel."""
        label = MAPPING_SOURCE_LABELS.get(self.source, self.source)
        if self.is_identity:
            return f"Positions: {label}."
        where = f"{self.pdb_id} chain {self.chain}".strip()
        return (
            f"Positions mapped via {label}"
            f"{f' for {where}' if where.strip() else ''}: "
            f"{len(self.pdb_to_uniprot)} residues aligned to UniProt numbering."
        )


def identity_mapping(pdb_id: str = "", chain: str = "") -> PositionMapping:
    """The "numbering already agrees" mapping, correct for AlphaFold models."""
    return PositionMapping({}, {}, "direct", pdb_id or "", chain or "")


def _invert(pdb_to_uniprot: dict[int, int]) -> dict[int, list[int]]:
    inverted: dict[int, list[int]] = {}
    for pdb_resi, uniprot_position in pdb_to_uniprot.items():
        inverted.setdefault(int(uniprot_position), []).append(int(pdb_resi))
    for numbers in inverted.values():
        numbers.sort()
    return inverted


def _same_accession(left: str, right: str) -> bool:
    """Compare accessions ignoring the isoform suffix (``P07949-2`` is ``P07949``)."""
    return (left or "").upper().split("-")[0] == (right or "").upper().split("-")[0]


def _author_number(bound: dict) -> int | None:
    """The author residue number at one end of a SIFTS segment, or None if absent.

    Deliberately does not fall back to ``residue_number``: that is label numbering, a
    different coordinate system.
    """
    try:
        return int(bound.get("author_residue_number"))
    except (TypeError, ValueError):
        return None


def parse_pdbe_uniprot_mappings(
    payload: dict, pdb_id: str, accession: str, chain: str | None = None
) -> dict[int, int]:
    """Expand the PDBe SIFTS API's mapping segments into per-residue numbers.

    Each segment states a UniProt span and the author-numbered span it corresponds to.
    Both are walked together one residue at a time; the shorter of the two bounds the
    walk, so a malformed segment truncates instead of inventing positions. Segments can
    run backwards, which is why the step is signed rather than assumed to be +1.
    """
    wanted_pdb = (pdb_id or "").strip()
    wanted_chain = (chain or "").strip().upper()[:1] if chain else None

    entry: dict = {}
    for key in (wanted_pdb, wanted_pdb.lower(), wanted_pdb.upper()):
        if key and isinstance(payload.get(key), dict):
            entry = payload[key]
            break
    uniprot_block = entry.get("UniProt") or entry.get("uniprot") or {}

    # Prefer the accession asked for; fall back to its isoforms, then to anything present.
    keys = [key for key in uniprot_block if key == accession]
    keys += [
        key
        for key in uniprot_block
        if key not in keys and _same_accession(str(key), accession)
    ]
    keys = keys or list(uniprot_block)

    mapping: dict[int, int] = {}
    for key in keys:
        for segment in (uniprot_block.get(key) or {}).get("mappings") or []:
            segment_chain = (
                segment.get("chain_id")
                or segment.get("author_chain_id")
                or segment.get("struct_asym_id")
                or ""
            ).strip().upper()[:1]
            if wanted_chain and segment_chain and segment_chain != wanted_chain:
                continue

            start, end = segment.get("start") or {}, segment.get("end") or {}
            try:
                uniprot_start = int(segment.get("unp_start") or segment.get("uniprot_start"))
                uniprot_end = int(segment.get("unp_end") or segment.get("uniprot_end"))
            except (TypeError, ValueError):
                continue

            # PDBe frequently reports author_residue_number as null at one end of a
            # segment -- 3 of the 5 segments across 2IVT and 1A3N do. The other end plus
            # the UniProt span is enough to reconstruct it, because a SIFTS segment is
            # colinear by construction. What must NOT be substituted is "residue_number":
            # that is label (entity) numbering, which for 2IVT runs 4-314 against author
            # numbering 703-1013, so using it would corrupt every position by ~700.
            span = abs(uniprot_end - uniprot_start)
            pdb_start = _author_number(start)
            pdb_end = _author_number(end)
            if pdb_start is None and pdb_end is None:
                continue
            if pdb_start is None:
                pdb_start = pdb_end - span
            if pdb_end is None:
                pdb_end = pdb_start + span

            pdb_step = 1 if pdb_end >= pdb_start else -1
            uniprot_step = 1 if uniprot_end >= uniprot_start else -1
            length = min(abs(pdb_end - pdb_start), abs(uniprot_end - uniprot_start)) + 1
            for offset in range(length):
                mapping.setdefault(
                    pdb_start + offset * pdb_step, uniprot_start + offset * uniprot_step
                )
    return mapping


def parse_sifts_map_from_cif(
    text: str, accession: str = "", chain: str | None = None
) -> dict[int, int]:
    """Read the SIFTS correspondence out of a PDBe "updated" mmCIF atom_site loop.

    The columns wanted are ``_atom_site.pdbx_sifts_xref_db_num`` (the UniProt position)
    and ``_atom_site.auth_seq_id`` (the author number). A plain RCSB mmCIF has the atom
    site loop but not the SIFTS columns, so an atom_site loop lacking them is skipped
    rather than treated as the answer -- an entry can have more than one such loop.
    """
    wanted_accession = (accession or "").strip().upper().split("-")[0]
    wanted_chain = (chain or "").strip().upper()[:1] if chain else None

    mapping: dict[int, int] = {}
    lines = (text or "").splitlines()
    index, total = 0, len(lines)

    while index < total:
        if lines[index].strip() != "loop_":
            index += 1
            continue
        index += 1

        headers: list[str] = []
        while index < total and lines[index].strip().startswith("_"):
            headers.append(lines[index].strip())
            index += 1
        if not any(header.startswith("_atom_site.") for header in headers):
            continue

        def column(*suffixes: str) -> int | None:
            for suffix in suffixes:
                for position, header in enumerate(headers):
                    if header.lower().endswith(suffix.lower()):
                        return position
            return None

        uniprot_column = column(
            "pdbx_sifts_xref_db_num",
            "pdbx_sifts_xref_db_res_num",
            "pdbx_sifts_xref_db_residue_number",
        )
        accession_column = column("pdbx_sifts_xref_db_acc", "pdbx_sifts_xref_db_accession")
        author_column = column("auth_seq_id") or column("label_seq_id")
        chain_column = column("auth_asym_id")
        if chain_column is None:
            chain_column = column("label_asym_id")

        if uniprot_column is None or author_column is None:
            # An atom_site loop without SIFTS columns: skip its rows and keep looking.
            while index < total and not lines[index].strip().startswith(("loop_", "_", "#")):
                index += 1
            continue

        while index < total:
            row = lines[index].strip()
            if not row or row.startswith("#") or row == "loop_" or row.startswith("_"):
                break
            index += 1
            fields = row.split()
            if len(fields) < len(headers):
                continue
            if accession_column is not None and wanted_accession:
                found = fields[accession_column].upper().split("-")[0]
                if found not in {wanted_accession, "?", ".", ""}:
                    continue
            if wanted_chain and chain_column is not None:
                if fields[chain_column].strip().upper()[:1] != wanted_chain:
                    continue
            author_raw, uniprot_raw = fields[author_column], fields[uniprot_column]
            if author_raw in {"?", "."} or uniprot_raw in {"?", "."}:
                continue
            try:
                mapping.setdefault(int(float(author_raw)), int(float(uniprot_raw)))
            except (TypeError, ValueError):
                continue
    return mapping


def offset_mapping(
    uniprot_range: tuple[int, int] | None, pdb_range: tuple[int, int] | None
) -> dict[int, int]:
    """Line a chain's observed residue numbers up against UniProt's stated range.

    The fallback when SIFTS is unavailable. It assumes both run without gaps, which is
    why it ranks below SIFTS and why the UI says so.
    """
    if not uniprot_range or not pdb_range:
        return {}
    try:
        uniprot_start, uniprot_end = int(uniprot_range[0]), int(uniprot_range[1])
        pdb_start, pdb_end = int(pdb_range[0]), int(pdb_range[1])
    except (TypeError, ValueError):
        return {}
    pdb_step = 1 if pdb_end >= pdb_start else -1
    uniprot_step = 1 if uniprot_end >= uniprot_start else -1
    length = min(abs(pdb_end - pdb_start), abs(uniprot_end - uniprot_start)) + 1
    return {
        pdb_start + offset * pdb_step: uniprot_start + offset * uniprot_step
        for offset in range(length)
    }


def remap_site_rows(
    rows: list[dict[str, object]], mapping: PositionMapping
) -> list[dict[str, object]]:
    """Rewrite each row's ``position`` into author numbering for the viewer.

    A UniProt position present in two chains yields one row per chain. A position SIFTS
    does not place in this structure is dropped rather than drawn at its UniProt number,
    because a mark in the wrong place is worse than a mark missing: the caller reports the
    count so the omission is visible.
    """
    if mapping.is_identity:
        return [dict(row) for row in rows]

    remapped: list[dict[str, object]] = []
    for row in rows:
        try:
            uniprot_position = int(row.get("position"))
        except (TypeError, ValueError):
            continue
        for author_position in mapping.to_pdb(uniprot_position):
            new_row = dict(row)
            new_row["position"] = author_position
            new_row["uniprot_position"] = uniprot_position
            remapped.append(new_row)
    return remapped


def remap_positions(positions: list[int], mapping: PositionMapping) -> list[int]:
    """Author numbers for a list of UniProt positions, dropping any not in the structure."""
    if mapping.is_identity:
        return [int(position) for position in positions]
    out: list[int] = []
    for position in positions:
        out.extend(mapping.to_pdb(position))
    return sorted(set(out))


# ---------------------------------------------------------------------------
# Bio2Byte colour scale
# ---------------------------------------------------------------------------

#: Interpretation bands per predictor. The boundaries match
#: common.mutation_effect.MutationEffectInference, which labels the same predictions in
#: its own inference step; a unit test pins the two together so they cannot drift.
B2B_THRESHOLDS: dict[str, list[dict[str, object]]] = {
    "backbone": [
        {"min": None, "max": 0.69, "label": "flexible"},
        {"min": 0.69, "max": 0.80, "label": "context dependent"},
        {"min": 0.80, "max": 1.00, "label": "rigid"},
        {"min": 1.00, "max": None, "label": "membrane spanning"},
    ],
    "disoMine": [
        {"min": None, "max": 0.5, "label": "ordered"},
        {"min": 0.5, "max": None, "label": "disordered"},
    ],
    "earlyFolding": [
        {"min": None, "max": 0.169, "label": "not early folding"},
        {"min": 0.169, "max": None, "label": "early folding"},
    ],
}


def pseudocolor(minimum: float, maximum: float, value: float) -> str:
    """Green (low) to red (high) via an HSV sweep. Returns ``#rrggbb``."""
    minimum, maximum = float(minimum), float(maximum)
    if maximum == minimum:
        hue = 120.0
    else:
        hue = ((maximum - float(value)) / (maximum - minimum)) * 120.0
    red, green, blue = colorsys.hsv_to_rgb(hue / 360.0, 1.0, 1.0)
    return "#%02x%02x%02x" % tuple(int(255 * channel) for channel in (red, green, blue))


def b2b_interpretation(metric: str, value: object) -> str | None:
    """The band label for a value, e.g. "rigid", or None for an unbanded metric."""
    bands = B2B_THRESHOLDS.get(metric)
    if not bands or value is None:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    for band in bands:
        low, high = band["min"], band["max"]
        # Half-open the other way round -- (low, high] -- so a value sitting exactly on a
        # boundary lands in the lower band. That matches
        # common.mutation_effect.MutationEffectInference, which tests `value > threshold`,
        # and a test pins the two sets of boundaries together.
        if (low is None or numeric > low) and (high is None or numeric <= high):
            return str(band["label"])
    return None


def numeric_b2b_columns(dataframe: pd.DataFrame) -> list[str]:
    """Predictor columns worth colouring by: numeric, and not bookkeeping."""
    if dataframe is None or dataframe.empty:
        return []
    columns = []
    for column in dataframe.columns:
        name = str(column)
        if name in {"Position"} or "runtime" in name.lower() or "execution_time" in name.lower():
            continue
        if pd.api.types.is_numeric_dtype(dataframe[column]):
            columns.append(name)
    return columns


def b2b_value_range(dataframe: pd.DataFrame, metric: str) -> tuple[float, float] | None:
    if dataframe is None or dataframe.empty or metric not in dataframe.columns:
        return None
    values = pd.to_numeric(dataframe[metric], errors="coerce").dropna()
    if values.empty:
        return None
    return float(values.min()), float(values.max())


def b2b_legend_html(metric: str, minimum: float, maximum: float, steps: int = 48) -> str:
    """A colour bar for one predictor, with its interpretation bands marked.

    The band boundaries are drawn on the bar because a raw 0-1 gradient says nothing
    about where "ordered" stops and "disordered" starts.
    """
    try:
        low, high = float(minimum), float(maximum)
    except (TypeError, ValueError):
        return ""
    if high <= low:
        high = low + 1e-6

    stops = [
        f"{pseudocolor(low, high, low + (high - low) * (step / steps))} {100.0 * step / steps:.1f}%"
        for step in range(steps + 1)
    ]
    gradient = "linear-gradient(to right, " + ", ".join(stops) + ")"

    def percent(value: float) -> float:
        return max(0.0, min(100.0, (value - low) / (high - low) * 100.0))

    markers = ""
    top = (
        f"<span style='position:absolute;left:0;top:0;'>{low:.3g}</span>"
        f"<span style='position:absolute;right:0;top:0;'>{high:.3g}</span>"
    )
    labels = ""
    bands = B2B_THRESHOLDS.get(metric) or []
    for boundary in sorted({band["max"] for band in bands if band["max"] is not None}):
        if low <= float(boundary) <= high:
            at = percent(float(boundary))
            markers += (
                f"<div style='position:absolute;left:{at:.2f}%;top:0;height:14px;width:0;"
                f"border-left:2px solid #222;'></div>"
            )
            top += (
                f"<span style='position:absolute;left:{at:.2f}%;top:0;"
                f"transform:translateX(-50%);font-size:10px;color:#222;'>{boundary:g}</span>"
            )
    for band in bands:
        start = low if band["min"] is None else float(band["min"])
        stop = high if band["max"] is None else float(band["max"])
        visible_start, visible_stop = max(start, low), min(stop, high)
        if visible_stop - visible_start > (high - low) * 0.04:
            centre = percent((visible_start + visible_stop) / 2.0)
            labels += (
                f"<span style='position:absolute;left:{centre:.2f}%;"
                f"transform:translateX(-50%);white-space:nowrap;'>{band['label']}</span>"
            )

    return (
        "<div style='margin:6px 0 24px 0;font-size:12px;max-width:600px;'>"
        f"<div style='font-weight:600;margin-bottom:4px;'>{html.escape(metric)} "
        "&mdash; colour scale</div>"
        f"<div style='position:relative;height:16px;color:#333;margin-bottom:4px;'>{top}</div>"
        "<div style='position:relative;height:16px;'>"
        f"<div style='height:14px;border:1px solid #ccc;border-radius:3px;background:{gradient};'></div>"
        f"{markers}</div>"
        f"<div style='position:relative;height:16px;margin-top:3px;color:#333;'>{labels}</div>"
        "</div>"
    )


class ChainRangeSelect(Select):
    def __init__(self, chain_id: str, start: int | None, end: int | None) -> None:
        self.chain_id = chain_id
        self.start = start
        self.end = end

    def accept_residue(self, residue):  # noqa: ANN001
        if residue.parent.id != self.chain_id:
            return False
        seq_number = residue.id[1]
        if self.start is not None and seq_number < self.start:
            return False
        if self.end is not None and seq_number > self.end:
            return False
        return True


class StructureOps:
    @staticmethod
    def validate_pdb_id(pdb_id: str) -> str:
        pdb_key = pdb_id.strip().upper()
        if len(pdb_key) != 4 or not pdb_key.isalnum():
            raise ValueError(
                f"Invalid PDB ID '{pdb_id}'. Expected a 4-character RCSB identifier such as 2IVT."
            )
        return pdb_key

    @staticmethod
    def bfactor_pdb(pdb_path: Path, dataframe: pd.DataFrame, value_col: str, out_path: Path, chain: str | None = None) -> Path:
        values = dataframe[value_col].tolist()
        position_to_value = {index + 1: float(value) for index, value in enumerate(values) if value is not None and value == value}

        def set_bfactor(line: str, bfactor: float) -> str:
            return line[:60] + f"{bfactor:6.2f}" + line[66:]

        with pdb_path.open("r", encoding="utf-8", errors="ignore") as source, out_path.open("w", encoding="utf-8") as target:
            for line in source:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    target.write(line)
                    continue
                residue_chain = line[21].strip()
                if chain and residue_chain and residue_chain != chain:
                    target.write(line)
                    continue
                try:
                    residue_pos = int(line[22:26].strip())
                except Exception:
                    target.write(line)
                    continue
                target.write(set_bfactor(line, position_to_value[residue_pos]) if residue_pos in position_to_value else line)
        return out_path

    @staticmethod
    def chain_range_from_pdb(pdb_path: Path, chain_id: str) -> tuple[int, int]:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("chain", str(pdb_path))
        model = next(structure.get_models())
        chain = model[chain_id]
        positions = [res.id[1] for res in chain.get_residues() if res.id[0] == " "]
        return min(positions), max(positions)

    @staticmethod
    def chain_ranges_from_pdb(pdb_path: Path) -> dict[str, tuple[int, int]]:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("chains", str(pdb_path))
        model = next(structure.get_models())
        chain_ranges: dict[str, tuple[int, int]] = {}
        for chain in model.get_chains():
            positions = [res.id[1] for res in chain.get_residues() if res.id[0] == " "]
            if positions:
                chain_ranges[chain.id] = (min(positions), max(positions))
        return chain_ranges

    @staticmethod
    def save_chain_segment(source_pdb: Path, target_pdb: Path, chain_id: str, start: int | None, end: int | None) -> Path:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("segment", str(source_pdb))
        model = next(structure.get_models())
        if chain_id not in model:
            available = ", ".join(chain.id for chain in model.get_chains()) or "(none)"
            raise ValueError(
                f"Chain '{chain_id}' not found in {source_pdb.name}. Available chains: {available}."
            )

        chain = model[chain_id]
        positions = [res.id[1] for res in chain.get_residues() if res.id[0] == " "]
        if not positions:
            raise ValueError(f"Chain '{chain_id}' in {source_pdb.name} does not contain standard residues.")

        chain_start = min(positions)
        chain_end = max(positions)
        if start is not None and end is not None and start > end:
            raise ValueError(f"Invalid range for chain '{chain_id}': start {start} is greater than end {end}.")
        if start is not None and start > chain_end:
            raise ValueError(f"Start {start} is outside chain '{chain_id}' range {chain_start}-{chain_end}.")
        if end is not None and end < chain_start:
            raise ValueError(f"End {end} is outside chain '{chain_id}' range {chain_start}-{chain_end}.")

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(target_pdb), select=ChainRangeSelect(chain_id, start, end))
        segment_positions = [
            int(line[22:26].strip())
            for line in target_pdb.read_text(encoding="utf-8", errors="ignore").splitlines()
            if line.startswith(("ATOM", "HETATM")) and line[21].strip() == chain_id and line[22:26].strip()
        ]
        if not segment_positions:
            raise ValueError(
                f"Chain/range selection produced an empty segment for chain '{chain_id}' "
                f"with range {start or chain_start}-{end or chain_end} in {source_pdb.name}."
            )
        return target_pdb

    @staticmethod
    def run_tmalign(pdb1: Path, pdb2: Path, out_dir: Path, out_name: str = "aligned") -> tuple[Path, str]:
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = ["TM-align", str(pdb1.resolve()), str(pdb2.resolve()), "-o", out_name]
        try:
            process = subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as error:
            stderr = (error.stderr or "").strip()
            stdout = (error.stdout or "").strip()
            details = stderr or stdout or f"exit code {error.returncode}"
            raise RuntimeError(
                f"TM-align failed for {pdb1.name} vs {pdb2.name}: {details}"
            ) from error
        candidates = [
            out_dir / out_name,
            out_dir / f"{out_name}.pdb",
            out_dir / "TM_sup.pdb",
            out_dir / f"{out_name}.sup",
            out_dir / f"{out_name}.sup_all",
            out_dir / f"{out_name}.sup_all_atm",
        ]
        aligned_path = next((candidate for candidate in candidates if candidate.exists()), None)
        if aligned_path is None:
            files = ", ".join(sorted(path.name for path in out_dir.iterdir()))
            raise RuntimeError(f"No TM-align output found in {out_dir}. Files: {files}")
        return aligned_path, process.stdout

    @staticmethod
    def build_rin_graph(pdb_path: Path, chain: str = "A", cutoff: float = 8.0, atom_name: str = "CA") -> nx.Graph:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        model = next(structure.get_models())
        chain_obj = model[chain] if chain in model else next(model.get_chains())

        residues: list[int] = []
        coordinates: list[list[float]] = []
        for residue in chain_obj.get_residues():
            if residue.id[0] != " ":
                continue
            residue_number = residue.id[1]
            atom = residue[atom_name] if atom_name in residue else residue["CA"] if "CA" in residue else None
            if atom is None:
                continue
            residues.append(residue_number)
            coordinates.append(list(atom.get_coord()))

        graph = nx.Graph()
        for residue in residues:
            graph.add_node(residue)

        if not coordinates:
            return graph

        kdtree = KDTree(coordinates)
        pairs = kdtree.query_pairs(cutoff)
        for left, right in pairs:
            residue_left = residues[left]
            residue_right = residues[right]
            distance = float(((kdtree.data[left] - kdtree.data[right]) ** 2).sum() ** 0.5)
            graph.add_edge(residue_left, residue_right, distance=distance)

        graph.graph["chain"] = chain_obj.id
        return graph

    @staticmethod
    def rin_to_pyvis_html(
        graph: nx.Graph,
        out_html: Path,
        ptm_positions: list[int] | None = None,
        mutation_positions: list[int] | None = None,
        b2b_frame: pd.DataFrame | None = None,
        b2b_metric: str | None = None,
    ) -> Path:
        """Render the network.

        By default a node's colour encodes its site status: PTM, disease variant, both,
        or neither. Pass ``b2b_frame`` and ``b2b_metric`` to switch the *fill* to a
        Bio2Byte value on the green-to-red scale while the *border* keeps carrying the
        site status -- so a residue's predicted flexibility and its modification state
        are legible at the same time, which is the whole point of the overlay.

        Values are matched on the frame's ``Position`` column rather than by row order:
        a PDB chain does not necessarily start at residue 1 or run without gaps, so a
        positional lookup would silently shift every value.
        """
        ptm_set = set(ptm_positions or [])
        mutation_set = set(mutation_positions or [])

        values_by_position: dict[int, float] = {}
        value_range: tuple[float, float] | None = None
        if (
            b2b_frame is not None
            and not b2b_frame.empty
            and b2b_metric
            and b2b_metric in b2b_frame.columns
        ):
            numeric = pd.to_numeric(b2b_frame[b2b_metric], errors="coerce")
            if "Position" in b2b_frame.columns:
                positions = pd.to_numeric(b2b_frame["Position"], errors="coerce")
            else:
                positions = pd.Series(range(1, len(b2b_frame) + 1))
            for position, value in zip(positions, numeric):
                if pd.notna(position) and pd.notna(value):
                    values_by_position[int(position)] = float(value)
            if values_by_position:
                span = min(values_by_position.values()), max(values_by_position.values())
                value_range = span
        network = pyvis.network.Network(height="650px", width="100%", bgcolor="#ffffff", font_color="#222")
        # Use CDN resources so the generated HTML does not depend on local /lib/bindings assets.
        network.cdn_resources = "remote"
        network.set_options(
            """
            {
              "nodes": {
                "borderWidth": 1
              },
              "interaction": {
                "hover": true,
                "selectConnectedEdges": true,
                "multiselect": false
              },
              "physics": {
                "stabilization": true
              }
            }
            """
        )
        network.from_nx(graph)
        for node in network.nodes:
            residue = int(node["id"])
            if residue in ptm_set and residue in mutation_set:
                status_color, size, status = "#9467bd", 40, "PTM + variant"
            elif residue in ptm_set:
                status_color, size, status = "#1f77b4", 36, "PTM"
            elif residue in mutation_set:
                status_color, size, status = "#d62728", 36, "Variant"
            else:
                status_color, size, status = "#b0b0b0", 30, "Other"
            node["size"] = size

            if value_range is None:
                background = border = status_color
                node["title"] = f"Residue {residue} | {status}"
            else:
                # Fill by value, border by status. A residue the predictor has no value
                # for is grey rather than being coloured as if it sat at the low end.
                value = values_by_position.get(residue)
                if value is None:
                    background = "#d9d9d9"
                    readout = "n/a"
                else:
                    background = pseudocolor(value_range[0], value_range[1], value)
                    readout = f"{value:.3f}"
                    band = b2b_interpretation(str(b2b_metric), value)
                    if band:
                        readout = f"{readout} ({band})"
                border = "#000000" if status == "Other" else status_color
                node["borderWidth"] = 1 if status == "Other" else 3
                node["title"] = f"Residue {residue} | {status} | {b2b_metric}: {readout}"

            node["color"] = {
                "background": background,
                "border": border,
                "highlight": {"background": background, "border": border},
                "hover": {"background": background, "border": border},
            }
        out_html.parent.mkdir(parents=True, exist_ok=True)
        network.write_html(str(out_html), open_browser=False, notebook=False)

        html_text = out_html.read_text(encoding="utf-8", errors="ignore")
        click_widget = """
<div id="rin-click-info" style="margin:8px 0 10px 0;padding:8px 10px;background:#f7f7f7;border:1px solid #ddd;border-radius:6px;font-family:monospace;font-size:12px;">
Click a node to inspect residue details.
</div>
<script type="text/javascript">
(() => {
  const info = document.getElementById("rin-click-info");
  const defaultMsg = "Click a node to inspect residue details.";
  const bindClick = () => {
    if (typeof network === "undefined" || !network) {
      window.setTimeout(bindClick, 80);
      return;
    }
    network.on("click", (params) => {
      if (!params.nodes || params.nodes.length === 0) {
        info.textContent = defaultMsg;
        return;
      }
      const node = params.nodes[0];
      let degreeText = "";
      try {
        if (typeof edges !== "undefined" && edges) {
          const degree = edges.get({
            filter: (edge) => edge.from === node || edge.to === node
          }).length;
          degreeText = ` | degree=${degree}`;
        }
      } catch (_err) {}
      info.textContent = `Selected residue: ${node}${degreeText}`;
    });
  };
  bindClick();
})();
</script>
"""
        if "</body>" in html_text:
            html_text = html_text.replace("</body>", click_widget + "\n</body>", 1)
            out_html.write_text(html_text, encoding="utf-8")

        return out_html


class StructureViewerBuilder:
    @staticmethod
    def ptm_html(pdb_text: str, accession: str, ptm_rows: list[dict[str, object]], chain: str | None = None) -> str:
        # SER/THR/TYR are the phospho-acceptors Scop3P reports. The others arrive with
        # the UniProt PTM source -- ASN for N-linked glycosylation, LYS for acetylation
        # and methylation, CYS for lipidation -- and without an entry here a third of
        # the marks on a protein like P07949 would be an undifferentiated "other".
        colors = {
            "SER": "#d62728",
            "THR": "#2ca02c",
            "TYR": "#ff7f0e",
            "ASN": "#1f77b4",
            "LYS": "#9467bd",
            "CYS": "#8c564b",
        }
        payload = {"ptms": ptm_rows, "chain": chain or ""}
        safe_accession = html.escape(accession)
        uid = uuid.uuid4().hex
        panel_id = f"panel_{uid}"
        viewport_id = f"viewport_{uid}"
        return f"""<div style=\"position:relative;width:100%;height:700px;\"> 
  <div id=\"{panel_id}\" style=\"position:absolute;top:10px;left:10px;z-index:10;background:rgba(255,255,255,.9);padding:8px;border-radius:8px;\"> 
    <b>{safe_accession}</b><br/>PTM viewer (PDB/AlphaFold)
  </div>
  <div id=\"{viewport_id}\" style=\"width:100%;height:100%;\"></div>
</div>
<script>
(() => {{
  const pdbText = {json.dumps(pdb_text)};
  const payload = {json.dumps(payload)};
  const colorMap = {json.dumps(colors)};
  const stageEl = document.getElementById('{viewport_id}');
  const panelEl = document.getElementById('{panel_id}');
  if (!stageEl || !panelEl) return;

  const paint = () => {{
    const stage = new window.NGL.Stage(stageEl, {{ backgroundColor: 'white' }});
    const blob = new Blob([pdbText], {{type:'text/plain'}});
    stage.loadFile(blob, {{ext:'pdb'}}).then(comp => {{
      comp.addRepresentation('cartoon', {{color:'lightgrey'}});
      for (const row of payload.ptms) {{
        const pos = row.position;
        const residue = String(row.residue || '').trim();
        if (!pos) continue;
        const color = colorMap[residue] || '#7b241c';
        const chainPart = payload.chain ? ` AND :${{payload.chain}}` : '';
        const sele = `${{pos}}${{chainPart}}`;
        comp.addRepresentation('ball+stick', {{sele: sele, color: color}});
      }}
      comp.autoView();
    }});
  }};

  if (window.NGL) {{ paint(); return; }}
  const s = document.createElement('script');
  s.src = 'https://unpkg.com/ngl@latest/dist/ngl.js';
  s.onload = paint;
  s.onerror = () => {{ panelEl.innerHTML += '<div style="color:#b00020">Could not load NGL assets.</div>'; }};
  document.head.appendChild(s);
}})();
</script>
"""

    @staticmethod
    def b2b_html(pdb_text: str, accession: str, metric: str) -> str:
        safe_accession = html.escape(accession)
        uid = uuid.uuid4().hex
        panel_id = f"panel_{uid}"
        viewport_id = f"viewport_{uid}"
        return f"""<div style=\"position:relative;width:100%;height:700px;\"> 
  <div id=\"{panel_id}\" style=\"position:absolute;top:10px;left:10px;z-index:10;background:rgba(255,255,255,.9);padding:8px;border-radius:8px;\"> 
    <b>{safe_accession}</b><br/>Bio2Byte metric: {html.escape(metric)}
  </div>
  <div id=\"{viewport_id}\" style=\"width:100%;height:100%;\"></div>
</div>
<script>
(() => {{
  const pdbText = {json.dumps(pdb_text)};
  const stageEl = document.getElementById('{viewport_id}');
  const panelEl = document.getElementById('{panel_id}');
  if (!stageEl || !panelEl) return;
  const paint = () => {{
    const stage = new window.NGL.Stage(stageEl, {{ backgroundColor: 'white' }});
    const blob = new Blob([pdbText], {{type:'text/plain'}});
    stage.loadFile(blob, {{ext:'pdb'}}).then(comp => {{
      comp.addRepresentation('cartoon', {{colorScheme: 'bfactor'}});
      comp.autoView();
    }});
  }};
  if (window.NGL) {{ paint(); return; }}
  const s = document.createElement('script');
  s.src = 'https://unpkg.com/ngl@latest/dist/ngl.js';
  s.onload = paint;
  s.onerror = () => {{ panelEl.innerHTML += '<div style="color:#b00020">Could not load NGL assets.</div>'; }};
  document.head.appendChild(s);
}})();
</script>
"""
