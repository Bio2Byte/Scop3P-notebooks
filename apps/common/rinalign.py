"""RINAlign: residue interaction network construction, diffing and alignment.

Ported from ``notebooks/RINAlign_align_and compare_networks.ipynb`` cells 2-5. The
science is unchanged; the differences from the notebook are deliberate and listed in
``docs/use-cases/rinalign.md``. In short:

* every HTTP call has a timeout, and structure downloads are cached in one
  per-session working directory instead of leaking a fresh ``mkdtemp`` per click;
* the three structure-discovery fallbacks log which one failed instead of hiding a
  network outage behind "no structures found";
* contact detection uses a KD-tree rather than an O(n^2) Python loop, which is the
  same technique ``common.structure_viz.build_rin_graph`` already uses. The result is
  identical: ``query_pairs(r)`` returns pairs at ``distance <= r``, matching the
  notebook's ``if d <= cutoff``.

The BioPython import guard from the notebook is dropped: biopython is a hard
dependency here (``requirements-shiny.txt``), and the notebook's fallback printed to
stdout, which would land in the structured log.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import networkx as nx
import numpy as np
import requests
from Bio.PDB import MMCIFParser, PDBParser
from scipy.optimize import linear_sum_assignment
from scipy.spatial import KDTree

from common.cache import memoize
from common.http_lookup import lookup as _http_lookup
from common.logging_utils import get_logger
from common.structure_labels import ALPHAFOLD_OPTION_LABEL, chain_option_label
from common.services import Scop3PClient

LOGGER = get_logger("scop3p.common.rinalign")


def http_lookup(url: str, **kwargs):
    """An annotation lookup under the policy shared with every other protocol."""
    return _http_lookup(url, logger=LOGGER, **kwargs)

UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb"
PDBE_BASE_URL = "https://www.ebi.ac.uk/pdbe/api/mappings"
PROTEINS_BASE_URL = "https://www.ebi.ac.uk/proteins/api"
ALPHAFOLD_BASE_URL = "https://alphafold.ebi.ac.uk/api/prediction"
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download"

# Modified residues are whitelisted so phospho-residues become real network nodes
# rather than being skipped as heteroatoms.
THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E",
    "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "PTR": "Y", "SEP": "S", "TPO": "T",
}

# Hydrophobic / aromatic / polar / positive / negative / small, used by restype_sim.
RESIDUE_GROUPS: Dict[str, str] = {
    "ALA": "H", "VAL": "H", "LEU": "H", "ILE": "H", "MET": "H",
    "PHE": "A", "TRP": "A", "TYR": "A", "PTR": "A",
    "SER": "P", "SEP": "P", "THR": "P", "TPO": "P", "ASN": "P", "GLN": "P",
    "ARG": "+", "LYS": "+", "HIS": "+",
    "ASP": "-", "GLU": "-",
    "GLY": "S", "CYS": "S", "PRO": "S",
}

# align_rins allocates a dense n1 x n2 float matrix and fills it from a Python double
# loop. At 2000 nodes that is 4M iterations and 32 MB, on the ASGI event loop, for
# every connected session. Refuse rather than freeze the worker.
MAX_ALIGNMENT_NODES = 1500


@dataclass(slots=True, frozen=True)
class StructureEntry:
    """One selectable structure: an AlphaFold model or a PDB entry chain."""

    pdb_id: str
    chain_id: str
    label: str
    source: str  # "AlphaFold" | "PDB"
    unp_start: Any = "?"
    unp_end: Any = "?"
    resolution: Optional[float] = None
    method: str = ""
    pdb_url: str = ""
    cif_url: str = ""

    @property
    def key(self) -> str:
        return f"{self.pdb_id}_{self.chain_id}"

    @property
    def display(self) -> str:
        """The short label the view builders put into plot titles."""
        return f"{self.pdb_id} {self.chain_id}".strip()


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------


class RINAlignService:
    """Network-facing half of RINAlign.

    ``workdir`` is one directory per session. The notebook called
    ``tempfile.mkdtemp()`` inside ``download_structure`` and always wrote to the same
    basename, so every Generate click leaked a directory; here downloads are named
    after the entry and reused.
    """

    def __init__(self, workdir: Path, timeout: int = 30) -> None:
        # Scop3P access goes through the shared client so the endpoint and its field
        # names are defined in exactly one place. This module used to call the pre-v1
        # URL itself, which the API now answers with the single-page-app HTML: the
        # broad except below turned that into "0 PTMs" for every protein rather than
        # an error, so the overlay silently did nothing.
        self.scop3p_client = Scop3PClient(timeout=timeout)
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout

    # -- protein metadata ---------------------------------------------------

    @memoize(name="uniprot.entry.info")
    def fetch_uniprot_info(self, accession: str) -> Dict[str, Any]:
        response = http_lookup(
            f"{UNIPROT_BASE_URL}/{accession}.json",
            headers={"Accept": "application/json"},
        )
        if response.status_code != 200:
            raise ValueError(f"UniProt ID not found ({response.status_code})")
        payload = response.json()

        info: Dict[str, Any] = {
            "accession": payload.get("primaryAccession", accession),
            "protein_name": "",
            "gene_name": "",
            "organism": "",
            "length": 0,
            "_raw": payload,
        }

        description = payload.get("proteinDescription", {})
        recommended = description.get("recommendedName", description.get("submittedName", [{}]))
        if isinstance(recommended, list):
            recommended = recommended[0] if recommended else {}
        full_name = recommended.get("fullName", {})
        info["protein_name"] = (
            full_name.get("value", "") if isinstance(full_name, dict) else str(full_name)
        )

        genes = payload.get("genes", [{}])
        if genes:
            gene_name = genes[0].get("geneName", {})
            info["gene_name"] = (
                gene_name.get("value", "") if isinstance(gene_name, dict) else str(gene_name)
            )

        info["organism"] = payload.get("organism", {}).get("scientificName", "")
        info["length"] = payload.get("sequence", {}).get("length", 0)
        return info

    # -- structure discovery ------------------------------------------------

    def fetch_pdb_structures(
        self, accession: str, uniprot_data: Optional[Dict[str, Any]] = None
    ) -> List[StructureEntry]:
        """Discover PDB entry chains covering this accession.

        Three strategies in precedence order: UniProt cross-references (already
        fetched, so free and most reliable), the PDBe ``best_structures`` endpoint,
        then the PDBe SIFTS ``mappings/uniprot`` endpoint. Each logs its own failure
        so a network outage is distinguishable from a protein with no PDB entries --
        the notebook's bare ``except: pass`` made those two cases identical.
        """
        entries: List[StructureEntry] = []
        seen: set[str] = set()

        if uniprot_data:
            entries.extend(self._from_uniprot_xrefs(uniprot_data, seen))

        if not entries:
            entries.extend(self._from_best_structures(accession, seen))

        if not entries:
            entries.extend(self._from_sifts_mappings(accession, seen))

        entries.sort(key=lambda entry: entry.resolution if entry.resolution else 999)
        return entries

    def _from_uniprot_xrefs(
        self, uniprot_data: Dict[str, Any], seen: set[str]
    ) -> List[StructureEntry]:
        entries: List[StructureEntry] = []
        for xref in uniprot_data.get("uniProtKBCrossReferences", []) or []:
            if xref.get("database") != "PDB":
                continue
            pdb_id = (xref.get("id") or "").upper()
            if not pdb_id:
                continue

            properties = {p.get("key"): p.get("value") for p in xref.get("properties", []) or []}
            method = properties.get("Method", "") or ""
            resolution_text = properties.get("Resolution", "") or ""
            chains_text = properties.get("Chains", "") or ""

            resolution: Optional[float] = None
            if resolution_text:
                try:
                    resolution = float(resolution_text.replace(" A", ""))
                except ValueError:
                    resolution = None

            # "A/B=1-100" or "A=10-200, B=10-200"
            chain_entries = [c.strip() for c in chains_text.split(",")] if chains_text else ["?"]
            for chain_entry in chain_entries:
                parts = chain_entry.split("=")
                chain_ids = parts[0].strip() if parts else "?"
                residue_range = parts[1].strip() if len(parts) > 1 else ""

                start: Any = "?"
                end: Any = "?"
                if "-" in residue_range:
                    bounds = residue_range.split("-")
                    try:
                        start, end = int(bounds[0]), int(bounds[1])
                    except (ValueError, IndexError):
                        start, end = "?", "?"

                for chain in chain_ids.split("/"):
                    chain = chain.strip()
                    if not chain:
                        continue
                    key = f"{pdb_id}_{chain}"
                    if key in seen:
                        continue
                    seen.add(key)

                    label = chain_option_label(
                        pdb_id,
                        chain,
                        unp_start=start,
                        unp_end=end,
                        method=method,
                        resolution=resolution,
                    )

                    entries.append(
                        StructureEntry(
                            pdb_id=pdb_id,
                            chain_id=chain,
                            label=label,
                            source="PDB",
                            unp_start=start,
                            unp_end=end,
                            resolution=resolution,
                            method=method,
                        )
                    )
        return entries

    def _from_best_structures(self, accession: str, seen: set[str]) -> List[StructureEntry]:
        try:
            response = http_lookup(f"{PDBE_BASE_URL}/best_structures/{accession}")
            if response.status_code != 200:
                return []
            payload = response.json()
        except Exception as error:  # noqa: BLE001 - logged, then we fall through
            LOGGER.warning(
                "best_structures lookup failed accession=%s error=%s",
                accession,
                error,
                extra={"event": "fetch_pdb_structures"},
            )
            return []

        entries: List[StructureEntry] = []
        for key in (accession.lower(), accession.upper(), accession):
            if key not in payload:
                continue
            records = payload[key]
            if not isinstance(records, list):
                break
            for record in records:
                pdb_id = (record.get("pdb_id") or "").upper()
                chain = record.get("chain_id", "?")
                if f"{pdb_id}_{chain}" in seen:
                    continue
                seen.add(f"{pdb_id}_{chain}")
                resolution = record.get("resolution")
                start = record.get("unp_start", "?")
                end = record.get("unp_end", "?")
                label = chain_option_label(
                    pdb_id, chain, unp_start=start, unp_end=end, resolution=resolution
                )
                entries.append(
                    StructureEntry(
                        pdb_id=pdb_id,
                        chain_id=chain,
                        label=label,
                        source="PDB",
                        unp_start=start,
                        unp_end=end,
                        resolution=resolution,
                    )
                )
            break
        return entries

    def _from_sifts_mappings(self, accession: str, seen: set[str]) -> List[StructureEntry]:
        try:
            response = http_lookup(f"{PDBE_BASE_URL}/uniprot/{accession}")
            if response.status_code != 200:
                return []
            payload = response.json()
        except Exception as error:  # noqa: BLE001 - logged, then we give up
            LOGGER.warning(
                "SIFTS mappings lookup failed accession=%s error=%s",
                accession,
                error,
                extra={"event": "fetch_pdb_structures"},
            )
            return []

        entries: List[StructureEntry] = []
        for _accession, value in payload.items():
            pdb_block = value.get("PDB", value) if isinstance(value, dict) else {}
            if not isinstance(pdb_block, dict):
                continue
            for pdb_id, chain_list in pdb_block.items():
                if not isinstance(chain_list, list):
                    continue
                for chain_info in chain_list:
                    if not isinstance(chain_info, dict):
                        continue
                    chain = chain_info.get("chain_id", chain_info.get("struct_asym_id", "?"))
                    key = f"{pdb_id.upper()}_{chain}"
                    if key in seen:
                        continue
                    seen.add(key)
                    start = _boundary(chain_info, "unp_start", "start")
                    end = _boundary(chain_info, "unp_end", "end")
                    resolution = chain_info.get("resolution")
                    label = chain_option_label(
                        pdb_id, chain, unp_start=start, unp_end=end, resolution=resolution
                    )
                    entries.append(
                        StructureEntry(
                            pdb_id=pdb_id.upper(),
                            chain_id=chain,
                            label=label,
                            source="PDB",
                            unp_start=start,
                            unp_end=end,
                            resolution=resolution,
                        )
                    )
        return entries

    def check_alphafold(self, accession: str) -> Optional[StructureEntry]:
        try:
            response = http_lookup(f"{ALPHAFOLD_BASE_URL}/{accession}")
            if response.status_code != 200:
                return None
            payload = response.json()
            if not payload:
                return None
        except Exception as error:  # noqa: BLE001 - absence is a normal outcome
            LOGGER.warning(
                "AlphaFold lookup failed accession=%s error=%s",
                accession,
                error,
                extra={"event": "check_alphafold"},
            )
            return None

        record = payload[0]
        return StructureEntry(
            pdb_id=f"AF-{accession}",
            chain_id="A",
            label=ALPHAFOLD_OPTION_LABEL,
            source="AlphaFold",
            pdb_url=record.get("pdbUrl", "") or "",
            cif_url=record.get("cifUrl", "") or "",
        )

    # -- structure download -------------------------------------------------

    def download_structure(self, entry: StructureEntry) -> Tuple[Path, str]:
        """Fetch coordinates for one entry, cached in the session working directory."""
        if entry.source == "AlphaFold":
            url = entry.pdb_url or entry.cif_url
            if not url:
                raise ValueError(f"No AlphaFold coordinate URL for {entry.pdb_id}.")
            suffix = ".pdb" if "pdb" in url.lower() else ".cif"
            attempts: List[Tuple[str, str]] = [(url, suffix)]
        else:
            # RCSB serves mmCIF for every entry and PDB only for small ones, so try
            # mmCIF first and fall back, as the notebook did.
            lower = entry.pdb_id.lower()
            attempts = [
                (f"{RCSB_DOWNLOAD_URL}/{lower}.cif", ".cif"),
                (f"{RCSB_DOWNLOAD_URL}/{lower}.pdb", ".pdb"),
            ]

        last_error = ""
        for url, suffix in attempts:
            target = self.workdir / f"{entry.pdb_id}_{entry.chain_id}{suffix}"
            if target.exists() and target.stat().st_size > 0:
                return target, target.read_text(encoding="utf-8", errors="replace")
            try:
                response = requests.get(url, timeout=self.timeout)
            except Exception as error:  # noqa: BLE001
                last_error = f"{type(error).__name__}: {error}"
                continue
            if response.status_code != 200:
                last_error = f"HTTP {response.status_code}"
                continue
            target.write_text(response.text, encoding="utf-8")
            return target, response.text

        raise ValueError(
            f"Download failed for {entry.pdb_id} chain {entry.chain_id} ({last_error})."
        )

    # -- PTM and variant positions -----------------------------------------

    def fetch_scop3p_ptm_positions(self, accession: str) -> set[int]:
        """Scop3P PTM positions. Scop3P mainly covers human phosphoproteins.

        An accession Scop3P does not cover returns an empty set: the API answers 200
        with an empty list for both an uncovered protein and an unknown accession, and
        neither is an error. A transport or endpoint problem raises
        :class:`~common.services.Scop3PApiError` instead, so a broken service is
        reported rather than being indistinguishable from a protein with no sites.
        """
        dataframe = self.scop3p_client.fetch_modifications(accession)
        if dataframe.empty or "position" not in dataframe.columns:
            return set()
        return _int_set(dataframe["position"].dropna().tolist())

    def fetch_uniprot_ptm_positions(self, accession: str) -> set[int]:
        """Single-residue UniProt PTM features."""
        try:
            response = http_lookup(
                f"{PROTEINS_BASE_URL}/features/{accession}",
                headers={"Accept": "application/json"},
                params=[("categories", "PTM")],
            )
            response.raise_for_status()
            payload = response.json()
        except Exception as error:  # noqa: BLE001
            LOGGER.warning(
                "UniProt PTM lookup failed accession=%s error=%s",
                accession,
                error,
                extra={"event": "fetch_uniprot_ptm_positions"},
            )
            return set()

        positions = []
        for feature in payload.get("features", []) or []:
            if feature.get("category") != "PTM":
                continue
            begin, end = feature.get("begin"), feature.get("end")
            # Ranges are skipped: a multi-residue feature has no single site to mark.
            if begin is None or str(begin) != str(end):
                continue
            positions.append(begin)
        return _int_set(positions)

    def fetch_uniprot_variant_positions(
        self, accession: str, disease_only: bool = True
    ) -> set[int]:
        try:
            response = http_lookup(
                f"{PROTEINS_BASE_URL}/variation/{accession}",
                headers={"Accept": "application/json"},
            )
            response.raise_for_status()
            payload = response.json()
        except Exception as error:  # noqa: BLE001
            LOGGER.warning(
                "UniProt variant lookup failed accession=%s error=%s",
                accession,
                error,
                extra={"event": "fetch_uniprot_variant_positions"},
            )
            return set()

        positions = []
        for feature in payload.get("features", []) or []:
            if feature.get("type") != "VARIANT":
                continue
            if disease_only and not any(
                association.get("disease") is True
                for association in (feature.get("association") or [])
            ):
                continue
            positions.append(feature.get("begin"))
        return _int_set(positions)


def _boundary(chain_info: Dict[str, Any], flat_key: str, nested_key: str) -> Any:
    """SIFTS records carry either ``unp_start`` or ``start.residue_number``."""
    if flat_key in chain_info:
        return chain_info[flat_key]
    nested = chain_info.get(nested_key)
    if isinstance(nested, dict):
        return nested.get("residue_number", "?")
    return "?"


def _int_set(values: Iterable[Any]) -> set[int]:
    out: set[int] = set()
    for value in values:
        if value is None:
            continue
        try:
            out.add(int(value))
        except (TypeError, ValueError):
            continue
    return out


# ---------------------------------------------------------------------------
# RIN construction
# ---------------------------------------------------------------------------


def get_cb(residue) -> Optional[np.ndarray]:
    """CB coordinates, falling back to CA for glycine and incomplete residues."""
    name = residue.get_resname().strip()
    if name == "GLY":
        return residue["CA"].get_vector().get_array() if "CA" in residue else None
    if "CB" in residue:
        return residue["CB"].get_vector().get_array()
    if "CA" in residue:
        return residue["CA"].get_vector().get_array()
    return None


def build_rin(
    filepath: str | Path,
    chain_id: Optional[str] = None,
    cutoff: float = 8.0,
) -> Tuple[nx.Graph, Dict[str, Dict[str, Any]]]:
    """Build a residue interaction network from a structure file.

    Nodes are residues keyed ``RESNAME_seqid``; edges join residues whose CB atoms
    (CA for glycine) lie within ``cutoff`` angstroms, excluding sequence-adjacent
    pairs within the same chain.

    Contacts come from a KD-tree rather than the notebook's O(n^2) Python loop. The
    edge set is identical -- ``KDTree.query_pairs(r)`` returns pairs with
    ``distance <= r``, the same comparison the notebook made -- but a 2000-residue
    model no longer blocks the event loop for every connected session.
    """
    path = Path(filepath)
    parser = MMCIFParser(QUIET=True) if path.suffix == ".cif" else PDBParser(QUIET=True)
    structure = parser.get_structure("p", str(path))

    residues: List[Dict[str, Any]] = []
    for chain in structure[0]:
        if chain_id and chain.get_id() != chain_id:
            continue
        for residue in chain:
            hetero_flag, sequence_id, _insertion = residue.get_id()
            name = residue.get_resname().strip()
            if hetero_flag != " " and name not in THREE_TO_ONE:
                continue
            if name not in THREE_TO_ONE:
                continue
            coordinates = get_cb(residue)
            if coordinates is None:
                continue
            residues.append(
                {
                    "nid": f"{name}_{sequence_id}",
                    "pos": sequence_id,
                    "rn": name,
                    "aa": THREE_TO_ONE.get(name, "X"),
                    "ch": chain.get_id(),
                    "cb": coordinates,
                }
            )

    graph = nx.Graph()
    for residue in residues:
        graph.add_node(
            residue["nid"],
            position=residue["pos"],
            resname=residue["rn"],
            one_letter=residue["aa"],
            chain=residue["ch"],
        )

    if len(residues) > 1:
        coordinates = np.asarray([residue["cb"] for residue in residues], dtype=float)
        tree = KDTree(coordinates)
        for left, right in tree.query_pairs(cutoff):
            first, second = residues[left], residues[right]
            if first["ch"] == second["ch"] and abs(first["pos"] - second["pos"]) <= 1:
                continue
            distance = float(np.linalg.norm(first["cb"] - second["cb"]))
            graph.add_edge(first["nid"], second["nid"], distance=round(distance, 2))

    return graph, {residue["nid"]: residue for residue in residues}


def rin_html(graph: nx.Graph, label: str) -> str:
    """Small stat block shown beside each structure selector."""
    return (
        "<div style='padding:8px;background:#f8f9fa;border-radius:6px;margin:4px 0;'>"
        f"<b>{label}</b><br>Nodes: <b>{graph.number_of_nodes()}</b> | "
        f"Edges: <b>{graph.number_of_edges()}</b> | "
        f"Density: <b>{nx.density(graph):.3f}</b></div>"
    )


# ---------------------------------------------------------------------------
# Same-protein diff
# ---------------------------------------------------------------------------


def diff_rins(graph_left: nx.Graph, graph_right: nx.Graph) -> Dict[str, Any]:
    """Compare two networks of the same protein by matching residue positions.

    Only positions present in both structures are compared, so a contact that is
    "lost" is genuinely absent rather than merely unresolved in one model.
    """
    positions_left = {
        graph_left.nodes[node].get("position"): node
        for node in graph_left
        if graph_left.nodes[node].get("position") is not None
    }
    positions_right = {
        graph_right.nodes[node].get("position"): node
        for node in graph_right
        if graph_right.nodes[node].get("position") is not None
    }

    matched = sorted(set(positions_left) & set(positions_right))
    only_left = sorted(set(positions_left) - set(positions_right))
    only_right = sorted(set(positions_right) - set(positions_left))

    mutations = [
        {
            "position": position,
            "left": graph_left.nodes[positions_left[position]].get("resname", "?"),
            "right": graph_right.nodes[positions_right[position]].get("resname", "?"),
        }
        for position in matched
        if graph_left.nodes[positions_left[position]].get("resname", "?")
        != graph_right.nodes[positions_right[position]].get("resname", "?")
    ]

    def edges_by_position(graph: nx.Graph, position_map: Dict[Any, str]) -> set:
        node_to_position = {node: position for position, node in position_map.items()}
        return {
            (
                min(node_to_position[u], node_to_position[v]),
                max(node_to_position[u], node_to_position[v]),
            )
            for u, v in graph.edges()
            if u in node_to_position and v in node_to_position
        }

    edges_left = edges_by_position(graph_left, positions_left)
    edges_right = edges_by_position(graph_right, positions_right)

    matched_set = set(matched)
    comparable_left = {e for e in edges_left if e[0] in matched_set and e[1] in matched_set}
    comparable_right = {e for e in edges_right if e[0] in matched_set and e[1] in matched_set}

    conserved = sorted(comparable_left & comparable_right)
    lost = sorted(comparable_left - comparable_right)
    gained = sorted(comparable_right - comparable_left)
    total = len(comparable_left | comparable_right)

    only_left_set = set(only_left)
    only_a_edges = sorted(
        {(min(i, j), max(i, j)) for i, j in edges_left if i in only_left_set or j in only_left_set}
    )

    impact = []
    for position in matched:
        lost_here = len([edge for edge in lost if position in edge])
        gained_here = len([edge for edge in gained if position in edge])
        impact.append(
            {
                "position": position,
                "resname_L": graph_left.nodes[positions_left[position]].get("resname", "?"),
                "resname_R": graph_right.nodes[positions_right[position]].get("resname", "?"),
                "deg_L": graph_left.degree(positions_left[position]),
                "deg_R": graph_right.degree(positions_right[position]),
                "lost": lost_here,
                "gained": gained_here,
                "net_change": gained_here - lost_here,
                "is_mutation": any(m["position"] == position for m in mutations),
            }
        )
    impact.sort(key=lambda row: abs(row["net_change"]), reverse=True)

    return {
        "conserved": conserved,
        "lost": lost,
        "gained": gained,
        "onlyA_edges": only_a_edges,
        "jaccard": len(conserved) / total if total > 0 else 1.0,
        "mutations": mutations,
        "residue_impact": impact,
        "matched_pos": matched,
        "only_left_pos": only_left,
        "only_right_pos": only_right,
        "edges_L": len(comparable_left),
        "edges_R": len(comparable_right),
        "pos_to_left": positions_left,
        "pos_to_right": positions_right,
    }


# ---------------------------------------------------------------------------
# Cross-protein graph alignment
# ---------------------------------------------------------------------------


def restype_sim(first: str, second: str) -> float:
    """1.0 for identical residues, 0.5 within a physicochemical group, else 0.0."""
    if first == second:
        return 1.0
    group_first = RESIDUE_GROUPS.get(first, "?")
    group_second = RESIDUE_GROUPS.get(second, "?")
    return 0.5 if group_first == group_second and group_first != "?" else 0.0


def wl_sigs(graph: nx.Graph, k: int = 3) -> Dict[str, List[str]]:
    """Weisfeiler-Lehman signatures: one label per refinement round, plus the seed.

    Each node ends up with ``k + 1`` labels describing its neighbourhood at
    increasing radius, which is what lets two graphs be compared structurally rather
    than only by residue identity.
    """
    labels = {node: graph.nodes[node].get("resname", f"D{graph.degree(node)}") for node in graph}
    history: Dict[str, List[str]] = {node: [labels[node]] for node in graph}
    for _round in range(k):
        updated = {}
        for node in graph:
            neighbours = sorted(labels[other] for other in graph.neighbors(node))
            updated[node] = hashlib.md5(
                (labels[node] + "|" + ",".join(neighbours)).encode()
            ).hexdigest()[:8]
        labels = updated
        for node in graph:
            history[node].append(labels[node])
    return history


def align_rins(graph_one: nx.Graph, graph_two: nx.Graph, wl_depth: int = 3) -> Dict[str, Any]:
    """Align two networks of different proteins and score the shared topology.

    Node similarity is ``0.25 * residue type + 0.25 * degree + 0.5 * WL agreement``,
    solved to a one-to-one node mapping by Hungarian assignment, then scored by the
    Jaccard overlap of the mapped edge sets.
    """
    nodes_one, nodes_two = list(graph_one), list(graph_two)
    count_one, count_two = len(nodes_one), len(nodes_two)

    if max(count_one, count_two) > MAX_ALIGNMENT_NODES:
        raise ValueError(
            f"Alignment is limited to {MAX_ALIGNMENT_NODES} residues per network "
            f"(got {count_one} and {count_two}). Pick a smaller chain, or lower the "
            "contact cutoff and compare a domain."
        )

    degrees_one = [graph_one.degree(node) for node in nodes_one]
    degrees_two = [graph_two.degree(node) for node in nodes_two]
    max_degree = max(max(degrees_one, default=0), max(degrees_two, default=0))

    signatures_one = wl_sigs(graph_one, wl_depth)
    signatures_two = wl_sigs(graph_two, wl_depth)

    similarity = np.zeros((count_one, count_two))
    for i in range(count_one):
        for j in range(count_two):
            score = 0.25 * restype_sim(
                graph_one.nodes[nodes_one[i]].get("resname", "UNK"),
                graph_two.nodes[nodes_two[j]].get("resname", "UNK"),
            )
            score += 0.25 * (
                1.0 - abs(degrees_one[i] - degrees_two[j]) / max_degree
                if max_degree > 0
                else 1.0
            )
            score += (
                0.5
                * sum(
                    1
                    for a, b in zip(signatures_one[nodes_one[i]], signatures_two[nodes_two[j]])
                    if a == b
                )
                / (wl_depth + 1)
            )
            similarity[i, j] = score

    rows, columns = linear_sum_assignment(-similarity)
    mapping = [
        (nodes_one[row], nodes_two[column], similarity[row, column])
        for row, column in zip(rows, columns)
    ]

    one_to_two = {a: b for a, b, _ in mapping}
    mapped_edges = {
        tuple(sorted([one_to_two[u], one_to_two[v]]))
        for u, v in graph_one.edges()
        if u in one_to_two and v in one_to_two
    }
    edges_two = {tuple(sorted([u, v])) for u, v in graph_two.edges()}

    conserved = sorted(mapped_edges & edges_two)
    union = len(mapped_edges | edges_two)
    return {
        "mapping": mapping,
        "conserved": conserved,
        "only_G1": sorted(mapped_edges - edges_two),
        "only_G2": sorted(edges_two - mapped_edges),
        "jaccard": len(conserved) / union if union > 0 else 0.0,
    }
