"""Annotations: structures, numbering, and sites.

Three jobs, in the order they have to happen:

    1. find the PDB entries covering a UniProt accession
    2. build the UniProt <-> author numbering map for the chosen entry
    3. fetch PTMs and disease variants, all UniProt-numbered

Step 2 is not optional plumbing. Scop3P and UniProt report positions in UniProt
numbering; PDB entries number residues however the depositor chose. Skipping the
map puts every mark on the wrong residue, and does so silently, which is worse
than failing. AlphaFold models are the one case where the two coincide.

The fetchers follow the working implementations in the Scop3P notebook rather
than reinventing them. Each returns plain data and raises ValueError with a
message meant for display.
"""

from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

# No custom User-Agent. The official Scop3P client and the Scop3P notebook both
# use a bare requests session, and identifying as something unfamiliar can earn
# an HTML block page instead of JSON -- which surfaces as a JSON decode error at
# character 0 and looks like a parsing bug rather than a refused request.
USER_AGENT = ""
TIMEOUT = 60

UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}.json"
# The proteins API is what the Scop3P app uses for PTM features, and its
# response shape differs from the UniProtKB entry endpoint: positions arrive as
# flat begin/end strings rather than nested location objects.
EBI_FEATURES = "https://www.ebi.ac.uk/proteins/api/features/{accession}"
PDBE_BEST = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{accession}"
PDBE_UNIPROT = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
# Current Scop3P REST API. Path parameters, and a response that is a bare JSON
# array rather than an object wrapping a "modifications" key.
SCOP3P_V1_MODS = "https://iomics.ugent.be/scop3p/api/v1/proteins/{accession}/modifications"
# Older query-parameter form, kept as a fallback.
SCOP3P_MODS = "https://iomics.ugent.be/scop3p/api/modifications"
UNIPROT_VARIANTS = "https://www.ebi.ac.uk/proteins/api/variation/{accession}"

# Site categories, coloured to match the wider workflow. A residue that is both
# modified and mutated gets its own colour rather than two overlapping marks,
# because the overlap is the interesting case: it is where a disease variant
# lands on a regulatory site.
CATEGORY_COLOURS = {
    "ptm": "#F0C808",       # modifications
    "variant": "#D7263D",   # disease mutations
    "both": "#2E9E4F",      # modified and mutated
}

CATEGORY_LABELS = {
    "ptm": "Modified",
    "variant": "Disease mutation",
    "both": "Modified and mutated",
}

# Residue palette from the Scop3P app, kept for the per-residue colouring mode.
SITE_COLOURS = {
    "SER": "#1F77B4",
    "THR": "#FF7F0E",
    "TYR": "#2CA02C",
    "PTM": "#7D5BA6",
    "VARIANT": "#7B241C",
}


@dataclass
class Site:
    position: int
    kind: str              # "ptm" | "variant"
    residue: str = ""
    name: str = ""
    source: str = ""
    evidence: str = ""
    detail: str = ""
    score: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "position": self.position,
            "kind": self.kind,
            "residue": self.residue,
            "name": self.name,
            "source": self.source,
            "evidence": self.evidence,
            "detail": self.detail,
            "score": self.score,
            "colour": self.colour(),
        }

    def colour(self) -> str:
        """Category colour. The "both" case is resolved later, once every site
        is known, since it depends on what else sits on the same residue."""
        return CATEGORY_COLOURS.get(self.kind, CATEGORY_COLOURS["ptm"])

    def residue_colour(self) -> str:
        """Per-residue colouring, as the Scop3P app uses."""
        if self.kind == "variant":
            return SITE_COLOURS["VARIANT"]
        return SITE_COLOURS.get((self.residue or "").upper(), SITE_COLOURS["PTM"])


@dataclass
class StructureRef:
    pdb_id: str
    chains: Dict[str, Optional[Tuple[int, int]]] = field(default_factory=dict)
    method: str = ""
    resolution: Optional[float] = None
    coverage: Optional[float] = None

    def label(self) -> str:
        bits = [self.pdb_id.upper()]
        if self.method:
            bits.append(self.method)
        if self.resolution:
            bits.append(f"{self.resolution:.2f} A")
        if self.coverage:
            bits.append(f"{self.coverage * 100:.0f}% cover")
        return " \u00b7 ".join(bits)


# --------------------------------------------------------------------------
# transport
# --------------------------------------------------------------------------

class NotFound(ValueError):
    """The service answered, and the answer was that it has no such record."""


def _get(url: str, accept: str = "application/json") -> bytes:
    """Fetch a URL, preferring requests when it is installed.

    urllib with a custom User-Agent is refused outright by some services, which
    surfaces as a bare 403 that looks identical to the resource being missing.
    requests is what the Scop3P notebook uses and what these endpoints are known
    to accept, so it is used when present.
    """
    headers = {"Accept": accept}
    if USER_AGENT:
        headers["User-Agent"] = USER_AGENT
    host = urllib.parse.urlparse(url).netloc

    try:
        import requests  # noqa: PLC0415 - optional, checked at call time
    except ImportError:
        requests = None

    if requests is not None:
        try:
            response = requests.get(url, headers=headers, timeout=TIMEOUT)
        except Exception as error:  # noqa: BLE001 - network failures vary
            raise ValueError(f"Could not reach {host}: {error}") from error
        if response.status_code == 404:
            raise NotFound(f"{host} has no record ({response.status_code}).")
        if response.status_code != 200:
            raise ValueError(f"{host} returned HTTP {response.status_code}.")
        return response.content

    request = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            return response.read()
    except urllib.error.HTTPError as error:
        if error.code == 404:
            raise NotFound(f"{host} has no record (404).") from error
        raise ValueError(f"{host} returned HTTP {error.code}.") from error
    except urllib.error.URLError as error:
        raise ValueError(f"Could not reach {host}: {error.reason}") from error


class NotJson(ValueError):
    """The service answered with something that is not JSON."""


def _get_json(url: str) -> Any:
    """Fetch and decode JSON, describing the body when it is not JSON.

    A bare "Expecting value: line 1 column 1" says nothing about what actually
    came back. A refusal page, a redirect to a login screen and an empty body
    all produce that message, and they need different responses.
    """
    raw = _get(url)
    text = raw.decode("utf-8", errors="replace").strip()
    if not text:
        raise NotJson(f"{urllib.parse.urlparse(url).netloc} returned an empty body.")
    try:
        return json.loads(text)
    except json.JSONDecodeError as error:
        opening = text[:120].replace("\n", " ")
        looks_like = "an HTML page" if text[:1] == "<" else "non-JSON text"
        raise NotJson(
            f"{urllib.parse.urlparse(url).netloc} returned {looks_like} "
            f"instead of JSON: {opening!r}"
        ) from error


# --------------------------------------------------------------------------
# 1. which structures cover this accession
# --------------------------------------------------------------------------

def parse_uniprot_xrefs(payload: Dict[str, Any]) -> List[StructureRef]:
    """Read PDB cross-references out of a UniProtKB entry.

    The Chains property arrives in several shapes -- ``A=1-200``,
    ``A=1-200, B=5-150``, ``A/B=1-120``, and multi-range forms -- so the split
    is deliberately tolerant.
    """
    refs: List[StructureRef] = []
    for entry in payload.get("uniProtKBCrossReferences", []) or []:
        if entry.get("database") != "PDB":
            continue
        properties = {
            item.get("key"): item.get("value")
            for item in (entry.get("properties") or [])
            if isinstance(item, dict)
        }
        raw_chains = properties.get("Chains") or properties.get("Chain") or ""
        chains: Dict[str, Optional[Tuple[int, int]]] = {}

        for part in [p.strip() for p in re.split(r",\s*(?=[A-Za-z0-9/]+=)", raw_chains) if p.strip()]:
            if "=" not in part:
                continue
            names, ranges = part.split("=", 1)
            starts, ends = [], []
            for match in re.finditer(r"(-?\d+)\s*-\s*(-?\d+)", ranges):
                starts.append(int(match.group(1)))
                ends.append(int(match.group(2)))
            for name in [c.strip().upper()[:1] for c in re.split(r"[/\s]+", names) if c.strip()]:
                chains[name] = (min(starts), max(ends)) if starts and ends else None

        resolution = None
        raw_resolution = properties.get("Resolution") or ""
        match = re.search(r"([\d.]+)", str(raw_resolution))
        if match:
            try:
                resolution = float(match.group(1))
            except ValueError:
                resolution = None

        refs.append(StructureRef(
            pdb_id=(entry.get("id") or "").upper(),
            chains=chains,
            method=properties.get("Method") or "",
            resolution=resolution,
        ))
    return refs


def parse_best_structures(payload: Dict[str, Any], accession: str) -> List[StructureRef]:
    """Read PDBe's coverage-ranked structure list."""
    records = payload.get(accession.upper()) or payload.get(accession.lower()) or []
    refs: List[StructureRef] = []
    for record in records:
        pdb_id = (record.get("pdb_id") or "").upper()
        if not pdb_id:
            continue
        chain = record.get("chain_id")
        start = record.get("unp_start")
        end = record.get("unp_end")
        span = (int(start), int(end)) if start is not None and end is not None else None
        refs.append(StructureRef(
            pdb_id=pdb_id,
            chains={chain: span} if chain else {},
            method=record.get("experimental_method") or "",
            resolution=record.get("resolution"),
            coverage=record.get("coverage"),
        ))
    return refs


def merge_refs(ranked: List[StructureRef], xrefs: List[StructureRef]) -> List[StructureRef]:
    """Prefer PDBe's ranking, but keep UniProt's fuller chain listing.

    PDBe returns one row per chain ordered by coverage, which is the ordering a
    user wants; UniProt lists every chain of an entry in one row, which is what
    the chain dropdown needs. Neither alone is sufficient.
    """
    by_id: Dict[str, StructureRef] = {}
    order: List[str] = []

    for ref in ranked:
        if ref.pdb_id not in by_id:
            by_id[ref.pdb_id] = StructureRef(
                pdb_id=ref.pdb_id, chains=dict(ref.chains), method=ref.method,
                resolution=ref.resolution, coverage=ref.coverage,
            )
            order.append(ref.pdb_id)
        else:
            by_id[ref.pdb_id].chains.update(ref.chains)

    for ref in xrefs:
        if ref.pdb_id in by_id:
            for chain, span in ref.chains.items():
                by_id[ref.pdb_id].chains.setdefault(chain, span)
            if not by_id[ref.pdb_id].method:
                by_id[ref.pdb_id].method = ref.method
            if by_id[ref.pdb_id].resolution is None:
                by_id[ref.pdb_id].resolution = ref.resolution
        else:
            by_id[ref.pdb_id] = ref
            order.append(ref.pdb_id)

    return [by_id[pdb_id] for pdb_id in order]


def fetch_structures(accession: str) -> List[StructureRef]:
    """List PDB entries for an accession, best coverage first."""
    accession = accession.strip().upper()
    ranked: List[StructureRef] = []
    xrefs: List[StructureRef] = []

    try:
        ranked = parse_best_structures(_get_json(PDBE_BEST.format(accession=accession)), accession)
    except ValueError:
        ranked = []

    try:
        xrefs = parse_uniprot_xrefs(_get_json(UNIPROT_ENTRY.format(accession=accession)))
    except ValueError:
        xrefs = []

    return merge_refs(ranked, xrefs)


# --------------------------------------------------------------------------
# 2. numbering
# --------------------------------------------------------------------------

def parse_pdbe_mapping(
    payload: Dict[str, Any], pdb_id: str, accession: str, chain: Optional[str] = None
) -> Dict[int, int]:
    """Build UniProt position -> author residue number from PDBe SIFTS.

    Mappings arrive as aligned blocks, each carrying a UniProt span and the
    author span it corresponds to. Offsets differ per block, so they cannot be
    collapsed into a single shift.
    """
    entry = (payload.get(pdb_id.lower()) or payload.get(pdb_id.upper()) or {})
    records = (entry.get("UniProt") or {}).get(accession.upper()) or {}
    mapping: Dict[int, int] = {}

    for block in records.get("mappings", []) or []:
        block_chain = block.get("chain_id") or block.get("struct_asym_id")
        if chain and block_chain and block_chain.upper() != chain.upper():
            continue
        unp_start = block.get("unp_start")
        unp_end = block.get("unp_end")
        start = (block.get("start") or {}).get("residue_number")
        if unp_start is None or unp_end is None or start is None:
            continue
        for offset in range(int(unp_end) - int(unp_start) + 1):
            mapping[int(unp_start) + offset] = int(start) + offset

    return mapping


def fetch_numbering(pdb_id: str, accession: str, chain: Optional[str] = None) -> Dict[int, int]:
    """UniProt position -> author residue number for one entry and chain."""
    payload = _get_json(PDBE_UNIPROT.format(pdb_id=pdb_id.lower()))
    return parse_pdbe_mapping(payload, pdb_id, accession, chain)


def identity_numbering(positions: List[int]) -> Dict[int, int]:
    """AlphaFold models are UniProt-numbered already, so the map is identity."""
    return {position: position for position in positions}


# --------------------------------------------------------------------------
# 3. sites
# --------------------------------------------------------------------------

# "Phosphoserine" and friends name the modified residue; the three-letter code
# has to be read back out of them.
_MODIFIED_RESIDUE = {
    "phosphoserine": "SER", "phosphothreonine": "THR", "phosphotyrosine": "TYR",
    "phosphohistidine": "HIS", "phosphoaspartate": "ASP", "serine": "SER",
    "threonine": "THR", "tyrosine": "TYR",
}


def _residue_code(*candidates: Any) -> str:
    for candidate in candidates:
        text = str(candidate or "").strip()
        if not text:
            continue
        if len(text) == 3 and text.isalpha():
            return text.upper()
        lowered = text.lower()
        for needle, code in _MODIFIED_RESIDUE.items():
            if needle in lowered:
                return code
    return ""


def _normalise_score(*candidates: Any) -> Optional[float]:
    """Return a 0-1 score.

    The v1 API reports ``best_probability`` as a percentage while the older
    ``functionalScore`` was already a fraction. Mixing the two scales would make
    the filter slider mean different things for different proteins, so anything
    above 1 is treated as a percentage.
    """
    for candidate in candidates:
        if candidate in (None, ""):
            continue
        try:
            value = float(candidate)
        except (TypeError, ValueError):
            continue
        return value / 100.0 if value > 1.0 else value
    return None


def _join(*parts: Any) -> str:
    seen: List[str] = []
    for part in parts:
        text = str(part or "").strip()
        if text and text.lower() != "none" and text not in seen:
            seen.append(text)
    return " \u00b7 ".join(seen)


def parse_scop3p(payload: Any) -> List[Site]:
    """Read Scop3P modification records.

    Handles the v1 response (a bare array, ``uniprot_position`` /
    ``modification_name`` / ``modified_residue`` / ``best_probability``) and the
    older wrapped form, because a deployment may be on either.
    """
    if not payload:
        return []
    if isinstance(payload, dict):
        records = payload.get("modifications") or payload.get("data") or []
    else:
        records = payload

    sites: List[Site] = []
    for record in records:
        if not isinstance(record, dict):
            continue

        position = record.get("uniprot_position", record.get("position"))
        try:
            position = int(position)
        except (TypeError, ValueError):
            continue

        tissue = str(record.get("grouped_tissue") or "")
        # Values arrive as "PXD024548=Pancreas"; only the tissue is meaningful
        # here, since the project is reported separately.
        if "=" in tissue:
            tissue = tissue.split("=", 1)[1]

        pubmed = str(record.get("pubmed") or record.get("reference") or "")
        first_pubmed = pubmed.split(";")[0].strip() if pubmed else ""
        citation_count = len([p for p in pubmed.split(";") if p.strip()])
        citation = (
            f"PMID {first_pubmed}" + (f" +{citation_count - 1} more" if citation_count > 1 else "")
            if first_pubmed else ""
        )

        sites.append(Site(
            position=position,
            kind="ptm",
            residue=_residue_code(
                record.get("modified_residue"),
                record.get("residue"),
                record.get("uniprot_residue_annotation"),
            ),
            name=str(record.get("modification_name") or record.get("name") or "modification"),
            source="Scop3P",
            evidence=str(
                record.get("evidence_terms") or record.get("evidence")
                or record.get("evidence_code") or ""
            ),
            detail=_join(
                citation,
                tissue,
                str(record.get("source") or ""),
                record.get("best_project"),
            ),
            score=_normalise_score(
                record.get("best_probability"), record.get("functionalScore")
            ),
        ))
    return sites


_PTM_TYPES = {"MOD_RES", "Modified residue", "CARBOHYD", "Glycosylation",
              "LIPID", "Lipidation", "CROSSLNK"}

# Residue implied by the feature wording, since neither endpoint states it.
_DESC_RESIDUE = [
    ("phosphoserine", "SER"), ("phosphothreonine", "THR"),
    ("phosphotyrosine", "TYR"), ("phosphohistidine", "HIS"),
]


def _residue_from_description(description: str) -> str:
    text = (description or "").lower()
    for needle, code in _DESC_RESIDUE:
        if needle in text:
            return code
    return ""


def parse_uniprot_ptms(payload: Dict[str, Any]) -> List[Site]:
    """Read single-residue PTM features.

    Handles both shapes this data arrives in: the proteins API, which gives flat
    ``begin``/``end`` strings and a ``category``, and the UniProtKB entry
    endpoint, which nests positions under ``location``. Reading only one shape
    silently returns nothing against the other, which looks like a protein with
    no modifications rather than a parsing mismatch.
    """
    sites: List[Site] = []
    for feature in payload.get("features", []) or []:
        category = feature.get("category")
        feature_type = feature.get("type")
        if category is not None:
            if category != "PTM":
                continue
        elif feature_type not in _PTM_TYPES:
            continue

        begin = feature.get("begin")
        end = feature.get("end")
        if begin is None:
            location = feature.get("location") or {}
            begin = (location.get("start") or {}).get("value")
            end = (location.get("end") or {}).get("value")
        if begin is None or (end is not None and str(end) != str(begin)):
            continue

        try:
            position = int(begin)
        except (TypeError, ValueError):
            continue

        description = str(feature.get("description") or feature_type or "PTM")
        sites.append(Site(
            position=position,
            kind="ptm",
            residue=_residue_from_description(description),
            # UniProt appends qualifiers after a semicolon; the site name is
            # the part before it.
            name=description.split(";")[0].strip() or "modified residue",
            source="UniProt",
            evidence="UniProt feature",
        ))
    return sites


def parse_variants(payload: Dict[str, Any], accession: str) -> List[Site]:
    """Read disease-associated variants from the EBI proteins API.

    UniProt is the source rather than Scop3P, which imports its mutations from
    UniProt in the first place: going through Scop3P would be a second-hand copy
    that can only be equal to or staler than this one. Scop3P is used for
    modifications, where it holds data of its own.

    A feature carrying several disease associations becomes one site with the
    diseases joined, rather than several marks stacked on one residue.
    """
    sites: List[Site] = []
    for feature in payload.get("features", []) or []:
        if feature.get("type") != "VARIANT":
            continue
        diseases = [
            association.get("name")
            for association in feature.get("association", []) or []
            if association.get("disease") is True and association.get("name")
        ]
        if not diseases:
            continue
        try:
            position = int(feature.get("begin"))
        except (TypeError, ValueError):
            continue
        wild = feature.get("wildType") or ""
        mutant = feature.get("mutatedType") or ""
        sites.append(Site(
            position=position,
            kind="variant",
            name=f"{wild}{position}{mutant}" if wild and mutant else "variant",
            source="UniProt",
            evidence=feature.get("consequenceType") or "",
            detail="; ".join(sorted(set(diseases))),
            score=None,
        ))
    return sites


def merge_sites(*groups: List[Site]) -> List[Site]:
    """Combine site lists, collapsing duplicates at the same position and kind.

    Scop3P and UniProt overlap heavily. Two marks on one residue would read as
    two findings, so the Scop3P record wins on terminology and the UniProt
    evidence is folded into it.
    """
    merged: Dict[Tuple[int, str], Site] = {}
    for group in groups:
        for site in group:
            key = (site.position, site.kind)
            existing = merged.get(key)
            if existing is None:
                merged[key] = site
                continue
            if existing.source != "Scop3P" and site.source == "Scop3P":
                site.detail = "; ".join(filter(None, [site.detail, existing.detail]))
                merged[key] = site
            else:
                extra = [existing.detail, site.detail]
                existing.detail = "; ".join(dict.fromkeys(filter(None, extra)))
    return sorted(merged.values(), key=lambda item: (item.position, item.kind))


def _fetch_scop3p_modifications(accession: str) -> Tuple[List[Site], str]:
    """Try the v1 endpoint, then the older query-parameter form.

    Returns the sites and a note, so a protein Scop3P simply does not cover
    reads as a fact about the protein rather than a broken service.
    """
    attempts = [
        ("v1", SCOP3P_V1_MODS.format(accession=accession)),
        ("legacy", f"{SCOP3P_MODS}?{urllib.parse.urlencode({'accession': accession})}"),
    ]

    first_problem = ""
    for label, url in attempts:
        try:
            return parse_scop3p(_get_json(url)), ""
        except NotFound:
            first_problem = first_problem or (
                f"Scop3P has no entry for {accession} (it covers human phosphoproteins)."
            )
        except NotJson as error:
            first_problem = first_problem or f"Scop3P ({label}) answered but not with JSON: {error}"
        except ValueError as error:
            first_problem = first_problem or f"Scop3P ({label}) request failed: {error}"

    return [], first_problem


def fetch_ptms(accession: str, include_uniprot: bool = True) -> Tuple[List[Site], List[str]]:
    """Scop3P modifications, optionally merged with UniProt MOD_RES features."""
    accession = accession.strip().upper()
    notes: List[str] = []
    scop3p: List[Site] = []
    uniprot: List[Site] = []

    scop3p, scop3p_note = _fetch_scop3p_modifications(accession)
    if scop3p_note:
        notes.append(scop3p_note)

    if include_uniprot:
        try:
            url = (
                EBI_FEATURES.format(accession=accession)
                + "?" + urllib.parse.urlencode({"categories": "PTM"})
            )
            uniprot = parse_uniprot_ptms(_get_json(url))
        except ValueError as error:
            notes.append(f"UniProt PTM features unavailable ({error}).")

    sites = merge_sites(scop3p, uniprot)
    if not sites and not notes:
        notes.append(f"No PTMs recorded for {accession}.")
    return sites, notes


def fetch_variants(accession: str) -> Tuple[List[Site], List[str]]:
    """Disease-associated variants from the EBI proteins API."""
    accession = accession.strip().upper()
    try:
        sites = parse_variants(
            _get_json(UNIPROT_VARIANTS.format(accession=accession)), accession
        )
    except ValueError as error:
        return [], [f"Disease variants unavailable ({error})."]
    if not sites:
        return [], [f"No disease variants recorded for {accession}."]
    return sites, []


def fetch_sites(accession: str, include_variants: bool = True) -> Tuple[List[Site], List[str]]:
    """Gather PTMs and variants. Returns (sites, notes about what failed)."""
    accession = accession.strip().upper()
    notes: List[str] = []
    scop3p: List[Site] = []
    uniprot: List[Site] = []
    variants: List[Site] = []

    try:
        url = f"{SCOP3P_MODS}?{urllib.parse.urlencode({'accession': accession})}"
        scop3p = parse_scop3p(_get_json(url))
    except ValueError as error:
        # Scop3P covers human phosphoproteins; anything else legitimately misses.
        notes.append(f"Scop3P unavailable ({error}). UniProt features used instead.")

    try:
        uniprot = parse_uniprot_ptms(_get_json(UNIPROT_ENTRY.format(accession=accession)))
    except ValueError as error:
        notes.append(f"UniProt PTM features unavailable ({error}).")

    if include_variants:
        try:
            variants = parse_variants(
                _get_json(UNIPROT_VARIANTS.format(accession=accession)), accession
            )
        except ValueError as error:
            notes.append(f"Disease variants unavailable ({error}).")

    return merge_sites(scop3p, uniprot, variants), notes


# --------------------------------------------------------------------------
# projection onto the diagram
# --------------------------------------------------------------------------

def attach_sites(
    elements: List[Dict[str, Any]],
    sites: List[Site],
    numbering: Optional[Dict[int, int]] = None,
) -> Dict[str, Any]:
    """Place sites onto elements and summarise the result.

    ``numbering`` converts UniProt positions to the structure's own numbering.
    Positions with no mapping are counted as unmapped rather than guessed at: a
    residue absent from the model is a real fact about the structure, and
    silently shifting it onto a neighbour would be a fabricated finding.
    """
    mapped: List[Dict[str, Any]] = []
    unmapped: List[Dict[str, Any]] = []

    # A residue carrying both a modification and a disease mutation is the case
    # worth seeing, so it is classified before anything is drawn.
    modified = {site.position for site in sites if site.kind == "ptm"}
    mutated = {site.position for site in sites if site.kind == "variant"}
    overlapping = modified & mutated

    for site in sites:
        record = site.to_dict()
        record["uniprot_position"] = site.position
        category = "both" if site.position in overlapping else site.kind
        record["category"] = category
        record["category_label"] = CATEGORY_LABELS[category]
        record["colour"] = CATEGORY_COLOURS[category]
        record["residue_colour"] = site.residue_colour()
        if numbering is not None:
            target = numbering.get(site.position)
            if target is None:
                unmapped.append(record)
                continue
            record["position"] = target
        mapped.append(record)

    by_element: Dict[str, List[Dict[str, Any]]] = {}
    for element in elements:
        element["sites"] = []
        element["site_counts"] = {"ptm": 0, "variant": 0}

    lookup = {
        element["id"]: element for element in elements
    }
    spans = [(element["start"], element["stop"], element["id"]) for element in elements]

    loose: List[Dict[str, Any]] = []
    for record in mapped:
        position = record["position"]
        placed = False
        for start, stop, element_id in spans:
            if start <= position <= stop:
                element = lookup[element_id]
                # Fraction along the element, used to position the mark.
                span = max(1, stop - start)
                record = dict(record)
                record["t"] = (position - start) / span
                element["sites"].append(record)
                element["site_counts"][record["kind"]] += 1
                by_element.setdefault(element_id, []).append(record)
                placed = True
                break
        if not placed:
            loose.append(record)

    counts = {"ptm": 0, "variant": 0}
    for record in mapped:
        counts[record["kind"]] = counts.get(record["kind"], 0) + 1

    categories = {"ptm": 0, "variant": 0, "both": 0}
    seen_positions: Dict[int, str] = {}
    for record in mapped:
        seen_positions[record["position"]] = record["category"]
    for category in seen_positions.values():
        categories[category] = categories.get(category, 0) + 1

    max_density = max(
        (len(element["sites"]) for element in elements), default=0
    )

    return {
        "counts": counts,
        "categories": categories,
        "category_colours": CATEGORY_COLOURS,
        "category_labels": CATEGORY_LABELS,
        "total": len(mapped),
        "in_coil": len(loose),
        "unmapped": len(unmapped),
        "max_density": max_density,
        "coil_sites": loose,
        "colours": SITE_COLOURS,
    }


PDBE_FILE = "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_updated.cif"


def structure_file_url(pdb_id: str) -> str:
    """The updated mmCIF for a PDB entry.

    The ``_updated`` file is the one carrying PDBe's SIFTS cross-references, so
    it stays consistent with the numbering map fetched alongside it.
    """
    return PDBE_FILE.format(pdb_id=pdb_id.lower())


def fetch_structure_file(pdb_id: str) -> str:
    return _get(structure_file_url(pdb_id), accept="text/plain").decode("utf-8", errors="replace")


def probe(accession: str = "P07949") -> None:
    """Report what each endpoint actually returns, one line per service.

    When annotations come back empty there are several possible reasons -- the
    service is unreachable, it refused the client, it has no record for this
    protein, or the response parsed to nothing -- and they need different fixes.
    This distinguishes them instead of collapsing them into "unavailable".
    """
    try:
        import requests
        transport = f"requests {requests.__version__}"
    except ImportError:
        transport = "urllib (requests not installed)"
    print(f"transport   : {transport}")
    print(f"accession   : {accession}\n")

    checks = [
        ("Scop3P v1 PTMs",
         SCOP3P_V1_MODS.format(accession=accession),
         lambda data: f"{len(parse_scop3p(data))} sites"),
        ("Scop3P legacy",
         f"{SCOP3P_MODS}?{urllib.parse.urlencode({'accession': accession})}",
         lambda data: f"{len(parse_scop3p(data))} sites"),
        ("UniProt PTMs",
         EBI_FEATURES.format(accession=accession) + "?categories=PTM",
         lambda data: f"{len(parse_uniprot_ptms(data))} sites"),
        ("Disease variants",
         UNIPROT_VARIANTS.format(accession=accession),
         lambda data: f"{len(parse_variants(data, accession))} disease variants"),
        ("PDBe structures",
         PDBE_BEST.format(accession=accession),
         lambda data: f"{len(parse_best_structures(data, accession))} entries"),
        ("UniProt xrefs",
         UNIPROT_ENTRY.format(accession=accession),
         lambda data: f"{len(parse_uniprot_xrefs(data))} PDB cross-references"),
    ]

    for label, url, summarise in checks:
        try:
            payload = _get_json(url)
        except NotFound as error:
            print(f"{label:17s} no record  -- {error}")
            _describe(url)
            continue
        except NotJson as error:
            print(f"{label:17s} NOT JSON   -- {error}")
            _describe(url)
            continue
        except ValueError as error:
            print(f"{label:17s} FAILED     -- {error}")
            _describe(url)
            continue
        try:
            print(f"{label:17s} ok         -- {summarise(payload)}")
        except Exception as error:  # noqa: BLE001
            shape = type(payload).__name__
            keys = list(payload)[:6] if isinstance(payload, dict) else "n/a"
            print(f"{label:17s} PARSE MISS -- {error} (got {shape}, keys={keys})")

    print("\nFAILED     the service refused the request or was unreachable")
    print("NOT JSON   something came back, but not JSON -- see the raw reply above")
    print("PARSE MISS JSON arrived, but its shape is not what the parser expects")


def _describe(url: str) -> None:
    """Print the raw reply, which is what actually identifies the problem."""
    try:
        import requests
        response = requests.get(url, timeout=TIMEOUT)
        body = (response.text or "").strip()
        print(f"                  status {response.status_code}, "
              f"content-type {response.headers.get('content-type', 'unknown')}, "
              f"{len(body)} chars")
        if body:
            print(f"                  starts: {body[:160]!r}")
    except Exception as error:  # noqa: BLE001
        print(f"                  raw probe also failed: {error}")
