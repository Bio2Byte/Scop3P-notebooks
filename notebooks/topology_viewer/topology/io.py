"""Structure loading, normalised across formats.

Everything upstream of the topology code funnels through :func:`load_structure`,
which flattens four input shapes into one :class:`Structure`:

    fetched AlphaFold mmCIF      annotated, UniProt-numbered
    fetched or uploaded PDB      annotated, author-numbered
    uploaded prediction mmCIF    usually no secondary structure at all
    uploaded prediction PDB      usually no secondary structure at all

The extension is a hint, never the decision: plenty of tools emit PDB-format
text under a ``.cif`` name.  Content sniffing wins.

Two numbering schemes are kept side by side.  ``seq`` is author numbering,
which is what a user reads off a paper and types into a search box, and
``label_seq`` is the mmCIF entity numbering, which is what Mol* selection
queries want.  For AlphaFold models they coincide; for real PDB entries they
routinely do not, and conflating them is the classic source of off-by-N
highlighting bugs.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from . import ss as ss_module

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H",
}


@dataclass
class Residue:
    chain: str
    seq: int
    comp_id: str
    aa: str
    x: float
    y: float
    z: float
    label_chain: Optional[str] = None
    label_seq: Optional[int] = None
    ins_code: str = ""
    bfactor: Optional[float] = None
    plddt: Optional[float] = None
    # SIFTS cross-reference carried in PDBe's _updated.cif files.
    uniprot_acc: Optional[str] = None
    uniprot_seq: Optional[int] = None

    @property
    def residue_id(self) -> str:
        return f"{self.chain}:{self.seq}{self.ins_code}".strip()

    @property
    def coords(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass
class SSRange:
    chain: str
    start: int
    stop: int
    kind: str  # "helix" | "strand"


@dataclass
class Structure:
    name: str
    fmt: str
    residues_by_chain: Dict[str, List[Residue]] = field(default_factory=dict)
    ss_by_chain: Dict[str, List[SSRange]] = field(default_factory=dict)
    ss_source: str = "none"
    has_sifts: bool = False
    entry_id: str = ""
    title: str = ""
    uniprot: str = ""
    has_plddt: bool = False
    raw_text: str = ""

    @property
    def chains(self) -> List[str]:
        return sorted(self.residues_by_chain)

    def chain_options(self) -> List[Tuple[str, str]]:
        """Dropdown entries, largest chain first so the default is the interesting one."""
        ordered = sorted(
            self.residues_by_chain.items(),
            key=lambda item: (-len(item[1]), item[0]),
        )
        return [
            (f"{chain} ({len(residues)} residues)", chain)
            for chain, residues in ordered
        ]

    def sifts_numbering(
        self, chain: str, accession: str = ""
    ) -> Dict[int, int]:
        """UniProt position -> author residue number, read from the file itself.

        Preferred over the PDBe mapping API: it needs no second request, and it
        cannot disagree with the coordinates being displayed, because it came
        from the same file. Returns an empty dict when the file carries no SIFTS
        columns, which is the signal to fall back to the API.
        """
        mapping: Dict[int, int] = {}
        wanted = (accession or "").strip().upper()
        for residue in self.residues_by_chain.get(chain, []):
            if residue.uniprot_seq is None:
                continue
            if wanted and residue.uniprot_acc and residue.uniprot_acc.upper() != wanted:
                continue
            mapping[residue.uniprot_seq] = residue.seq
        return mapping

    def sifts_accessions(self) -> List[str]:
        """Every UniProt accession this structure cross-references."""
        found = {
            residue.uniprot_acc
            for residues in self.residues_by_chain.values()
            for residue in residues
            if residue.uniprot_acc
        }
        return sorted(found)

    def chains_for_accession(self, accession: str) -> List[str]:
        """Chains that actually map to this accession."""
        wanted = (accession or "").strip().upper()
        if not wanted:
            return []
        return [
            chain for chain, residues in self.residues_by_chain.items()
            if any(
                r.uniprot_acc and r.uniprot_acc.upper() == wanted for r in residues
            )
        ]

    def default_chain(self) -> str:
        if not self.residues_by_chain:
            return ""
        return max(self.residues_by_chain, key=lambda c: len(self.residues_by_chain[c]))


# --------------------------------------------------------------------------
# format sniffing
# --------------------------------------------------------------------------

def sniff_format(text: str, name: str = "") -> str:
    """Decide mmCIF vs PDB from content, falling back to the extension."""
    head = text[:20000]
    if re.search(r"^\s*data_", head, re.MULTILINE):
        return "mmcif"
    if re.search(r"^(loop_|_atom_site\.)", head, re.MULTILINE):
        return "mmcif"
    if re.search(r"^(ATOM  |HETATM|HEADER|MODEL |CRYST1|EXPDTA)", head, re.MULTILINE):
        return "pdb"
    suffix = Path(name).suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return "mmcif"
    if suffix in {".pdb", ".ent"}:
        return "pdb"
    return "unknown"


# --------------------------------------------------------------------------
# mmCIF
# --------------------------------------------------------------------------

def _cif_tokens(text: str) -> List[str]:
    """Tokenise mmCIF, honouring quotes and semicolon text blocks."""
    tokens: List[str] = []
    for raw_line in text.splitlines():
        if raw_line.startswith(";"):
            # Multi-line value; the opening line carries the first chunk.
            tokens.append(raw_line[1:].strip())
            continue
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        position = 0
        length = len(line)
        while position < length:
            char = line[position]
            if char.isspace():
                position += 1
                continue
            if char in "'\"":
                close = line.find(char, position + 1)
                while close != -1 and close + 1 < length and not line[close + 1].isspace():
                    close = line.find(char, close + 1)
                if close == -1:
                    tokens.append(line[position + 1 :])
                    break
                tokens.append(line[position + 1 : close])
                position = close + 1
                continue
            end = position
            while end < length and not line[end].isspace():
                end += 1
            tokens.append(line[position:end])
            position = end
    return tokens


def _parse_cif(text: str) -> Tuple[Dict[str, List[Dict[str, str]]], Dict[str, str]]:
    """Return loop rows keyed by category, plus standalone key/value items.

    Single-copy categories are written as plain key/value pairs rather than
    loops, so a protein with exactly one helix has no ``loop_`` for
    ``_struct_conf`` at all.  Both shapes have to be read.
    """
    tokens = _cif_tokens(text)
    loops: Dict[str, List[Dict[str, str]]] = {}
    items: Dict[str, str] = {}

    index = 0
    total = len(tokens)
    while index < total:
        token = tokens[index]

        if token.lower() == "loop_":
            index += 1
            headers: List[str] = []
            while index < total and tokens[index].startswith("_"):
                headers.append(tokens[index])
                index += 1
            if not headers:
                continue
            category = headers[0].split(".", 1)[0]
            width = len(headers)
            while index < total:
                peek = tokens[index]
                if peek.startswith("_") or peek.lower() in {"loop_"} or peek.lower().startswith("data_"):
                    break
                row = tokens[index : index + width]
                if len(row) < width:
                    break
                loops.setdefault(category, []).append(dict(zip(headers, row)))
                index += width
            continue

        if token.startswith("_") and index + 1 < total and not tokens[index + 1].startswith("_"):
            items[token] = tokens[index + 1]
            index += 2
            continue

        index += 1

    # Present standalone items as a one-row loop so callers have a single shape.
    grouped: Dict[str, Dict[str, str]] = {}
    for key, value in items.items():
        if "." not in key:
            continue
        category = key.split(".", 1)[0]
        grouped.setdefault(category, {})[key] = value
    for category, row in grouped.items():
        loops.setdefault(category, []).append(row)

    return loops, items


def _to_int(value: Any) -> Optional[int]:
    if value in {None, "", ".", "?"}:
        return None
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def _to_float(value: Any) -> Optional[float]:
    if value in {None, "", ".", "?"}:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _conf_kind(conf_type: str) -> Optional[str]:
    value = (conf_type or "").upper()
    if value.startswith("HELX"):
        return "helix"
    if value.startswith("STRN") or value.startswith("SHEET"):
        return "strand"
    return None


def parse_mmcif(text: str, name: str) -> Structure:
    loops, items = _parse_cif(text)
    structure = Structure(name=name, fmt="mmcif", raw_text=text)

    structure.entry_id = (
        items.get("_entry.id")
        or (loops.get("_entry", [{}])[0].get("_entry.id") if loops.get("_entry") else "")
        or Path(name).stem
    )
    for row in loops.get("_struct", []):
        structure.title = row.get("_struct.title", "") or structure.title

    # AlphaFold stores per-residue confidence separately from the B-factor.
    plddt_by_label_seq: Dict[int, float] = {}
    for row in loops.get("_ma_qa_metric_local", []):
        seq = _to_int(row.get("_ma_qa_metric_local.label_seq_id"))
        value = _to_float(row.get("_ma_qa_metric_local.metric_value"))
        if seq is not None and value is not None:
            plddt_by_label_seq[seq] = value

    for row in loops.get("_struct_ref", []):
        database = (row.get("_struct_ref.db_name") or "").upper()
        if database in {"UNP", "UNIPROT"}:
            structure.uniprot = row.get("_struct_ref.pdbx_db_accession", "") or structure.uniprot

    atom_rows = loops.get("_atom_site", [])
    if not atom_rows:
        raise ValueError(
            "No atom records found in this mmCIF. If the file came from a viewer export, "
            "re-save it with coordinates included."
        )

    seen: set = set()
    first_model: Optional[str] = None
    for row in atom_rows:
        model = row.get("_atom_site.pdbx_PDB_model_num")
        if first_model is None:
            first_model = model
        if model != first_model:
            continue

        atom_name = row.get("_atom_site.label_atom_id") or row.get("_atom_site.auth_atom_id")
        if atom_name != "CA":
            continue

        altloc = row.get("_atom_site.label_alt_id", ".")
        if altloc not in {".", "?", "", "A"}:
            continue

        chain = (
            row.get("_atom_site.auth_asym_id")
            or row.get("_atom_site.label_asym_id")
            or "A"
        )
        label_chain = row.get("_atom_site.label_asym_id") or chain
        seq = _to_int(row.get("_atom_site.auth_seq_id"))
        label_seq = _to_int(row.get("_atom_site.label_seq_id"))
        if seq is None:
            seq = label_seq
        if seq is None:
            continue

        x = _to_float(row.get("_atom_site.Cartn_x"))
        y = _to_float(row.get("_atom_site.Cartn_y"))
        z = _to_float(row.get("_atom_site.Cartn_z"))
        if x is None or y is None or z is None:
            continue

        ins = row.get("_atom_site.pdbx_PDB_ins_code", ".")
        ins_code = "" if ins in {".", "?", None} else str(ins)

        key = (chain, seq, ins_code)
        if key in seen:
            continue
        seen.add(key)

        # PDBe's updated mmCIF carries the SIFTS mapping per atom, so the
        # UniProt correspondence arrives with the coordinates instead of needing
        # a second request that could disagree with them.
        sifts_acc = row.get("_atom_site.pdbx_sifts_xref_db_acc")
        sifts_num = _to_int(row.get("_atom_site.pdbx_sifts_xref_db_num"))
        if sifts_acc in {".", "?", ""}:
            sifts_acc = None

        comp = row.get("_atom_site.label_comp_id") or row.get("_atom_site.auth_comp_id") or "UNK"
        bfactor = _to_float(row.get("_atom_site.B_iso_or_equiv"))
        plddt = plddt_by_label_seq.get(label_seq) if label_seq is not None else None

        structure.residues_by_chain.setdefault(chain, []).append(
            Residue(
                chain=chain,
                seq=seq,
                comp_id=comp,
                aa=AA3_TO_1.get(comp.upper(), "X"),
                x=x, y=y, z=z,
                label_chain=label_chain,
                label_seq=label_seq,
                ins_code=ins_code,
                bfactor=bfactor,
                plddt=plddt if plddt is not None else bfactor,
                uniprot_acc=sifts_acc,
                uniprot_seq=sifts_num,
            )
        )

    structure.has_sifts = any(
        residue.uniprot_seq is not None
        for residues in structure.residues_by_chain.values()
        for residue in residues
    )
    structure.has_plddt = bool(plddt_by_label_seq) or _looks_like_plddt(structure)

    # Secondary structure as annotated, if the writer bothered.
    label_to_auth: Dict[Tuple[str, int], int] = {}
    for chain, residues in structure.residues_by_chain.items():
        for residue in residues:
            if residue.label_seq is not None:
                label_to_auth[(residue.label_chain or chain, residue.label_seq)] = residue.seq

    ranges: List[SSRange] = []
    for row in loops.get("_struct_conf", []):
        kind = _conf_kind(row.get("_struct_conf.conf_type_id", ""))
        if kind is None:
            continue
        chain = (
            row.get("_struct_conf.beg_auth_asym_id")
            or row.get("_struct_conf.beg_label_asym_id")
        )
        start = _to_int(row.get("_struct_conf.beg_auth_seq_id"))
        stop = _to_int(row.get("_struct_conf.end_auth_seq_id"))
        if start is None or stop is None:
            label_chain = row.get("_struct_conf.beg_label_asym_id")
            label_start = _to_int(row.get("_struct_conf.beg_label_seq_id"))
            label_stop = _to_int(row.get("_struct_conf.end_label_seq_id"))
            if label_chain and label_start is not None and label_stop is not None:
                start = label_to_auth.get((label_chain, label_start))
                stop = label_to_auth.get((label_chain, label_stop))
                chain = chain or label_chain
        if not chain or start is None or stop is None or stop < start:
            continue
        ranges.append(SSRange(chain=chain, start=start, stop=stop, kind=kind))

    # Sheets live in their own category rather than _struct_conf.
    for row in loops.get("_struct_sheet_range", []):
        chain = (
            row.get("_struct_sheet_range.beg_auth_asym_id")
            or row.get("_struct_sheet_range.beg_label_asym_id")
        )
        start = _to_int(row.get("_struct_sheet_range.beg_auth_seq_id"))
        stop = _to_int(row.get("_struct_sheet_range.end_auth_seq_id"))
        if not chain or start is None or stop is None or stop < start:
            continue
        ranges.append(SSRange(chain=chain, start=start, stop=stop, kind="strand"))

    for entry in ranges:
        structure.ss_by_chain.setdefault(entry.chain, []).append(entry)
    if ranges:
        structure.ss_source = "file (_struct_conf)"

    return structure


# --------------------------------------------------------------------------
# PDB
# --------------------------------------------------------------------------

def parse_pdb(text: str, name: str) -> Structure:
    structure = Structure(name=name, fmt="pdb", raw_text=text)
    structure.entry_id = Path(name).stem

    ranges: List[SSRange] = []
    seen: set = set()
    in_later_model = False

    for line in text.splitlines():
        record = line[:6]

        if record == "MODEL ":
            index = line[10:14].strip()
            in_later_model = index not in {"", "1"}
            continue
        if record == "ENDMDL":
            in_later_model = True
            continue
        if in_later_model:
            continue

        if record == "HEADER":
            structure.title = line[10:50].strip() or structure.title
            candidate = line[62:66].strip()
            if candidate:
                structure.entry_id = candidate
            continue

        if record == "TITLE ":
            structure.title = (structure.title + " " + line[10:80].strip()).strip()
            continue

        if record == "DBREF ":
            if line[26:32].strip().upper() in {"UNP", "UNIPROT"}:
                structure.uniprot = line[33:41].strip() or structure.uniprot
            continue

        if record == "HELIX ":
            chain = line[19:20].strip()
            start = _to_int(line[21:25].strip())
            stop = _to_int(line[33:37].strip())
            if chain and start is not None and stop is not None and stop >= start:
                ranges.append(SSRange(chain=chain, start=start, stop=stop, kind="helix"))
            continue

        if record == "SHEET ":
            chain = line[21:22].strip()
            start = _to_int(line[22:26].strip())
            stop = _to_int(line[33:37].strip())
            if chain and start is not None and stop is not None and stop >= start:
                ranges.append(SSRange(chain=chain, start=start, stop=stop, kind="strand"))
            continue

        if record not in {"ATOM  ", "HETATM"}:
            continue

        if line[12:16].strip() != "CA":
            continue

        altloc = line[16:17].strip()
        if altloc not in {"", "A"}:
            continue

        comp = line[17:20].strip() or "UNK"
        if record == "HETATM" and comp.upper() not in AA3_TO_1:
            continue

        chain = line[21:22].strip() or "A"
        seq = _to_int(line[22:26].strip())
        if seq is None:
            continue
        ins_code = line[26:27].strip()

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue

        key = (chain, seq, ins_code)
        if key in seen:
            continue
        seen.add(key)

        bfactor = _to_float(line[60:66].strip()) if len(line) >= 66 else None

        structure.residues_by_chain.setdefault(chain, []).append(
            Residue(
                chain=chain,
                seq=seq,
                comp_id=comp,
                aa=AA3_TO_1.get(comp.upper(), "X"),
                x=x, y=y, z=z,
                label_chain=chain,
                label_seq=seq,
                ins_code=ins_code,
                bfactor=bfactor,
                plddt=bfactor,
            )
        )

    if not structure.residues_by_chain:
        raise ValueError(
            "No CA atoms found in this PDB file. Check that it contains ATOM records "
            "rather than only header or connectivity lines."
        )

    for entry in ranges:
        structure.ss_by_chain.setdefault(entry.chain, []).append(entry)
    if ranges:
        structure.ss_source = "file (HELIX/SHEET)"

    structure.has_plddt = _looks_like_plddt(structure)
    return structure


def _looks_like_plddt(structure: Structure) -> bool:
    """Guess whether the B-factor column is really per-residue confidence.

    Prediction tools reuse the B-factor slot for pLDDT. Values confined to
    0-100 with a high mean are the tell; crystallographic B-factors spread
    wider and start lower. Only affects how the colour ramp is labelled.
    """
    values = [
        residue.bfactor
        for residues in structure.residues_by_chain.values()
        for residue in residues
        if residue.bfactor is not None
    ]
    if len(values) < 5:
        return False
    inside = sum(1 for value in values if 0.0 <= value <= 100.0)
    return inside == len(values) and (sum(values) / len(values)) > 40.0


# --------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------

def load_structure(text: str, name: str = "structure", *, prefer_ss: str = "auto") -> Structure:
    """Parse any supported coordinate file and guarantee a secondary structure.

    Raises ``ValueError`` with a message intended for display when the text is
    not a structure file at all.
    """
    if not text or not text.strip():
        raise ValueError("This file is empty.")

    fmt = sniff_format(text, name)
    if fmt == "mmcif":
        structure = parse_mmcif(text, name)
    elif fmt == "pdb":
        structure = parse_pdb(text, name)
    else:
        raise ValueError(
            "Unrecognised file. Upload a coordinate file in PDB or mmCIF format "
            "(.pdb, .ent, .cif, .mmcif)."
        )

    for chain, residues in structure.residues_by_chain.items():
        residues.sort(key=lambda residue: (residue.seq, residue.ins_code))

    ensure_secondary_structure(structure, prefer=prefer_ss)
    return structure


def ensure_secondary_structure(structure: Structure, *, prefer: str = "auto") -> None:
    """Fill in secondary structure for any chain the file left unannotated."""
    sources: set = set()

    for chain, residues in structure.residues_by_chain.items():
        annotated = structure.ss_by_chain.get(chain)
        if annotated:
            # Keep only ranges that actually land on observed residues.
            observed = {residue.seq for residue in residues}
            kept = [
                entry for entry in annotated
                if any(seq in observed for seq in range(entry.start, entry.stop + 1))
            ]
            if kept:
                structure.ss_by_chain[chain] = kept
                sources.add(structure.ss_source)
                continue

        coords = [residue.coords for residue in residues]
        codes, provenance = ss_module.compute(
            coords,
            pdb_text=structure.raw_text if structure.fmt == "pdb" else "",
            chain=chain,
            prefer=prefer,
        )
        sources.add(provenance)
        structure.ss_by_chain[chain] = _codes_to_ranges(chain, residues, codes)

    if sources:
        structure.ss_source = " / ".join(sorted(sources))


def _codes_to_ranges(chain: str, residues: List[Residue], codes: List[str]) -> List[SSRange]:
    """Collapse a per-residue code string into contiguous element ranges."""
    ranges: List[SSRange] = []
    current: Optional[str] = None
    start_index = 0

    for index, code in enumerate(codes):
        if code != current:
            if current in {"H", "E"}:
                ranges.append(
                    SSRange(
                        chain=chain,
                        start=residues[start_index].seq,
                        stop=residues[index - 1].seq,
                        kind="helix" if current == "H" else "strand",
                    )
                )
            current = code
            start_index = index

    if current in {"H", "E"} and codes:
        ranges.append(
            SSRange(
                chain=chain,
                start=residues[start_index].seq,
                stop=residues[len(codes) - 1].seq,
                kind="helix" if current == "H" else "strand",
            )
        )
    return ranges
