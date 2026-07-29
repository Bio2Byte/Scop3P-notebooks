"""Secondary-structure elements and the contacts between them.

Turns a parsed chain into the element list the layout engines consume, and
works out which strands pair into sheets.  Contact detection is purely
geometric -- CA-CA proximity plus a direction test -- so it behaves identically
whether the secondary structure came from the file or was derived from
coordinates.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .io import Residue, SSRange, Structure

# Two strands are paired when this many residue pairs sit within the cutoff.
_CONTACT_CUTOFF = 7.0
_MIN_SEQUENCE_SEPARATION = 3
_MIN_CONTACT_PAIRS = 2


def _distance(a: Residue, b: Residue) -> float:
    return math.sqrt(
        (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
    )


def _element_vector(
    element: Dict[str, Any], residue_by_seq: Dict[int, Residue]
) -> Tuple[float, float, float]:
    start = residue_by_seq.get(element["start"])
    stop = residue_by_seq.get(element["stop"])
    if not start or not stop:
        return (0.0, 1.0, 0.0)
    return (stop.x - start.x, stop.y - start.y, stop.z - start.z)


def _centroid(
    element: Dict[str, Any], residue_by_seq: Dict[int, Residue]
) -> Tuple[float, float, float]:
    members = [
        residue_by_seq[seq]
        for seq in range(element["start"], element["stop"] + 1)
        if seq in residue_by_seq
    ]
    if not members:
        return (0.0, 0.0, 0.0)
    count = float(len(members))
    return (
        sum(residue.x for residue in members) / count,
        sum(residue.y for residue in members) / count,
        sum(residue.z for residue in members) / count,
    )


def build_elements(
    structure: Structure, chain: str
) -> Tuple[List[Dict[str, Any]], List[Residue]]:
    """Produce ordered helix/strand elements for one chain."""
    residues = structure.residues_by_chain.get(chain, [])
    if not residues:
        raise ValueError(f"Chain {chain} has no residues.")

    residue_by_seq = {residue.seq: residue for residue in residues}
    ranges: Sequence[SSRange] = sorted(
        structure.ss_by_chain.get(chain, []), key=lambda entry: (entry.start, entry.stop)
    )

    elements: List[Dict[str, Any]] = []
    counters = {"helix": 0, "strand": 0}

    for entry in ranges:
        members = [
            residue_by_seq[seq]
            for seq in range(entry.start, entry.stop + 1)
            if seq in residue_by_seq
        ]
        if not members:
            continue
        counters[entry.kind] += 1
        element_id = ("H" if entry.kind == "helix" else "S") + str(counters[entry.kind])
        confidences = [r.plddt for r in members if r.plddt is not None]
        elements.append(
            {
                "id": element_id,
                "type": entry.kind,
                "chain": chain,
                "start": members[0].seq,
                "stop": members[-1].seq,
                "length": len(members),
                "sequence": "".join(residue.aa for residue in members),
                "ss_code": "H" if entry.kind == "helix" else "E",
                "ss_name": "Alpha helix" if entry.kind == "helix" else "Beta strand",
                "residue_numbers": [residue.seq for residue in members],
                "residue_ids": [residue.residue_id for residue in members],
                "label_seq_numbers": [residue.label_seq for residue in members],
                "confidence_mean": (
                    round(sum(confidences) / len(confidences), 2) if confidences else None
                ),
            }
        )

    elements.sort(key=lambda element: (element["start"], element["stop"]))
    return elements, residues


def strand_contacts(
    elements: List[Dict[str, Any]], residues: List[Residue]
) -> List[Dict[str, Any]]:
    """Find beta-sheet pairings between strands.

    Orientation comes from the dot product of the two strand vectors, which is
    a coarse test but agrees with hydrogen-bond topology for ordinary sheets.
    """
    residue_by_seq = {residue.seq: residue for residue in residues}
    strands = [element for element in elements if element["type"] == "strand"]
    contacts: List[Dict[str, Any]] = []

    for index, left in enumerate(strands):
        left_residues = [
            residue_by_seq[seq]
            for seq in range(left["start"], left["stop"] + 1)
            if seq in residue_by_seq
        ]
        for right in strands[index + 1 :]:
            right_residues = [
                residue_by_seq[seq]
                for seq in range(right["start"], right["stop"] + 1)
                if seq in residue_by_seq
            ]

            pairs: List[Tuple[int, int]] = []
            closest: Optional[float] = None
            for residue_a in left_residues:
                for residue_b in right_residues:
                    if abs(residue_a.seq - residue_b.seq) < _MIN_SEQUENCE_SEPARATION:
                        continue
                    separation = _distance(residue_a, residue_b)
                    if closest is None or separation < closest:
                        closest = separation
                    if separation <= _CONTACT_CUTOFF:
                        pairs.append((residue_a.seq, residue_b.seq))

            if len(pairs) < _MIN_CONTACT_PAIRS:
                continue

            left_vector = _element_vector(left, residue_by_seq)
            right_vector = _element_vector(right, residue_by_seq)
            alignment = sum(a * b for a, b in zip(left_vector, right_vector))

            contacts.append(
                {
                    "source": left["id"],
                    "target": right["id"],
                    "count": len(pairs),
                    "orientation": "parallel" if alignment >= 0 else "antiparallel",
                    "min_distance": round(closest or 0.0, 2),
                    "pairs": pairs[:8],
                }
            )

    return sorted(
        contacts, key=lambda item: (-item["count"], item["source"], item["target"])
    )


def annotate_geometry(
    elements: List[Dict[str, Any]], residues: List[Residue]
) -> None:
    """Attach 3D centroid and axis vector to each element, in place.

    The sheet layout uses these to place helices by where they actually sit in
    space rather than by where they fall in the sequence.
    """
    residue_by_seq = {residue.seq: residue for residue in residues}
    for element in elements:
        element["centroid"] = _centroid(element, residue_by_seq)
        element["axis"] = _element_vector(element, residue_by_seq)


def assign_sheets(
    elements: List[Dict[str, Any]], contacts: List[Dict[str, Any]]
) -> None:
    """Label each strand with the sheet it belongs to, in place.

    Which strands pair into which sheet is the single most useful thing a
    topology diagram can say beyond helix-versus-strand, and it is already
    implied by the contact graph. Surfacing it as a label lets the renderer
    colour by sheet, so a reader can see the sheets without tracing every
    pairing line by hand.
    """
    strand_ids = [e["id"] for e in elements if e["type"] == "strand"]
    known = set(strand_ids)

    adjacency: Dict[str, List[str]] = {item: [] for item in strand_ids}
    for contact in contacts:
        if contact["source"] in known and contact["target"] in known:
            adjacency[contact["source"]].append(contact["target"])
            adjacency[contact["target"]].append(contact["source"])

    by_id = {element["id"]: element for element in elements}
    seen: set = set()
    groups: List[List[str]] = []

    for strand_id in strand_ids:
        if strand_id in seen:
            continue
        stack = [strand_id]
        seen.add(strand_id)
        group: List[str] = []
        while stack:
            current = stack.pop()
            group.append(current)
            for neighbour in adjacency[current]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        groups.append(group)

    # Number sheets in order of first appearance along the chain, so the labels
    # read in the same direction as the sequence.
    groups.sort(key=lambda group: min(by_id[item]["start"] for item in group))

    for element in elements:
        element["sheet"] = None
        element["sheet_index"] = None

    for index, group in enumerate(groups, start=1):
        # A lone strand pairs with nothing, so calling it a sheet would be a lie.
        label = f"B{index}" if len(group) > 1 else None
        for item in group:
            by_id[item]["sheet"] = label
            by_id[item]["sheet_index"] = index if label else None
