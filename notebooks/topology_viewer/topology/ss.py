"""Secondary-structure assignment.

Structures arrive from four places and only some of them carry secondary
structure:

    AlphaFold DB / RCSB mmCIF   -> _struct_conf records
    AlphaFold DB / RCSB PDB     -> HELIX / SHEET records
    ColabFold / Boltz / Chai    -> nothing; coordinates only
    hand-built or stripped      -> nothing

When the file says nothing we derive it from CA coordinates.  Three backends
are tried in order of fidelity; all of them return the same per-residue string
of ``H`` (helix), ``E`` (strand) and ``C`` (coil), so downstream code never
learns which one ran.

    mkdssp     reference Kabsch-Sander, needs the binary on PATH
    biotite    P-SEA, pip-installable, CA only
    builtin    P-SEA reimplemented here, no dependencies at all

The built-in path exists so the app keeps working on a bare Voila deployment.
Callers get the provenance back alongside the assignment and should surface it,
because P-SEA and DSSP routinely disagree about where an element stops by a
residue or two.
"""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

# P-SEA thresholds (Labesse et al. 1997).  Distances in angstrom between CA
# atoms i and i+n; angles in degrees.
_HELIX_D2 = (5.5, 0.5)
_HELIX_D3 = (5.3, 0.5)
_HELIX_D4 = (6.4, 0.6)
_HELIX_ANGLE = (89.0, 12.0)
_HELIX_DIHEDRAL = (50.0, 20.0)

_STRAND_D2 = (6.7, 0.6)
_STRAND_D3 = (9.9, 0.9)
_STRAND_D4 = (12.4, 1.1)
_STRAND_ANGLE = (124.0, 14.0)
_STRAND_DIHEDRAL = (-170.0, 45.0)

# Consecutive CA atoms sit ~3.8 A apart.  Anything past this is a chain break
# and no element may span it.
_CHAIN_BREAK = 4.5

# Elements below these lengths are noise rather than structure.
_MIN_HELIX = 4
_MIN_STRAND = 2

Point = Tuple[float, float, float]


def _within(value: Optional[float], spec: Tuple[float, float]) -> bool:
    if value is None:
        return False
    centre, tolerance = spec
    return abs(value - centre) <= tolerance


def _distance(a: Point, b: Point) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _subtract(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Point, b: Point) -> Point:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Point, b: Point) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a: Point) -> float:
    return math.sqrt(_dot(a, a))


def _angle(a: Point, b: Point, c: Point) -> Optional[float]:
    """Angle at ``b`` in degrees."""
    left = _subtract(a, b)
    right = _subtract(c, b)
    scale = _norm(left) * _norm(right)
    if scale < 1e-6:
        return None
    cosine = max(-1.0, min(1.0, _dot(left, right) / scale))
    return math.degrees(math.acos(cosine))


def _dihedral(a: Point, b: Point, c: Point, d: Point) -> Optional[float]:
    """Torsion about the ``b``-``c`` axis in degrees, signed, in (-180, 180]."""
    b1 = _subtract(b, a)
    b2 = _subtract(c, b)
    b3 = _subtract(d, c)
    n1 = _cross(b1, b2)
    n2 = _cross(b2, b3)
    b2_len = _norm(b2)
    if b2_len < 1e-6 or _norm(n1) < 1e-6 or _norm(n2) < 1e-6:
        return None
    m1 = _cross(n1, (b2[0] / b2_len, b2[1] / b2_len, b2[2] / b2_len))
    x = _dot(n1, n2)
    y = _dot(m1, n2)
    return math.degrees(math.atan2(y, x))


def _segments(coords: Sequence[Point]) -> List[Tuple[int, int]]:
    """Split an index range at chain breaks, returning inclusive spans."""
    if not coords:
        return []
    spans: List[Tuple[int, int]] = []
    start = 0
    for index in range(1, len(coords)):
        if _distance(coords[index - 1], coords[index]) > _CHAIN_BREAK:
            spans.append((start, index - 1))
            start = index
    spans.append((start, len(coords) - 1))
    return spans


def _prune(assignment: List[str], code: str, minimum: int) -> None:
    """Erase runs of ``code`` shorter than ``minimum``, in place."""
    run_start = None
    for index in range(len(assignment) + 1):
        active = index < len(assignment) and assignment[index] == code
        if active and run_start is None:
            run_start = index
        elif not active and run_start is not None:
            if index - run_start < minimum:
                for position in range(run_start, index):
                    assignment[position] = "C"
            run_start = None


def assign_psea(coords: Sequence[Point]) -> List[str]:
    """Assign secondary structure from CA coordinates alone.

    A faithful-enough P-SEA: a residue window votes helix or strand when either
    its CA-CA distance profile or its angle/dihedral profile matches the
    reference geometry.  Helix wins ties because a helical distance profile is
    the more specific signal of the two.
    """
    count = len(coords)
    assignment = ["C"] * count
    if count < 5:
        return assignment

    for span_start, span_end in _segments(coords):
        length = span_end - span_start + 1
        if length < 5:
            continue

        helix_votes = [False] * count
        strand_votes = [False] * count

        for index in range(span_start, span_end + 1):
            d2 = _distance(coords[index], coords[index + 2]) if index + 2 <= span_end else None
            d3 = _distance(coords[index], coords[index + 3]) if index + 3 <= span_end else None
            d4 = _distance(coords[index], coords[index + 4]) if index + 4 <= span_end else None

            angle = None
            if span_start <= index - 1 and index + 1 <= span_end:
                angle = _angle(coords[index - 1], coords[index], coords[index + 1])

            dihedral = None
            if span_start <= index - 1 and index + 2 <= span_end:
                dihedral = _dihedral(
                    coords[index - 1], coords[index], coords[index + 1], coords[index + 2]
                )

            helix_by_distance = (
                _within(d2, _HELIX_D2) and _within(d3, _HELIX_D3) and _within(d4, _HELIX_D4)
            )
            helix_by_shape = _within(angle, _HELIX_ANGLE) and _within(dihedral, _HELIX_DIHEDRAL)
            if helix_by_distance or helix_by_shape:
                for offset in range(-1, 4):
                    position = index + offset
                    if span_start <= position <= span_end:
                        helix_votes[position] = True

            strand_by_distance = (
                _within(d2, _STRAND_D2) and _within(d3, _STRAND_D3) and _within(d4, _STRAND_D4)
            )
            strand_by_shape = _within(angle, _STRAND_ANGLE) and _within(
                dihedral, _STRAND_DIHEDRAL
            )
            if strand_by_distance or strand_by_shape:
                for offset in range(-1, 3):
                    position = index + offset
                    if span_start <= position <= span_end:
                        strand_votes[position] = True

        for index in range(span_start, span_end + 1):
            if helix_votes[index]:
                assignment[index] = "H"
            elif strand_votes[index]:
                assignment[index] = "E"

    _prune(assignment, "H", _MIN_HELIX)
    _prune(assignment, "E", _MIN_STRAND)
    return assignment


def _assign_biotite(coords: Sequence[Point]) -> Optional[List[str]]:
    """P-SEA via biotite, when it happens to be installed."""
    try:
        import numpy as np
        from biotite.structure import AtomArray, annotate_sse
    except Exception:
        return None
    try:
        atoms = AtomArray(len(coords))
        atoms.coord = np.array(coords, dtype=float)
        atoms.chain_id = np.full(len(coords), "A")
        atoms.res_id = np.arange(1, len(coords) + 1)
        atoms.res_name = np.full(len(coords), "GLY")
        atoms.atom_name = np.full(len(coords), "CA")
        atoms.element = np.full(len(coords), "C")
        codes = annotate_sse(atoms)
        if len(codes) != len(coords):
            return None
        return [{"a": "H", "b": "E"}.get(str(code), "C") for code in codes]
    except Exception:
        return None


def _assign_mkdssp(pdb_text: str, chain: str) -> Optional[List[str]]:
    """Reference DSSP, when the binary is on PATH.

    Needs full backbone atoms, so it is only offered the original file text
    rather than the CA-only view the other backends use.
    """
    binary = shutil.which("mkdssp") or shutil.which("dssp")
    if not binary or not pdb_text:
        return None
    try:
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "model.pdb"
            source.write_text(pdb_text, encoding="utf-8")
            result = subprocess.run(
                [binary, "--output-format", "dssp", str(source)],
                capture_output=True,
                text=True,
                timeout=60,
            )
            if result.returncode != 0 or not result.stdout:
                return None
            return _read_dssp_column(result.stdout, chain)
    except Exception:
        return None


def _read_dssp_column(text: str, chain: str) -> Optional[List[str]]:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.startswith("  #  RESIDUE"):
            body = lines[index + 1 :]
            break
    else:
        return None

    codes: List[str] = []
    for line in body:
        if len(line) < 17 or line[13] == "!":
            continue
        if chain and line[11] != chain:
            continue
        code = line[16]
        if code in "HGI":
            codes.append("H")
        elif code in "EB":
            codes.append("E")
        else:
            codes.append("C")
    return codes or None


def compute(
    coords: Sequence[Point],
    *,
    pdb_text: str = "",
    chain: str = "",
    prefer: str = "auto",
) -> Tuple[List[str], str]:
    """Derive secondary structure, returning the codes and their provenance.

    ``prefer`` pins a backend for testing; ``auto`` walks them in fidelity
    order.  The built-in path never fails, so this always returns something.
    """
    if prefer in {"auto", "mkdssp"}:
        codes = _assign_mkdssp(pdb_text, chain)
        if codes and len(codes) == len(coords):
            return codes, "mkdssp"
        if prefer == "mkdssp":
            return assign_psea(coords), "builtin (mkdssp unavailable)"

    if prefer in {"auto", "biotite"}:
        codes = _assign_biotite(coords)
        if codes:
            return codes, "biotite P-SEA"
        if prefer == "biotite":
            return assign_psea(coords), "builtin (biotite unavailable)"

    return assign_psea(coords), "built-in P-SEA"
