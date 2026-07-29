"""Projection layout.

Sheet grouping and sequence order both throw away most of what a topology
diagram is supposed to show. PDBe's arrangement reads well because elements sit
near each other on the page when they sit near each other in space: strands of
one sheet end up adjacent because they *are* adjacent, helices land beside
whatever they pack against, and connectors stay short as a consequence rather
than by routing effort.

So the arrangement is derived rather than composed:

    1. fit the plane the chain occupies (PCA over element centroids)
    2. project centroids onto it
    3. orient so element axes run vertically, the way they are drawn
    4. snap to columns, which is what gives the tidy PDBe look
    5. relax overlaps within each column

Directions come from the projected N-to-C axis, so an arrow points the way the
element actually runs in space.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

Vector = Tuple[float, float, float]


def _sub(a: Vector, b: Vector) -> Vector:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Vector, b: Vector) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a: Vector) -> float:
    return math.sqrt(_dot(a, a))


def _jacobi_eigen(matrix: List[List[float]], sweeps: int = 60):
    """Eigen decomposition of a symmetric 3x3, by cyclic Jacobi rotation.

    Written out rather than pulled from numpy so the package keeps working on a
    bare Voila deployment with no scientific stack installed.
    """
    size = 3
    a = [row[:] for row in matrix]
    vectors = [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]

    for _ in range(sweeps):
        off_diagonal = sum(a[i][j] ** 2 for i in range(size) for j in range(size) if i != j)
        if off_diagonal < 1e-12:
            break
        for p in range(size - 1):
            for q in range(p + 1, size):
                if abs(a[p][q]) < 1e-14:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
                sign = 1.0 if theta >= 0 else -1.0
                t = sign / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c

                for k in range(size):
                    akp = a[k][p]
                    akq = a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(size):
                    apk = a[p][k]
                    aqk = a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(size):
                    vkp = vectors[k][p]
                    vkq = vectors[k][q]
                    vectors[k][p] = c * vkp - s * vkq
                    vectors[k][q] = s * vkp + c * vkq

    values = [a[i][i] for i in range(size)]
    columns = [tuple(vectors[row][col] for row in range(size)) for col in range(size)]
    order = sorted(range(size), key=lambda i: -values[i])
    return [values[i] for i in order], [columns[i] for i in order]


def principal_axes(points: Sequence[Vector]):
    """Return (centre, axes) with axes ordered by descending variance."""
    count = len(points)
    if count == 0:
        return (0.0, 0.0, 0.0), [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]

    centre = (
        sum(p[0] for p in points) / count,
        sum(p[1] for p in points) / count,
        sum(p[2] for p in points) / count,
    )
    covariance = [[0.0] * 3 for _ in range(3)]
    for point in points:
        delta = _sub(point, centre)
        for i in range(3):
            for j in range(3):
                covariance[i][j] += delta[i] * delta[j]
    for i in range(3):
        for j in range(3):
            covariance[i][j] /= max(1, count)

    _, axes = _jacobi_eigen(covariance)
    return centre, axes


def project_elements(
    elements: List[Dict[str, Any]],
) -> Optional[Dict[str, Tuple[float, float, int]]]:
    """Project elements onto their best-fit plane.

    Returns id -> (u, v, direction) in structure units, or None when the
    elements carry no geometry to project.
    """
    points = [element.get("centroid") for element in elements]
    if any(point is None for point in points) or len(points) < 3:
        return None

    centre, axes = principal_axes(points)

    # The chain's dominant direction becomes the page's vertical, because
    # elements are drawn as vertical bars. The second axis becomes horizontal,
    # so neighbouring strands of a sheet land in neighbouring columns.
    vertical, horizontal = axes[0], axes[1]

    # Strand axes should agree with the page vertical; if they mostly oppose it,
    # flip so arrows read downward as the chain progresses.
    agreement = 0.0
    for element in elements:
        axis = element.get("axis")
        if axis:
            agreement += _dot(axis, vertical)
    if agreement < 0:
        vertical = (-vertical[0], -vertical[1], -vertical[2])

    projected: Dict[str, Tuple[float, float, int]] = {}
    for element in elements:
        delta = _sub(element["centroid"], centre)
        u = _dot(delta, horizontal)
        v = _dot(delta, vertical)

        axis = element.get("axis") or (0.0, 0.0, 0.0)
        along = _dot(axis, vertical)
        direction = 1 if along >= 0 else -1
        projected[element["id"]] = (u, v, direction)

    return projected


def choose_column_count(
    element_count: int,
    mean_height: float,
    pitch: float,
    gap: float,
    target_aspect: float,
) -> int:
    """How many columns to spread across, so the result matches the frame.

    Solved rather than guessed. With ``n`` elements over ``c`` columns the
    canvas is about ``c * pitch`` wide and ``(n / c) * (mean_height + gap)``
    tall, so setting that ratio equal to the target gives
    ``c = sqrt(target * n * (mean_height + gap) / pitch)``.
    """
    if element_count <= 1:
        return 1
    numerator = target_aspect * element_count * (mean_height + gap)
    columns = math.sqrt(max(1.0, numerator / max(1e-6, pitch)))
    return max(1, min(element_count, int(round(columns))))


def assign_columns(
    projected: Dict[str, Tuple[float, float, int]],
    order: List[str],
    column_count: int,
) -> Dict[str, int]:
    """Snap projected positions onto integer columns.

    Free projection leaves elements at arbitrary offsets, which reads as
    scattered. Columns are what make the PDBe diagram look deliberate, and they
    also give connectors clean vertical runs to follow.
    """
    if not projected:
        return {}
    if column_count <= 1:
        return {item: 0 for item in order}

    values = [projected[item][0] for item in order]
    low, high = min(values), max(values)
    span = high - low

    if span < 1e-6:
        return {item: index % column_count for index, item in enumerate(order)}

    scale = (column_count - 1) / span
    return {
        item: int(round((projected[item][0] - low) * scale)) for item in order
    }


def uniform_scale(
    projected: Dict[str, Tuple[float, float, int]],
    order: List[str],
    column_count: int,
    pitch: float,
    max_height: float,
) -> float:
    """One scale factor for both axes.

    Stretching the vertical to fill a target height destroys the thing that
    makes a projection worth drawing: it is no longer a projection, just a
    scatter with the proportions thrown away. A single factor keeps it
    faithful.

    The factor is bounded by both axes and the smaller wins, because deriving
    it from the horizontal alone explodes whenever elements happen to be
    clustered in that direction: a small horizontal spread forces a large scale,
    which then throws the vertical off the canvas.
    """
    us = [projected[item][0] for item in order]
    vs = [projected[item][1] for item in order]
    u_span = max(us) - min(us)
    v_span = max(vs) - min(vs)

    candidates = []
    if u_span > 1e-6 and column_count > 1:
        candidates.append(((column_count - 1) * pitch) / u_span)
    if v_span > 1e-6:
        candidates.append(max_height / v_span)

    return min(candidates) if candidates else pitch


def scale_vertical(
    projected: Dict[str, Tuple[float, float, int]],
    order: List[str],
    scale: float,
) -> Dict[str, float]:
    """Map projected depth onto the page at the given scale."""
    values = [projected[item][1] for item in order]
    low = min(values)
    return {item: (projected[item][1] - low) * scale for item in order}


def relax_column(
    entries: List[Tuple[str, float, float]],
    gap: float,
) -> Dict[str, float]:
    """Space elements within one column so none overlap, preserving their order.

    A single pass from the top is enough because the entries are pre-sorted:
    each element is pushed just far enough to clear the one above it.
    """
    placed: Dict[str, float] = {}
    previous_bottom: Optional[float] = None

    for element_id, centre, height in entries:
        top = centre - height / 2.0
        if previous_bottom is not None and top < previous_bottom + gap:
            centre = previous_bottom + gap + height / 2.0
        placed[element_id] = centre
        previous_bottom = centre + height / 2.0

    return placed
