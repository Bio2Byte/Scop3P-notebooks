"""Layout engines.

Two ways of arranging the same elements, both emitting one schema so the
renderer and the annotation layer never learn which one ran:

    sheet        beta-sheet grouping preserved; the real topology diagram
    serpentine   strict N-to-C reading order wrapped across rows

Schema::

    {
      "mode": "sheet" | "serpentine",
      "elements":   [{id, kind, x, y, h, direction, path, start, stop}],
      "connectors": [{source, target, start, stop, path, lane, side}],
      "termini":    [{type, x, y, resnum}],
      "extents":    [min_x, min_y, max_x, max_y],
    }

Coordinates are unitless SVG user space; the renderer fits them to the
viewport.

Three defects in the original layout are addressed here.  Helix direction was
read off the last character of the element id, which is arbitrary; it now
follows N-to-C continuity.  Helices were positioned next to whichever strand
was nearest *in sequence*, which has nothing to do with space; they are now
placed by fitting the 3D coordinates onto the 2D sheet frame.  Connectors all
collapsed onto a shared midline and drew on top of each other; they are now
assigned separate lanes by interval colouring.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

# Sheet layout metrics.
_STRAND_PITCH = 74.0
_SHEET_GAP = 90.0
_ROW_Y = 260.0
_STAGGER = 8.0
_SHAFT = 12.0
_HEAD = 28.0
_MIN_HEIGHT = 72.0
_MAX_HEIGHT = 168.0

# How far off the strand row a helix sits, and how much the in-plane fit is
# allowed to shift it from that band.
_HELIX_CLEARANCE = 96.0
_HELIX_DRIFT = 140.0
_HELIX_STACK = 118.0

# Vertical gap between packed bands of sheets.
_BAND_GAP = 120.0

# Connector routing.
_LANE_SPACING = 16.0
_LANE_PADDING = 14.0

# How far a loop rises above (or drops below) the elements it joins, and how
# close two loops must be vertically before they are treated as colliding.
_LOOP_CLEARANCE = 20.0
_BAND_TOLERANCE = 70.0
_TERMINAL_STUB = 44.0

# Serpentine metrics.
_SERP_ROW_HEIGHT = 210.0
_SERP_MARGIN = 70.0
_SERP_GAP = 52.0


# --------------------------------------------------------------------------
# shared geometry
# --------------------------------------------------------------------------

def _element_height(length: int) -> float:
    return round(max(_MIN_HEIGHT, min(_MAX_HEIGHT, 34.0 + length * 7.0)), 2)


def _arrow_path(x: float, y: float, height: float, direction: int) -> List[float]:
    """Strand arrow, pointing along ``direction`` (+1 down, -1 up)."""
    if direction > 0:
        return [
            x - _SHAFT, y,
            x - _SHAFT, y + height - _HEAD,
            x - _HEAD, y + height - _HEAD,
            x, y + height,
            x + _HEAD, y + height - _HEAD,
            x + _SHAFT, y + height - _HEAD,
            x + _SHAFT, y,
        ]
    return [
        x + _SHAFT, y + height,
        x + _SHAFT, y + _HEAD,
        x + _HEAD, y + _HEAD,
        x, y,
        x - _HEAD, y + _HEAD,
        x - _SHAFT, y + _HEAD,
        x - _SHAFT, y + height,
    ]


def _endpoints(x: float, y: float, height: float, direction: int) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Return (N-terminal point, C-terminal point) for a placed element."""
    top = (x, y)
    bottom = (x, y + height)
    return (top, bottom) if direction > 0 else (bottom, top)


# --------------------------------------------------------------------------
# sheet grouping
# --------------------------------------------------------------------------

def _components(strand_ids: List[str], contacts: Sequence[Dict[str, Any]]) -> List[List[str]]:
    known = set(strand_ids)
    adjacency: Dict[str, List[str]] = {item: [] for item in strand_ids}
    for contact in contacts:
        if contact["source"] in known and contact["target"] in known:
            adjacency[contact["source"]].append(contact["target"])
            adjacency[contact["target"]].append(contact["source"])

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
    return groups


def _order_component(
    group: List[str],
    by_id: Dict[str, Dict[str, Any]],
    contacts: Sequence[Dict[str, Any]],
) -> List[str]:
    """Order strands across a sheet, starting from an edge strand."""
    if len(group) <= 1:
        return list(group)

    members = set(group)
    weights: Dict[Tuple[str, str], int] = {}
    degree: Dict[str, int] = {item: 0 for item in group}
    for contact in contacts:
        if contact["source"] in members and contact["target"] in members:
            key = tuple(sorted((contact["source"], contact["target"])))
            weights[key] = contact["count"]
            degree[contact["source"]] += 1
            degree[contact["target"]] += 1

    current = min(group, key=lambda item: (degree[item] != 1, by_id[item]["start"]))
    order = [current]
    remaining = members - {current}

    while remaining:
        best = max(
            remaining,
            key=lambda candidate: (
                weights.get(tuple(sorted((current, candidate))), 0),
                -abs(by_id[current]["start"] - by_id[candidate]["start"]),
            ),
        )
        order.append(best)
        remaining.discard(best)
        current = best
    return order


def _propagate_directions(
    order: List[str], contacts: Sequence[Dict[str, Any]]
) -> Dict[str, int]:
    """Set strand directions so paired strands agree with their sheet geometry."""
    members = set(order)
    neighbours: Dict[str, List[Tuple[str, str]]] = {item: [] for item in order}
    for contact in contacts:
        if contact["source"] in members and contact["target"] in members:
            neighbours[contact["source"]].append((contact["target"], contact["orientation"]))
            neighbours[contact["target"]].append((contact["source"], contact["orientation"]))

    directions: Dict[str, int] = {}
    for seed in order:
        if seed in directions:
            continue
        directions[seed] = 1
        queue = [seed]
        while queue:
            current = queue.pop(0)
            for neighbour, orientation in neighbours[current]:
                if neighbour in directions:
                    continue
                directions[neighbour] = (
                    directions[current] if orientation == "parallel" else -directions[current]
                )
                queue.append(neighbour)
    return directions


# --------------------------------------------------------------------------
# fitting helices onto the sheet frame
# --------------------------------------------------------------------------

from .projection import (
    assign_columns,
    choose_column_count,
    project_elements,
    relax_column,
    scale_vertical,
    uniform_scale,
)

Vector = Tuple[float, float, float]


def _sub(a: Vector, b: Vector) -> Vector:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot3(a: Vector, b: Vector) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross3(a: Vector, b: Vector) -> Vector:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _unit(a: Vector) -> Optional[Vector]:
    length = _dot3(a, a) ** 0.5
    if length < 1e-6:
        return None
    return (a[0] / length, a[1] / length, a[2] / length)


class SheetFrame:
    """An orthonormal frame built from the sheet itself.

    A general affine fit from 3D onto the diagram is ill-conditioned here,
    because strand centroids are close to coplanar and the coefficient normal
    to that plane is then unconstrained -- it throws helices arbitrarily far
    off. Instead the sheet supplies its own axes: ``along`` runs down the
    strands, ``across`` runs between neighbouring strands, and ``normal``
    completes the set. In-plane position is fitted (well-conditioned), while
    the out-of-plane component only chooses which side of the sheet a helix
    sits on, at a fixed clearance.
    """

    def __init__(self, origin: Vector, along: Vector, across: Vector, normal: Vector):
        self.origin = origin
        self.along = along
        self.across = across
        self.normal = normal
        self.scale_x = 1.0
        self.offset_x = 0.0
        self.scale_y = 1.0
        self.offset_y = 0.0

    def local(self, point: Vector) -> Tuple[float, float, float]:
        delta = _sub(point, self.origin)
        return (
            _dot3(delta, self.along),
            _dot3(delta, self.across),
            _dot3(delta, self.normal),
        )

    def to_screen(self, point: Vector) -> Tuple[float, float, float]:
        along, across, off_plane = self.local(point)
        return (
            self.scale_x * across + self.offset_x,
            self.scale_y * along + self.offset_y,
            off_plane,
        )


def _separate_helices(helix_positions: Dict[str, Tuple[float, float]]) -> None:
    """Push apart helices that the projection landed on top of each other."""
    entries = sorted(helix_positions.items(), key=lambda item: (item[1][1], item[1][0]))
    for index in range(len(entries)):
        element_id, (x, y) = entries[index]
        for _, (other_x, other_y) in entries[:index]:
            if abs(x - other_x) < _STRAND_PITCH and abs(y - other_y) < _HELIX_STACK:
                y = other_y + (_HELIX_STACK if y >= other_y else -_HELIX_STACK)
        entries[index] = (element_id, (x, y))
        helix_positions[element_id] = (x, y)


def _fit_linear(samples: List[Tuple[float, float]]) -> Tuple[float, float]:
    """One-dimensional least squares, returning (scale, offset)."""
    count = len(samples)
    if count < 2:
        return (1.0, samples[0][1] - samples[0][0] if samples else 0.0)
    mean_x = sum(item[0] for item in samples) / count
    mean_y = sum(item[1] for item in samples) / count
    variance = sum((item[0] - mean_x) ** 2 for item in samples)
    if variance < 1e-9:
        return (0.0, mean_y)
    covariance = sum((item[0] - mean_x) * (item[1] - mean_y) for item in samples)
    scale = covariance / variance
    return (scale, mean_y - scale * mean_x)


def _build_frame(
    strand_ids: List[str],
    by_id: Dict[str, Dict[str, Any]],
    positions: Dict[str, Tuple[float, float, int]],
) -> Optional[SheetFrame]:
    """Derive the sheet frame and calibrate it against the placed strands."""
    usable = [
        strand_id
        for strand_id in strand_ids
        if strand_id in positions and by_id[strand_id].get("centroid")
    ]
    if len(usable) < 3:
        return None

    # Average strand axis, with signs aligned so they do not cancel.
    reference = _unit(by_id[usable[0]]["axis"])
    if reference is None:
        return None
    accumulated = [0.0, 0.0, 0.0]
    for strand_id in usable:
        axis = _unit(by_id[strand_id]["axis"])
        if axis is None:
            continue
        sign = 1.0 if _dot3(axis, reference) >= 0 else -1.0
        for index in range(3):
            accumulated[index] += sign * axis[index]
    along = _unit(tuple(accumulated))
    if along is None:
        return None

    ordered = sorted(usable, key=lambda item: positions[item][0])
    first = by_id[ordered[0]]["centroid"]
    last = by_id[ordered[-1]]["centroid"]
    spread = _sub(last, first)

    # Remove any component running along the strands so the axes stay square.
    projection = _dot3(spread, along)
    across = _unit(
        (
            spread[0] - projection * along[0],
            spread[1] - projection * along[1],
            spread[2] - projection * along[2],
        )
    )
    if across is None:
        return None

    normal = _unit(_cross3(along, across))
    if normal is None:
        return None

    origin = by_id[ordered[0]]["centroid"]
    frame = SheetFrame(origin, along, across, normal)

    along_samples: List[Tuple[float, float]] = []
    across_samples: List[Tuple[float, float]] = []
    for strand_id in usable:
        local_along, local_across, _ = frame.local(by_id[strand_id]["centroid"])
        screen_x, screen_y, _ = positions[strand_id]
        across_samples.append((local_across, screen_x))
        along_samples.append((local_along, screen_y))

    frame.scale_x, frame.offset_x = _fit_linear(across_samples)
    frame.scale_y, frame.offset_y = _fit_linear(along_samples)
    return frame


# --------------------------------------------------------------------------
# connector routing
# --------------------------------------------------------------------------

def _assign_lanes(spans: List[Tuple[int, float, float, float]]) -> Dict[int, int]:
    """Greedy interval colouring, per band, so connectors do not overlap.

    Two connectors only conflict when they overlap horizontally *and* sit at a
    similar height. Colouring globally forces separate lanes on connectors that
    are nowhere near each other, which pushes them further and further out.
    """
    lanes: Dict[int, int] = {}
    by_band: Dict[int, List[Tuple[int, float, float]]] = {}

    for key, start, end, reference in spans:
        band = int(round(reference / _BAND_TOLERANCE))
        by_band.setdefault(band, []).append((key, start, end))

    for band_items in by_band.values():
        occupied: List[float] = []
        for key, start, end in sorted(band_items, key=lambda item: min(item[1], item[2])):
            low, high = min(start, end), max(start, end)
            for lane, last_end in enumerate(occupied):
                if low >= last_end + _LANE_PADDING:
                    occupied[lane] = high
                    lanes[key] = lane
                    break
            else:
                occupied.append(high)
                lanes[key] = len(occupied) - 1
    return lanes


def _route_connectors(
    placed: List[Dict[str, Any]],
    bounds: Tuple[float, float],
) -> List[Dict[str, Any]]:
    """Route loops locally, hugging the elements they join.

    The previous version sent every connector out to a lane beyond the whole
    diagram's bounding box, so a loop between two neighbours travelled the full
    height of the figure to get there and back. Here a connector rises only just
    above (or drops just below) the two elements it actually joins.

    Which side to use is not a free choice: the C-terminal end of one element
    and the N-terminal end of the next each sit at a definite end, decided by
    the direction the element runs. When both are at the top the loop belongs
    above; when both are at the bottom it belongs below; and when they disagree
    the path has to cross between them anyway, so it goes straight through the
    gap rather than around the outside.
    """
    pending: List[Dict[str, Any]] = []
    spans: Dict[str, List[Tuple[int, float, float, float]]] = {"above": [], "below": []}

    for index in range(len(placed) - 1):
        left = placed[index]
        right = placed[index + 1]
        exit_point = left["c_point"]
        entry_point = right["n_point"]

        # Where each endpoint physically sits on its element.
        exit_side = "below" if left["direction"] > 0 else "above"
        entry_side = "above" if right["direction"] > 0 else "below"

        record = {
            "source": left["id"],
            "target": right["id"],
            "start": left["stop"] + 1,
            "stop": right["start"] - 1,
            "from": exit_point,
            "to": entry_point,
            "key": index,
        }

        if exit_side == entry_side:
            record["side"] = exit_side
            reference = (
                min(exit_point[1], entry_point[1])
                if exit_side == "above"
                else max(exit_point[1], entry_point[1])
            )
            record["reference"] = reference
            spans[exit_side].append((index, exit_point[0], entry_point[0], reference))
        else:
            # Ends face opposite ways, so the loop crosses the gap directly.
            record["side"] = "through"
            record["reference"] = (exit_point[1] + entry_point[1]) / 2.0

        pending.append(record)

    lanes = {
        "above": _assign_lanes(spans["above"]),
        "below": _assign_lanes(spans["below"]),
    }

    connectors: List[Dict[str, Any]] = []
    for record in pending:
        key = record["key"]
        side = record["side"]
        start_point = record["from"]
        end_point = record["to"]

        if side == "through":
            lane = 0
            mid_y = record["reference"]
            path = [
                start_point[0], start_point[1],
                start_point[0], mid_y,
                end_point[0], mid_y,
                end_point[0], end_point[1],
            ]
        else:
            lane = lanes[side].get(key, 0)
            offset = _LOOP_CLEARANCE + lane * _LANE_SPACING
            lane_y = (
                record["reference"] - offset
                if side == "above"
                else record["reference"] + offset
            )
            path = [
                start_point[0], start_point[1],
                start_point[0], lane_y,
                end_point[0], lane_y,
                end_point[0], end_point[1],
            ]

        has_residues = record["start"] <= record["stop"]
        connectors.append(
            {
                "source": record["source"],
                "target": record["target"],
                "start": record["start"] if has_residues else -1,
                "stop": record["stop"] if has_residues else -1,
                "path": [round(value, 2) for value in path],
                "lane": lane,
                "side": side,
            }
        )
    return connectors


def _extents(
    placed: List[Dict[str, Any]], connectors: List[Dict[str, Any]], termini: List[Dict[str, Any]]
) -> List[float]:
    xs: List[float] = []
    ys: List[float] = []
    for element in placed:
        xs.extend(element["path"][0::2])
        ys.extend(element["path"][1::2])
    for connector in connectors:
        xs.extend(connector["path"][0::2])
        ys.extend(connector["path"][1::2])
    for terminus in termini:
        xs.append(terminus["x"])
        ys.append(terminus["y"])
    if not xs or not ys:
        return [0.0, 0.0, 100.0, 100.0]
    return [min(xs), min(ys), max(xs), max(ys)]


def _finish(
    placed: List[Dict[str, Any]],
    residues: Sequence[Any],
    mode: str,
) -> Dict[str, Any]:
    """Shared tail: route connectors, add termini, compute extents."""
    if not placed:
        return {
            "mode": mode,
            "elements": [],
            "connectors": [],
            "termini": [],
            "extents": [0.0, 0.0, 100.0, 100.0],
        }

    tops = [element["y"] for element in placed]
    bottoms = [element["y"] + element["h"] for element in placed]
    connectors = _route_connectors(placed, (min(tops), max(bottoms)))

    first = placed[0]
    last = placed[-1]
    first_seq = residues[0].seq if residues else first["start"]
    last_seq = residues[-1].seq if residues else last["stop"]

    n_point = first["n_point"]
    c_point = last["c_point"]
    n_offset = -_TERMINAL_STUB if first["direction"] > 0 else _TERMINAL_STUB
    c_offset = _TERMINAL_STUB if last["direction"] > 0 else -_TERMINAL_STUB

    termini = [
        {
            "type": "N",
            "x": round(n_point[0], 2),
            "y": round(n_point[1] + n_offset, 2),
            "anchor": [round(n_point[0], 2), round(n_point[1], 2)],
            "resnum": first_seq,
        },
        {
            "type": "C",
            "x": round(c_point[0], 2),
            "y": round(c_point[1] + c_offset, 2),
            "anchor": [round(c_point[0], 2), round(c_point[1], 2)],
            "resnum": last_seq,
        },
    ]

    elements_out = [
        {
            "id": element["id"],
            "kind": element["kind"],
            "x": round(element["x"], 2),
            "y": round(element["y"], 2),
            "h": round(element["h"], 2),
            "direction": element["direction"],
            "start": element["start"],
            "stop": element["stop"],
            "path": [round(value, 2) for value in element["path"]],
            "n_point": [round(v, 2) for v in element["n_point"]],
            "c_point": [round(v, 2) for v in element["c_point"]],
        }
        for element in placed
    ]

    return {
        "mode": mode,
        "elements": elements_out,
        "connectors": connectors,
        "termini": termini,
        "extents": [round(value, 2) for value in _extents(elements_out, connectors, termini)],
    }


def _place(
    element: Dict[str, Any], x: float, y: float, direction: int
) -> Dict[str, Any]:
    height = _element_height(element["length"])
    n_point, c_point = _endpoints(x, y, height, direction)
    return {
        "id": element["id"],
        "kind": element["type"],
        "x": x,
        "y": y,
        "h": height,
        "direction": direction,
        "start": element["start"],
        "stop": element["stop"],
        "path": _arrow_path(x, y, height, direction),
        "n_point": n_point,
        "c_point": c_point,
    }


# --------------------------------------------------------------------------
# sheet layout
# --------------------------------------------------------------------------

def _assign_helices(
    helices: List[Dict[str, Any]],
    groups: List[List[str]],
    by_id: Dict[str, Dict[str, Any]],
) -> Dict[int, List[str]]:
    """Attach each helix to the sheet it packs against in space.

    A helix belongs beside the sheet it actually touches, which is a question
    about 3D distance, not about which strand happens to be nearest in the
    sequence.
    """
    assignment: Dict[int, List[str]] = {index: [] for index in range(len(groups))}
    if not groups:
        return assignment

    centres: List[Optional[Vector]] = []
    for group in groups:
        points = [by_id[item].get("centroid") for item in group]
        points = [point for point in points if point]
        if not points:
            centres.append(None)
            continue
        count = float(len(points))
        centres.append((
            sum(p[0] for p in points) / count,
            sum(p[1] for p in points) / count,
            sum(p[2] for p in points) / count,
        ))

    for helix in helices:
        centroid = helix.get("centroid")
        best_index = 0
        if centroid:
            best_distance = None
            for index, centre in enumerate(centres):
                if centre is None:
                    continue
                delta = _sub(centroid, centre)
                distance = _dot3(delta, delta)
                if best_distance is None or distance < best_distance:
                    best_distance = distance
                    best_index = index
        else:
            # No geometry: fall back to the sheet nearest in sequence.
            best_index = min(
                range(len(groups)),
                key=lambda i: min(abs(by_id[m]["start"] - helix["start"]) for m in groups[i]),
            )
        assignment[best_index].append(helix["id"])

    return assignment


def _build_block(
    order: List[str],
    helix_ids: List[str],
    by_id: Dict[str, Dict[str, Any]],
    contacts: Sequence[Dict[str, Any]],
) -> Dict[str, Any]:
    """Lay out one sheet and its helices in local coordinates.

    Local origin is the left edge of the strand row, with the row centred on
    y = 0, so blocks can be translated freely when they are packed.
    """
    directions = _propagate_directions(order, contacts)
    local: Dict[str, Tuple[float, float, int]] = {}

    for index, strand_id in enumerate(order):
        local[strand_id] = (
            index * _STRAND_PITCH,
            (index % 2) * _STAGGER,
            directions.get(strand_id, 1),
        )

    frame = _build_frame(order, by_id, local)

    row_span = max(1, len(order) - 1) * _STRAND_PITCH
    for position, helix_id in enumerate(helix_ids):
        helix = by_id[helix_id]
        centroid = helix.get("centroid")
        clearance = _HELIX_CLEARANCE + _element_height(helix["length"]) / 2.0

        if frame and centroid:
            x, in_plane_y, off_plane = frame.to_screen(centroid)
            y = clearance if off_plane >= 0 else -clearance
            y += max(-_HELIX_DRIFT, min(_HELIX_DRIFT, in_plane_y)) * 0.12
            x = max(-_STRAND_PITCH, min(row_span + _STRAND_PITCH, x))
        else:
            # Spread evenly above the row when there is nothing to fit against.
            slots = max(1, len(helix_ids))
            x = (row_span * position) / slots if slots > 1 else row_span / 2.0
            y = -clearance
        local[helix_id] = (x, y, 0)

    helix_positions = {item: (local[item][0], local[item][1]) for item in helix_ids}
    _separate_helices(helix_positions)
    for helix_id, (x, y) in helix_positions.items():
        local[helix_id] = (x, y, 0)

    xs: List[float] = []
    ys: List[float] = []
    for item, (x, y, _) in local.items():
        half = _element_height(by_id[item]["length"]) / 2.0
        xs.extend([x - _HEAD, x + _HEAD])
        ys.extend([y - half, y + half])

    return {
        "local": local,
        "min_x": min(xs), "max_x": max(xs),
        "min_y": min(ys), "max_y": max(ys),
        "width": max(xs) - min(xs),
        "height": max(ys) - min(ys),
        "order": min(by_id[item]["start"] for item in local),
    }


def _shelf_pack(
    blocks: List[Dict[str, Any]],
    target_aspect: float,
    gap_x: float,
    gap_y: float,
) -> List[Tuple[float, float]]:
    """Pack blocks into bands, choosing the band count that best fits the frame.

    Laying every sheet on one row is what makes a multi-sheet protein render as
    a long thin strip that has to be scaled down to nothing to fit the panel.
    Sweeping the band width and scoring the resulting aspect keeps the diagram
    close to the shape of the space it is drawn into.
    """
    if not blocks:
        return []

    widths = [block["width"] for block in blocks]
    total_width = sum(widths) + gap_x * max(0, len(blocks) - 1)

    best_offsets: List[Tuple[float, float]] = []
    best_score: Optional[float] = None

    for band_count in range(1, len(blocks) + 1):
        limit = total_width / band_count
        bands: List[List[int]] = []
        current: List[int] = []
        current_width = 0.0

        for index, width in enumerate(widths):
            addition = width if not current else gap_x + width
            if current and current_width + addition > limit:
                bands.append(current)
                current = [index]
                current_width = width
            else:
                current.append(index)
                current_width += addition
        if current:
            bands.append(current)

        band_widths = [
            sum(widths[i] for i in band) + gap_x * max(0, len(band) - 1) for band in bands
        ]
        band_heights = [max(blocks[i]["height"] for i in band) for band in bands]
        canvas_width = max(band_widths)
        canvas_height = sum(band_heights) + gap_y * max(0, len(bands) - 1)
        if canvas_height <= 0:
            continue

        aspect = canvas_width / canvas_height
        score = abs(math.log(aspect / target_aspect))
        if best_score is not None and score >= best_score:
            continue

        offsets: List[Tuple[float, float]] = [(0.0, 0.0)] * len(blocks)
        cursor_y = 0.0
        for band, band_width, band_height in zip(bands, band_widths, band_heights):
            # Centre each band so the diagram stays visually balanced.
            cursor_x = (canvas_width - band_width) / 2.0
            for index in band:
                offsets[index] = (cursor_x, cursor_y + band_height / 2.0)
                cursor_x += blocks[index]["width"] + gap_x
            cursor_y += band_height + gap_y

        best_score = score
        best_offsets = offsets

    return best_offsets


def layout_sheet(
    elements: List[Dict[str, Any]],
    residues: Sequence[Any],
    contacts: Sequence[Dict[str, Any]],
    target_aspect: float = 1.3,
) -> Dict[str, Any]:
    """Arrange strands into sheets, hang helices off them, and pack to frame."""
    if not elements:
        return _finish([], residues, "sheet")

    by_id = {element["id"]: element for element in elements}
    strands = [element for element in elements if element["type"] == "strand"]
    helices = [element for element in elements if element["type"] == "helix"]

    groups = _components([item["id"] for item in strands], contacts)
    groups.sort(key=lambda group: min(by_id[item]["start"] for item in group))

    if not groups:
        # All-helix chain: no sheet to organise around, so fall back to the
        # sequence layout, which is the honest representation of that topology.
        return layout_serpentine(elements, residues, contacts, target_aspect=target_aspect)

    helix_map = _assign_helices(helices, groups, by_id)

    blocks = [
        _build_block(_order_component(group, by_id, contacts), helix_map[index], by_id, contacts)
        for index, group in enumerate(groups)
    ]

    offsets = _shelf_pack(blocks, target_aspect, _SHEET_GAP, _BAND_GAP)

    positions: Dict[str, Tuple[float, float, int]] = {}
    for block, (offset_x, offset_y) in zip(blocks, offsets):
        shift_x = offset_x - block["min_x"]
        for item, (x, y, direction) in block["local"].items():
            positions[item] = (x + shift_x, y + offset_y, direction)

    # Direction by N-to-C continuity: pick the orientation that leaves the
    # shortest hop from the previous element's C-terminal exit point.
    placed: List[Dict[str, Any]] = []
    previous_exit: Optional[Tuple[float, float]] = None

    for element in elements:
        centre_x, centre_y, direction = positions.get(element["id"], (0.0, 0.0, 1))
        y = centre_y - _element_height(element["length"]) / 2.0
        candidates = [direction] if element["type"] == "strand" and direction else [1, -1]

        best = None
        for candidate in candidates:
            trial = _place(element, centre_x, y, candidate)
            if previous_exit is None:
                cost = 0.0
            else:
                entry = trial["n_point"]
                cost = abs(entry[0] - previous_exit[0]) + abs(entry[1] - previous_exit[1])
            if best is None or cost < best[0]:
                best = (cost, trial)

        assert best is not None
        placed.append(best[1])
        previous_exit = best[1]["c_point"]

    return _finish(placed, residues, "sheet")


# --------------------------------------------------------------------------
# serpentine layout
# --------------------------------------------------------------------------

def layout_serpentine(
    elements: List[Dict[str, Any]],
    residues: Sequence[Any],
    contacts: Sequence[Dict[str, Any]] = (),
    per_row: Optional[int] = None,
    target_aspect: float = 1.3,
) -> Dict[str, Any]:
    """Lay elements out strictly N to C, wrapping in boustrophedon rows.

    Sheet grouping is sacrificed, but connectors stay short and local and the
    horizontal axis is monotone in residue number, which leaves clean room for
    per-residue annotation.

    The row width is chosen to fit the frame rather than from a fixed rule: a
    guess like sqrt(count) produces a wide strip for one protein and a tall
    column for the next, and both waste most of the panel once scaled to fit.
    """
    if not elements:
        return _finish([], residues, "serpentine")

    count = len(elements)
    column_pitch = _STRAND_PITCH + _SERP_GAP

    if per_row is None:
        best_choice, best_score = count, None
        for candidate in range(1, count + 1):
            rows = math.ceil(count / candidate)
            width = max(1, min(candidate, count) - 1) * column_pitch + _HEAD * 2
            height = rows * _SERP_ROW_HEIGHT
            score = abs(math.log((width / height) / target_aspect))
            if best_score is None or score < best_score:
                best_score, best_choice = score, candidate
        per_row = best_choice

    per_row = max(1, per_row)

    placed: List[Dict[str, Any]] = []
    for index, element in enumerate(elements):
        row = index // per_row
        column = index % per_row
        if row % 2 == 1:
            column = per_row - 1 - column

        x = _SERP_MARGIN + column * column_pitch
        height = _element_height(element["length"])
        y = _SERP_MARGIN + row * _SERP_ROW_HEIGHT + (_MAX_HEIGHT - height) / 2.0

        direction = 1 if (index % 2 == 0) else -1
        placed.append(_place(element, x, y, direction))

    return _finish(placed, residues, "serpentine")


# --------------------------------------------------------------------------
# projection layout
# --------------------------------------------------------------------------

_PROJ_PITCH = 68.0
_PROJ_GAP = 30.0
_PROJ_SCALE = 9.0


def _project_with_columns(
    elements: List[Dict[str, Any]],
    projected: Dict[str, Tuple[float, float, int]],
    column_count: int,
    target_aspect: float,
) -> Tuple[List[Dict[str, Any]], float]:
    """Place elements for a fixed column count, and report the aspect achieved."""
    order = [element["id"] for element in elements]
    by_id = {element["id"]: element for element in elements}

    columns = assign_columns(projected, order, column_count)
    heights = {item: _element_height(by_id[item]["length"]) for item in order}

    max_height = (column_count * _PROJ_PITCH) / max(0.2, target_aspect)
    scale = uniform_scale(
        projected, order, column_count, _PROJ_PITCH, max_height
    )
    depths = scale_vertical(projected, order, scale)

    grouped: Dict[int, List[Tuple[str, float, float]]] = {}
    for element_id in order:
        grouped.setdefault(columns[element_id], []).append(
            (element_id, depths[element_id], heights[element_id])
        )

    centres: Dict[str, float] = {}
    for entries in grouped.values():
        entries.sort(key=lambda item: item[1])
        centres.update(relax_column(entries, _PROJ_GAP))

    offset = _SERP_MARGIN - min(centres[item] - heights[item] / 2.0 for item in order)

    placed: List[Dict[str, Any]] = []
    for element in elements:
        element_id = element["id"]
        x = _SERP_MARGIN + columns[element_id] * _PROJ_PITCH
        y = centres[element_id] - heights[element_id] / 2.0 + offset
        placed.append(_place(element, x, y, projected[element_id][2]))

    width = max(item["x"] for item in placed) - min(item["x"] for item in placed) + _HEAD * 2
    top = min(item["y"] for item in placed)
    bottom = max(item["y"] + item["h"] for item in placed)
    aspect = width / max(1e-6, bottom - top)
    return placed, aspect


def layout_projection(
    elements: List[Dict[str, Any]],
    residues: Sequence[Any],
    contacts: Sequence[Dict[str, Any]] = (),
    target_aspect: float = 1.3,
) -> Dict[str, Any]:
    """Arrange elements by projecting the real 3D packing onto a plane.

    Elements finish near each other on the page when they are near each other in
    space, which is what makes the PDBe diagram legible: sheets come out
    adjacent because they are adjacent, helices sit beside what they pack
    against, and connectors stay short because the things they join are close.
    """
    if not elements:
        return _finish([], residues, "projection")

    projected = project_elements(elements)
    if projected is None:
        # No geometry to project. Sequence order is the honest fallback.
        return layout_serpentine(elements, residues, contacts, target_aspect=target_aspect)

    # The projection does not spread evenly, so a column count derived from
    # averages overshoots whenever several elements land in one column. Sweeping
    # and scoring the achieved shape is more reliable than predicting it.
    best_placed: Optional[List[Dict[str, Any]]] = None
    best_score: Optional[float] = None

    for column_count in range(1, len(elements) + 1):
        placed, aspect = _project_with_columns(
            elements, projected, column_count, target_aspect
        )
        score = abs(math.log(max(1e-6, aspect) / target_aspect))
        if best_score is None or score < best_score:
            best_score, best_placed = score, placed

    return _finish(best_placed or [], residues, "projection")


# --------------------------------------------------------------------------
# spatial arrangement
# --------------------------------------------------------------------------

_SEG_ROW_HEIGHT = 232.0
_SEG_COL_PITCH = 118.0


def layout_spatial(
    elements: List[Dict[str, Any]],
    residues: Sequence[Any],
    contacts: Sequence[Dict[str, Any]] = (),
    target_aspect: float = 1.3,
) -> Dict[str, Any]:
    """Segmented rows, every row reading left to right.

    Unlike the serpentine, rows are not reversed on alternate passes. Snaking
    keeps the connector between rows short, but it costs the reader the one
    thing a topology diagram should never make them work for: every other row
    runs backwards, so following the chain means constantly re-checking which
    way to read. Here every row starts at the left, and the wrap between rows is
    paid for with a single visible sweep instead.

    Element directions alternate so the C-terminal end of one lands opposite the
    N-terminal end of the next, which lets the loop between them stay a short
    hop rather than a climb around the element.
    """
    if not elements:
        return _finish([], residues, "spatial")

    count = len(elements)

    per_row, best_score = count, None
    for candidate in range(1, count + 1):
        rows = math.ceil(count / candidate)
        width = max(1, min(candidate, count)) * _SEG_COL_PITCH
        height = rows * _SEG_ROW_HEIGHT
        score = abs(math.log((width / height) / target_aspect))
        if best_score is None or score < best_score:
            best_score, per_row = score, candidate
    per_row = max(1, per_row)

    placed: List[Dict[str, Any]] = []
    for index, element in enumerate(elements):
        row = index // per_row
        column = index % per_row

        x = _SERP_MARGIN + column * _SEG_COL_PITCH
        height = _element_height(element["length"])
        y = _SERP_MARGIN + row * _SEG_ROW_HEIGHT + (_MAX_HEIGHT - height) / 2.0
        direction = 1 if (index % 2 == 0) else -1
        placed.append(_place(element, x, y, direction))

    result = _finish(placed, residues, "spatial")

    # Row bands, so a reader can see where one segment ends and the next begins.
    segments: List[Dict[str, Any]] = []
    for row in range(math.ceil(count / per_row)):
        members = placed[row * per_row : (row + 1) * per_row]
        if not members:
            continue
        segments.append({
            "index": row + 1,
            "label": f"segment {row + 1}",
            "x": min(item["x"] for item in members) - _HEAD - 18,
            "width": (max(item["x"] for item in members)
                      - min(item["x"] for item in members)) + _HEAD * 2 + 36,
            "y": min(item["y"] for item in members) - 34,
            "height": (max(item["y"] + item["h"] for item in members)
                       - min(item["y"] for item in members)) + 56,
        })
    result["segments"] = segments
    return result


LAYOUTS = {
    "spatial": layout_spatial,
    "projection": layout_projection,
    "sheet": layout_sheet,
    "serpentine": layout_serpentine,
}


def build_layout(
    mode: str,
    elements: List[Dict[str, Any]],
    residues: Sequence[Any],
    contacts: Sequence[Dict[str, Any]],
    target_aspect: float = 1.3,
) -> Dict[str, Any]:
    """Build one layout. ``target_aspect`` is width/height of the frame it will
    be drawn into, so the packing can aim at the shape of the actual panel."""
    engine = LAYOUTS.get(mode, layout_sheet)
    return engine(elements, residues, contacts, target_aspect=target_aspect)
