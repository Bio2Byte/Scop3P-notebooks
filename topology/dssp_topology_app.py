from __future__ import annotations

import json
import math
import re
import urllib.error
import urllib.request
import uuid
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


DEFAULT_DSSP_PATH = (
    r"C:\Users\drpat\OneDrive - UGent\Desktop\PD\Scop3P_update2025"
    r"\Scop3P_AF_DSSP\DSSP_flatfiles\P07949.dssp"
)

FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d+|\d+)(?:[Ee][-+]?\d+)?")

SS_NAMES = {
    "H": "Alpha helix",
    "G": "3-10 helix",
    "I": "Pi helix",
    "E": "Beta strand",
    "B": "Beta bridge",
    "T": "Turn",
    "S": "Bend",
    "C": "Coil",
}

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
}


def _parse_int(value: str) -> Optional[int]:
    value = value.strip()
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _ss_kind(ss: str) -> str:
    if ss in {"H", "G", "I"}:
        return "helix"
    if ss in {"E", "B"}:
        return "strand"
    if ss == "T":
        return "turn"
    if ss == "S":
        return "bend"
    return "coil"


def _clean_header_value(line: str) -> str:
    return line.rstrip(" .").strip()


def _extract_metadata(lines: Iterable[str]) -> Dict[str, Any]:
    metadata: Dict[str, Any] = {}
    for line in lines:
        if line.startswith("HEADER"):
            metadata["header"] = _clean_header_value(line[10:])
        elif line.startswith("COMPND"):
            molecule_match = re.search(r"MOLECULE:\s*([^;]+)", line)
            chain_match = re.search(r"CHAIN:\s*([^;]+)", line)
            if molecule_match:
                metadata["molecule"] = molecule_match.group(1).strip()
            if chain_match:
                metadata["declared_chain"] = chain_match.group(1).strip()
        elif "TOTAL NUMBER OF RESIDUES" in line:
            total = _parse_int(line[:6])
            if total is not None:
                metadata["reported_residue_count"] = total
    return metadata


def _parse_residue_line(line: str) -> Optional[Dict[str, Any]]:
    if len(line) < 17:
        return None

    dssp_index = _parse_int(line[0:5])
    if dssp_index is None:
        return None

    residue_number = line[5:10].strip()
    if not residue_number:
        return None

    insertion_code = line[10:11].strip()
    chain = line[11:12].strip() or "."
    aa = line[13:14].strip() or "X"
    if aa == "!":
        return None

    ss = line[16:17].strip() or "C"
    accessibility = _parse_int(line[34:38]) if len(line) >= 38 else None
    bp1 = _parse_int(line[25:29]) if len(line) >= 29 else None
    bp2 = _parse_int(line[29:33]) if len(line) >= 33 else None

    coords = [float(match.group(0)) for match in FLOAT_RE.finditer(line)]
    x, y, z = (coords[-3:] if len(coords) >= 3 else (None, None, None))

    residue_id = f"{chain}:{residue_number}{insertion_code}"
    resseq_int = _parse_int(residue_number)

    return {
        "dssp_index": dssp_index,
        "residue_number": residue_number,
        "resseq_int": resseq_int,
        "insertion_code": insertion_code,
        "chain": chain,
        "aa": aa,
        "ss": ss,
        "ss_name": SS_NAMES.get(ss, "Coil"),
        "ss_kind": _ss_kind(ss),
        "accessibility": accessibility,
        "bp1": bp1,
        "bp2": bp2,
        "x": x,
        "y": y,
        "z": z,
        "residue_id": residue_id,
    }


def parse_dssp(text: str) -> Dict[str, Any]:
    """Parse a DSSP flatfile into residues, secondary-structure elements, and links."""
    lines = text.splitlines()
    header_index = next(
        (i for i, line in enumerate(lines) if line.lstrip().startswith("#  RESIDUE")),
        None,
    )
    if header_index is None:
        raise ValueError("Could not find the DSSP residue table header.")

    residues = [
        residue
        for line in lines[header_index + 1 :]
        if (residue := _parse_residue_line(line)) is not None
    ]
    if not residues:
        raise ValueError("No DSSP residue rows were parsed from the file.")

    elements = _build_elements(residues)
    links = _build_beta_links(residues, elements)
    chains = sorted({residue["chain"] for residue in residues})
    metadata = _extract_metadata(lines[:header_index])

    return {
        "metadata": metadata,
        "residues": residues,
        "elements": elements,
        "links": links,
        "stats": {
            "residue_count": len(residues),
            "element_count": len(elements),
            "helix_count": sum(1 for element in elements if element["type"] == "helix"),
            "strand_count": sum(1 for element in elements if element["type"] == "strand"),
            "chain_count": len(chains),
            "chains": chains,
            "beta_link_count": len(links),
        },
    }


def _can_extend_element(element: Dict[str, Any], residue: Dict[str, Any]) -> bool:
    last = element["residues"][-1]
    return (
        residue["chain"] == element["chain"]
        and residue["ss"] == element["ss_code"]
        and residue["dssp_index"] == last["dssp_index"] + 1
    )


def _build_elements(residues: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    elements: List[Dict[str, Any]] = []
    active: Optional[Dict[str, Any]] = None
    counters = {"helix": 0, "strand": 0}

    def close_active() -> None:
        nonlocal active
        if active is None:
            return
        residues_in_element = active.pop("residues")
        counters[active["type"]] += 1
        prefix = "H" if active["type"] == "helix" else "S"
        element_id = f"{prefix}{counters[active['type']]}"
        first = residues_in_element[0]
        last = residues_in_element[-1]
        active.update(
            {
                "id": element_id,
                "label": element_id,
                "start_dssp": first["dssp_index"],
                "end_dssp": last["dssp_index"],
                "start_residue": first["residue_number"] + first["insertion_code"],
                "end_residue": last["residue_number"] + last["insertion_code"],
                "start_residue_id": first["residue_id"],
                "end_residue_id": last["residue_id"],
                "length": len(residues_in_element),
                "sequence": "".join(residue["aa"] for residue in residues_in_element),
                "residue_ids": [residue["residue_id"] for residue in residues_in_element],
                "residue_indices": [
                    residue["dssp_index"] for residue in residues_in_element
                ],
                "residue_numbers": [
                    residue["residue_number"] + residue["insertion_code"]
                    for residue in residues_in_element
                ],
                "accessibility_mean": _mean(
                    residue["accessibility"] for residue in residues_in_element
                ),
            }
        )
        elements.append(active)
        active = None

    for residue in residues:
        if residue["ss_kind"] not in {"helix", "strand"}:
            close_active()
            continue

        if active is not None and _can_extend_element(active, residue):
            active["residues"].append(residue)
            continue

        close_active()
        active = {
            "type": residue["ss_kind"],
            "ss_code": residue["ss"],
            "ss_name": residue["ss_name"],
            "chain": residue["chain"],
            "residues": [residue],
        }

    close_active()
    return elements


def _mean(values: Iterable[Optional[int]]) -> Optional[float]:
    clean = [value for value in values if value is not None]
    if not clean:
        return None
    return round(sum(clean) / len(clean), 1)


def _build_beta_links(
    residues: List[Dict[str, Any]], elements: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    element_by_dssp: Dict[int, str] = {}
    element_by_resseq: Dict[Tuple[str, int], str] = {}

    for element in elements:
        if element["type"] != "strand":
            continue
        for dssp_index in element["residue_indices"]:
            element_by_dssp[dssp_index] = element["id"]
        for residue_number in element["residue_numbers"]:
            resseq = _parse_int(residue_number)
            if resseq is not None:
                element_by_resseq[(element["chain"], resseq)] = element["id"]

    link_counts: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for residue in residues:
        source = element_by_dssp.get(residue["dssp_index"])
        if not source:
            continue
        for partner in (residue["bp1"], residue["bp2"]):
            if not partner:
                continue
            target = element_by_dssp.get(partner) or element_by_resseq.get(
                (residue["chain"], partner)
            )
            if not target or target == source:
                continue
            source_id, target_id = sorted((source, target))
            key = (source_id, target_id)
            if key not in link_counts:
                link_counts[key] = {
                    "source": source_id,
                    "target": target_id,
                    "count": 0,
                    "examples": [],
                }
            link_counts[key]["count"] += 1
            if len(link_counts[key]["examples"]) < 6:
                link_counts[key]["examples"].append(residue["residue_id"])

    return sorted(
        link_counts.values(),
        key=lambda link: (-link["count"], link["source"], link["target"]),
    )


def topology_from_dssp(
    text: str,
    name: str = "uploaded.dssp",
    afdb_accession: Optional[str] = None,
) -> Dict[str, Any]:
    topology = parse_dssp(text)
    topology["name"] = name
    if afdb_accession:
        topology["afdb_accession"] = afdb_accession.strip()
    return topology


def _guess_afdb_accession(name: str) -> str:
    stem = Path(name).stem.upper()
    if stem.startswith("AF-"):
        parts = stem.split("-")
        if len(parts) >= 2:
            return parts[1]
    match = re.search(r"\b[A-Z][A-Z0-9]{5,9}\b", stem)
    return match.group(0) if match else ""


def _fetch_json(url: str) -> Any:
    request = urllib.request.Request(url, headers={"User-Agent": "Scop3P-topology/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def _fetch_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "Scop3P-topology/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8", errors="replace")


def fetch_alphafold_cif(accession: str) -> Tuple[Dict[str, Any], str]:
    accession = accession.strip()
    if accession.upper().startswith("AF-"):
        accession = _guess_afdb_accession(accession) or accession
    if not accession:
        raise ValueError("Enter a UniProt or AlphaFold DB accession first.")
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    payload = _fetch_json(api_url)
    record = payload[0] if isinstance(payload, list) and payload else payload
    if not record:
        raise ValueError(f"No AlphaFold DB prediction found for {accession}.")
    cif_url = record.get("cifUrl")
    if not cif_url:
        raise ValueError(f"AlphaFold DB did not return a CIF URL for {accession}.")
    return record, _fetch_text(cif_url)


def _cif_tokens(text: str) -> List[str]:
    tokens: List[str] = []
    i = 0
    n = len(text)
    at_line_start = True
    while i < n:
        char = text[i]
        if char in " \t\r\n":
            at_line_start = char in "\r\n"
            i += 1
            continue
        if char == "#":
            while i < n and text[i] not in "\r\n":
                i += 1
            at_line_start = True
            continue
        if char == ";" and at_line_start:
            i += 1
            start = i
            end = text.find("\n;", i)
            if end == -1:
                tokens.append(text[start:].strip())
                break
            tokens.append(text[start:end].strip())
            i = end + 2
            at_line_start = False
            continue
        if char in {"'", '"'}:
            quote = char
            i += 1
            start = i
            value_parts: List[str] = []
            while i < n:
                if text[i] == quote and (i + 1 == n or text[i + 1] in " \t\r\n#"):
                    value_parts.append(text[start:i])
                    i += 1
                    break
                i += 1
            else:
                value_parts.append(text[start:i])
            tokens.append("".join(value_parts))
            at_line_start = False
            continue
        start = i
        while i < n and text[i] not in " \t\r\n#":
            i += 1
        tokens.append(text[start:i])
        at_line_start = False
    return tokens


def _parse_cif_loops(text: str) -> Dict[str, List[Dict[str, str]]]:
    tokens = _cif_tokens(text)
    loops: Dict[str, List[Dict[str, str]]] = {}
    i = 0
    while i < len(tokens):
        if tokens[i].lower() != "loop_":
            i += 1
            continue
        i += 1
        headers: List[str] = []
        while i < len(tokens) and tokens[i].startswith("_"):
            headers.append(tokens[i])
            i += 1
        if not headers:
            continue
        category = headers[0].split(".", 1)[0]
        width = len(headers)
        while i < len(tokens):
            token = tokens[i]
            if token.lower() == "loop_" or token.lower().startswith("data_"):
                break
            if token.startswith("_"):
                break
            row = tokens[i : i + width]
            if len(row) < width:
                break
            loops.setdefault(category, []).append(dict(zip(headers, row)))
            i += width
    return loops


def _to_int(value: Any) -> Optional[int]:
    if value in {None, "", ".", "?"}:
        return None
    try:
        return int(str(value))
    except ValueError:
        return None


def _to_float(value: Any) -> Optional[float]:
    if value in {None, "", ".", "?"}:
        return None
    try:
        return float(str(value))
    except ValueError:
        return None


def _conf_kind(conf_type: str) -> Optional[str]:
    value = conf_type.upper()
    if value.startswith("HELX"):
        return "helix"
    if value.startswith("STRN"):
        return "strand"
    return None


def _distance(a: Dict[str, Any], b: Dict[str, Any]) -> float:
    return math.sqrt(
        (a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2
    )


def _dot(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _vector_for_element(element: Dict[str, Any], residue_by_seq: Dict[int, Dict[str, Any]]) -> Tuple[float, float, float]:
    start = residue_by_seq.get(element["start"])
    stop = residue_by_seq.get(element["stop"])
    if not start or not stop:
        return (0.0, 1.0, 0.0)
    return (stop["x"] - start["x"], stop["y"] - start["y"], stop["z"] - start["z"])


def _strand_contacts(
    elements: List[Dict[str, Any]],
    residue_by_seq: Dict[int, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    strands = [element for element in elements if element["type"] == "strand"]
    contacts: List[Dict[str, Any]] = []
    for i, left in enumerate(strands):
        left_residues = [
            residue_by_seq[seq]
            for seq in range(left["start"], left["stop"] + 1)
            if seq in residue_by_seq
        ]
        for right in strands[i + 1 :]:
            right_residues = [
                residue_by_seq[seq]
                for seq in range(right["start"], right["stop"] + 1)
                if seq in residue_by_seq
            ]
            examples: List[str] = []
            count = 0
            min_distance = None
            for residue_a in left_residues:
                for residue_b in right_residues:
                    if abs(residue_a["seq"] - residue_b["seq"]) < 3:
                        continue
                    distance = _distance(residue_a, residue_b)
                    if min_distance is None or distance < min_distance:
                        min_distance = distance
                    if distance <= 7.0:
                        count += 1
                        if len(examples) < 6:
                            examples.append(f"{residue_a['chain']}:{residue_a['seq']}-{residue_b['seq']}")
            if count:
                orientation = (
                    "parallel"
                    if _dot(
                        _vector_for_element(left, residue_by_seq),
                        _vector_for_element(right, residue_by_seq),
                    )
                    >= 0
                    else "antiparallel"
                )
                contacts.append(
                    {
                        "source": left["id"],
                        "target": right["id"],
                        "count": count,
                        "orientation": orientation,
                        "min_distance": round(min_distance or 0.0, 2),
                        "examples": examples,
                    }
                )
    return sorted(contacts, key=lambda item: (-item["count"], item["source"], item["target"]))


def _arrow_path(x: float, y: float, height: float, direction: int) -> List[float]:
    shaft = 12.0
    head = 30.0
    if direction > 0:
        return [
            x - shaft,
            y,
            x - shaft,
            y + height - head,
            x - head,
            y + height - head,
            x,
            y + height,
            x + head,
            y + height - head,
            x + shaft,
            y + height - head,
            x + shaft,
            y,
        ]
    return [
        x + shaft,
        y + height,
        x + shaft,
        y + head,
        x + head,
        y + head,
        x,
        y,
        x - head,
        y + head,
        x - shaft,
        y + head,
        x - shaft,
        y + height,
    ]


def _endpoint(element: Dict[str, Any], which: str) -> Tuple[float, float]:
    x = element["layout_x"]
    y = element["layout_y"]
    height = element["layout_h"]
    direction = element["direction"]
    if which == "start":
        return (x, y) if direction > 0 else (x, y + height)
    return (x, y + height) if direction > 0 else (x, y)


def _coil_path(start: Tuple[float, float], end: Tuple[float, float]) -> List[float]:
    mid_y = (start[1] + end[1]) / 2
    return [start[0], start[1], start[0], mid_y, end[0], mid_y, end[0], end[1]]


def _sheet_components(strands: List[Dict[str, Any]], contacts: List[Dict[str, Any]]) -> List[List[str]]:
    strand_ids = {strand["id"] for strand in strands}
    adjacency: Dict[str, List[str]] = {strand_id: [] for strand_id in strand_ids}
    for contact in contacts:
        if contact["source"] in strand_ids and contact["target"] in strand_ids:
            adjacency[contact["source"]].append(contact["target"])
            adjacency[contact["target"]].append(contact["source"])
    seen: set[str] = set()
    components: List[List[str]] = []
    for strand in strands:
        if strand["id"] in seen:
            continue
        stack = [strand["id"]]
        seen.add(strand["id"])
        component: List[str] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adjacency[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(component)
    return components


def _order_sheet_component(
    component: List[str],
    elements_by_id: Dict[str, Dict[str, Any]],
    contacts: List[Dict[str, Any]],
) -> List[str]:
    if len(component) <= 1:
        return component
    weights: Dict[Tuple[str, str], int] = {}
    degree: Dict[str, int] = {item: 0 for item in component}
    for contact in contacts:
        if contact["source"] in degree and contact["target"] in degree:
            key = tuple(sorted((contact["source"], contact["target"])))
            weights[key] = contact["count"]
            degree[contact["source"]] += 1
            degree[contact["target"]] += 1
    current = min(component, key=lambda item: (degree[item] != 1, elements_by_id[item]["start"]))
    order = [current]
    remaining = set(component) - {current}
    while remaining:
        candidates = [
            (
                weights.get(tuple(sorted((current, candidate))), 0),
                -abs(elements_by_id[current]["start"] - elements_by_id[candidate]["start"]),
                candidate,
            )
            for candidate in remaining
        ]
        candidates.sort(reverse=True)
        current = candidates[0][2]
        order.append(current)
        remaining.remove(current)
    return order


def _build_pdbe_topology_paths(
    elements: List[Dict[str, Any]],
    residues: List[Dict[str, Any]],
    contacts: List[Dict[str, Any]],
) -> Dict[str, Any]:
    residue_by_seq = {residue["seq"]: residue for residue in residues}
    elements_by_id = {element["id"]: element for element in elements}
    strands = [element for element in elements if element["type"] == "strand"]
    helices = [element for element in elements if element["type"] == "helix"]
    components = _sheet_components(strands, contacts)
    components.sort(key=lambda comp: min(elements_by_id[item]["start"] for item in comp))

    strand_positions: Dict[str, Tuple[float, float, int]] = {}
    y_base = 220.0
    x_cursor = 70.0
    for component in components:
        order = _order_sheet_component(component, elements_by_id, contacts)
        width = max(1, len(order))
        direction_by_id: Dict[str, int] = {order[0]: 1}
        for index, strand_id in enumerate(order):
            if strand_id not in direction_by_id:
                direction_by_id[strand_id] = 1 if index % 2 == 0 else -1
            strand_positions[strand_id] = (
                x_cursor + index * 72.0,
                y_base + (index % 2) * 8.0,
                direction_by_id[strand_id],
            )
            for contact in contacts:
                if strand_id not in {contact["source"], contact["target"]}:
                    continue
                neighbor = contact["target"] if contact["source"] == strand_id else contact["source"]
                if neighbor not in component or neighbor in direction_by_id:
                    continue
                direction_by_id[neighbor] = (
                    direction_by_id[strand_id]
                    if contact["orientation"] == "parallel"
                    else -direction_by_id[strand_id]
                )
        x_cursor += width * 72.0 + 76.0

    sequence_min = min((element["start"] for element in elements), default=1)
    sequence_max = max((element["stop"] for element in elements), default=1)
    sequence_span = max(1, sequence_max - sequence_min)
    for element in elements:
        if element["type"] == "strand":
            x, y, direction = strand_positions.get(
                element["id"],
                (70.0 + (element["start"] - sequence_min) * 3.0, y_base, 1),
            )
        else:
            midpoint = (element["start"] + element["stop"]) / 2
            x = 70.0 + ((midpoint - sequence_min) / sequence_span) * max(260.0, x_cursor - 120.0)
            nearest_strand = min(
                strands,
                key=lambda strand: abs(((strand["start"] + strand["stop"]) / 2) - midpoint),
                default=None,
            )
            if nearest_strand and nearest_strand["id"] in strand_positions:
                _, strand_y, _ = strand_positions[nearest_strand["id"]]
                y = strand_y - 150.0 if element["start"] < nearest_strand["start"] else strand_y + 160.0
            else:
                y = 80.0 + (len(helices) % 3) * 95.0
            direction = 1 if element["id"].endswith(("1", "3", "5", "7", "9")) else -1
        element["layout_x"] = round(x, 3)
        element["layout_y"] = round(y, 3)
        element["layout_h"] = round(max(70.0, min(150.0, 34.0 + element["length"] * 7.5)), 3)
        element["direction"] = direction

    strands_data: List[Dict[str, Any]] = []
    helices_data: List[Dict[str, Any]] = []
    for element in elements:
        if element["type"] == "strand":
            strands_data.append(
                {
                    "start": element["start"],
                    "stop": element["stop"],
                    "path": _arrow_path(
                        element["layout_x"],
                        element["layout_y"],
                        element["layout_h"],
                        element["direction"],
                    ),
                }
            )
        else:
            start = _endpoint(element, "start")
            end = _endpoint(element, "end")
            helices_data.append(
                {
                    "start": element["start"],
                    "stop": element["stop"],
                    "path": [start[0] - 18.0, start[1], end[0] + 18.0, end[1]],
                    "majoraxis": 18.0,
                    "minoraxis": 9.0,
                }
            )

    coils_data: List[Dict[str, Any]] = []
    first_residue = residues[0]["seq"]
    last_residue = residues[-1]["seq"]
    if elements:
        first = elements[0]
        n_start = (_endpoint(first, "start")[0], _endpoint(first, "start")[1] + (-46 if first["direction"] > 0 else 46))
        n_has_residues = first["start"] > first_residue
        coils_data.append(
            {
                "start": first_residue if n_has_residues else -1,
                "stop": first["start"] - 1 if n_has_residues else -1,
                "path": _coil_path(n_start, _endpoint(first, "start")),
            }
        )
        for left, right in zip(elements, elements[1:]):
            start = left["stop"] + 1
            stop = right["start"] - 1
            coils_data.append(
                {
                    "start": start if start <= stop else -1,
                    "stop": stop if start <= stop else -1,
                    "path": _coil_path(_endpoint(left, "end"), _endpoint(right, "start")),
                }
            )
        last = elements[-1]
        c_end = (_endpoint(last, "end")[0], _endpoint(last, "end")[1] + (46 if last["direction"] > 0 else -46))
        c_has_residues = last["stop"] < last_residue
        coils_data.append(
            {
                "start": last["stop"] + 1 if c_has_residues else -1,
                "stop": last_residue if c_has_residues else -1,
                "path": _coil_path(_endpoint(last, "end"), c_end),
            }
        )

    all_xy = []
    for item in strands_data + helices_data + coils_data:
        all_xy.extend(item["path"])
    xs = all_xy[0::2] or [0.0, 100.0]
    ys = all_xy[1::2] or [0.0, 100.0]
    return {
        "strands": strands_data,
        "coils": coils_data,
        "terms": [
            {"type": "N", "resnum": str(first_residue), "start": -1, "stop": -1, "path": []},
            {"type": "C", "resnum": str(last_residue), "start": -1, "stop": -1, "path": []},
        ],
        "helices": helices_data,
        "extents": [min(xs), min(ys), max(xs), max(ys)],
    }


def topology_from_alphafold_cif(
    text: str,
    name: str = "alphafold.cif",
    afdb_metadata: Optional[Dict[str, Any]] = None,
    chain_id: Optional[str] = None,
) -> Dict[str, Any]:
    loops = _parse_cif_loops(text)
    atom_rows = loops.get("_atom_site", [])
    qa_rows = loops.get("_ma_qa_metric_local", [])
    plddt_by_seq: Dict[int, float] = {}
    for row in qa_rows:
        seq = _to_int(row.get("_ma_qa_metric_local.label_seq_id"))
        value = _to_float(row.get("_ma_qa_metric_local.metric_value"))
        if seq is not None and value is not None:
            plddt_by_seq[seq] = value

    residues_by_key: Dict[Tuple[str, int], Dict[str, Any]] = {}
    for row in atom_rows:
        atom_name = row.get("_atom_site.label_atom_id") or row.get("_atom_site.auth_atom_id")
        if atom_name != "CA":
            continue
        chain = (
            row.get("_atom_site.label_asym_id")
            or row.get("_atom_site.auth_asym_id")
            or "A"
        )
        seq = _to_int(row.get("_atom_site.label_seq_id") or row.get("_atom_site.auth_seq_id"))
        if seq is None:
            continue
        x = _to_float(row.get("_atom_site.Cartn_x"))
        y = _to_float(row.get("_atom_site.Cartn_y"))
        z = _to_float(row.get("_atom_site.Cartn_z"))
        if x is None or y is None or z is None:
            continue
        comp = row.get("_atom_site.label_comp_id") or row.get("_atom_site.auth_comp_id") or "UNK"
        b_factor = _to_float(row.get("_atom_site.B_iso_or_equiv"))
        residues_by_key[(chain, seq)] = {
            "seq": seq,
            "chain": chain,
            "aa": AA3_TO_1.get(comp.upper(), "X"),
            "comp_id": comp,
            "x": x,
            "y": y,
            "z": z,
            "plddt": plddt_by_seq.get(seq, b_factor),
            "residue_id": f"{chain}:{seq}",
        }

    if not residues_by_key:
        raise ValueError("No CA atom coordinates were found in the mmCIF file.")

    chain_counts: Dict[str, int] = {}
    for chain, _ in residues_by_key:
        chain_counts[chain] = chain_counts.get(chain, 0) + 1
    selected_chain = chain_id or max(chain_counts, key=chain_counts.get)
    residues = sorted(
        [residue for (chain, _), residue in residues_by_key.items() if chain == selected_chain],
        key=lambda residue: residue["seq"],
    )
    residue_by_seq = {residue["seq"]: residue for residue in residues}

    conf_rows = loops.get("_struct_conf", [])
    elements: List[Dict[str, Any]] = []
    counters = {"helix": 0, "strand": 0}
    for row in conf_rows:
        kind = _conf_kind(row.get("_struct_conf.conf_type_id", ""))
        if kind not in {"helix", "strand"}:
            continue
        beg_chain = row.get("_struct_conf.beg_label_asym_id") or row.get("_struct_conf.beg_auth_asym_id")
        end_chain = row.get("_struct_conf.end_label_asym_id") or row.get("_struct_conf.end_auth_asym_id")
        if beg_chain != selected_chain or end_chain != selected_chain:
            continue
        start = _to_int(row.get("_struct_conf.beg_label_seq_id") or row.get("_struct_conf.beg_auth_seq_id"))
        stop = _to_int(row.get("_struct_conf.end_label_seq_id") or row.get("_struct_conf.end_auth_seq_id"))
        if start is None or stop is None or stop < start:
            continue
        if not any(seq in residue_by_seq for seq in range(start, stop + 1)):
            continue
        counters[kind] += 1
        element_id = ("H" if kind == "helix" else "S") + str(counters[kind])
        sequence = "".join(residue_by_seq[seq]["aa"] for seq in range(start, stop + 1) if seq in residue_by_seq)
        elements.append(
            {
                "id": element_id,
                "label": element_id,
                "type": kind,
                "chain": selected_chain,
                "start": start,
                "stop": stop,
                "start_residue": str(start),
                "end_residue": str(stop),
                "start_residue_id": f"{selected_chain}:{start}",
                "end_residue_id": f"{selected_chain}:{stop}",
                "length": stop - start + 1,
                "sequence": sequence,
                "ss_code": "H" if kind == "helix" else "E",
                "ss_name": "Alpha helix" if kind == "helix" else "Beta strand",
                "residue_indices": list(range(start, stop + 1)),
                "residue_ids": [f"{selected_chain}:{seq}" for seq in range(start, stop + 1)],
                "residue_numbers": [str(seq) for seq in range(start, stop + 1)],
            }
        )
    elements.sort(key=lambda element: (element["start"], element["stop"]))
    if not elements:
        raise ValueError("No helix or strand records were found in _struct_conf.")

    contacts = _strand_contacts(elements, residue_by_seq)
    pdbe_topology = _build_pdbe_topology_paths(elements, residues, contacts)
    entry_id = (
        str(afdb_metadata.get("entryId") or afdb_metadata.get("modelEntityId"))
        if afdb_metadata
        else Path(name).stem
    )
    entry_key = entry_id.lower()
    sequence = (
        str(afdb_metadata.get("sequence") or afdb_metadata.get("uniprotSequence"))
        if afdb_metadata
        else "".join(residue["aa"] for residue in residues)
    )
    uniprot = str(afdb_metadata.get("uniprotAccession") or "") if afdb_metadata else ""
    metadata = {
        "source": "AlphaFold DB mmCIF" if afdb_metadata else "Uploaded mmCIF",
        "entry_id": entry_id,
        "chain": selected_chain,
        "uniprot": uniprot,
        "cif_url": afdb_metadata.get("cifUrl") if afdb_metadata else None,
        "model_created": afdb_metadata.get("modelCreatedDate") if afdb_metadata else None,
        "latest_version": afdb_metadata.get("latestVersion") if afdb_metadata else None,
    }

    mapping_records: Dict[str, Any] = {}
    if uniprot:
        mapping_records["UniProt"] = {
            uniprot: {
                "mappings": [
                    {
                        "entity_id": "1",
                        "chain_id": selected_chain,
                        "start": {"residue_number": residues[0]["seq"]},
                        "end": {"residue_number": residues[-1]["seq"]},
                    }
                ]
            }
        }
    low_confidence_residues = []
    for residue in residues:
        plddt = residue.get("plddt")
        if plddt is None or plddt >= 70:
            continue
        label = "pLDDT <50" if plddt < 50 else "pLDDT 50-70"
        low_confidence_residues.append(
            {"residue_number": residue["seq"], "outlier_types": [label]}
        )

    pdbe_api_data = [
        {entry_key: [{"entity_id": "1", "sequence": sequence}]},
        {entry_key: mapping_records},
        {entry_key: {"1": {selected_chain: pdbe_topology}}},
        {
            entry_key: {
                "molecules": [
                    {
                        "entity_id": "1",
                        "chains": [
                            {
                                "chain_id": selected_chain,
                                "models": [{"model_id": 1, "residues": low_confidence_residues}],
                            }
                        ],
                    }
                ]
            }
        },
        {
            entry_key: {
                "molecules": [
                    {
                        "entity_id": "1",
                        "chains": [
                            {
                                "chain_id": selected_chain,
                                "observed": [
                                    {
                                        "start": {"residue_number": residues[0]["seq"]},
                                        "end": {"residue_number": residues[-1]["seq"]},
                                    }
                                ],
                            }
                        ],
                    }
                ]
            }
        },
    ]

    return {
        "name": name,
        "metadata": metadata,
        "residues": residues,
        "elements": elements,
        "links": contacts,
        "pdbe_api_data": pdbe_api_data,
        "pdbe_entry_id": entry_key,
        "pdbe_entity_id": "1",
        "pdbe_chain_id": selected_chain,
        "afdb_accession": uniprot,
        "cif_url": metadata["cif_url"],
        "stats": {
            "residue_count": len(residues),
            "element_count": len(elements),
            "helix_count": sum(1 for element in elements if element["type"] == "helix"),
            "strand_count": sum(1 for element in elements if element["type"] == "strand"),
            "chain_count": 1,
            "chains": [selected_chain],
            "beta_link_count": len(contacts),
        },
    }


def _legacy_topology_html(topology: Dict[str, Any]) -> str:
    root_id = "dssp-topology-" + uuid.uuid4().hex
    data_id = root_id + "-data"
    title = topology.get("name") or "DSSP topology"
    metadata = topology.get("metadata", {})
    stats = topology["stats"]
    subtitle = metadata.get("molecule") or metadata.get("header") or "Uploaded DSSP"
    json_blob = json.dumps(topology, separators=(",", ":")).replace("</", "<\\/")

    css = """
#__ROOT_ID__ {
  --ink: #1f2933;
  --muted: #5d6978;
  --line: #d4dce6;
  --panel: #f7f9fc;
  --helix: #d94b62;
  --helix-alt: #f59f43;
  --helix-pi: #8b5cf6;
  --strand: #2878b8;
  --strand-alt: #21a1a1;
  --strand-link: #5b8bd8;
  color: var(--ink);
  font-family: Inter, "Segoe UI", Arial, sans-serif;
  line-height: 1.35;
}
#__ROOT_ID__ .topology-shell {
  border: 1px solid var(--line);
  border-radius: 8px;
  overflow: hidden;
  background: white;
}
#__ROOT_ID__ .topology-header {
  display: grid;
  grid-template-columns: minmax(260px, 1fr) auto;
  gap: 12px;
  align-items: center;
  padding: 14px 16px;
  border-bottom: 1px solid var(--line);
  background: var(--panel);
}
#__ROOT_ID__ .title {
  font-size: 17px;
  font-weight: 700;
}
#__ROOT_ID__ .subtitle {
  color: var(--muted);
  font-size: 13px;
  margin-top: 2px;
}
#__ROOT_ID__ .stats {
  color: var(--muted);
  font-size: 12px;
  margin-top: 6px;
}
#__ROOT_ID__ .controls {
  display: flex;
  flex-wrap: wrap;
  justify-content: flex-end;
  gap: 6px;
}
#__ROOT_ID__ button,
#__ROOT_ID__ input[type="text"] {
  border: 1px solid #c6d0dc;
  border-radius: 6px;
  background: #ffffff;
  color: var(--ink);
  font: inherit;
  font-size: 12px;
  min-height: 30px;
}
#__ROOT_ID__ button {
  cursor: pointer;
  padding: 4px 9px;
}
#__ROOT_ID__ button:hover {
  border-color: #8ca1b8;
  background: #eef4fb;
}
#__ROOT_ID__ input[type="text"] {
  width: 150px;
  padding: 4px 8px;
}
#__ROOT_ID__ label.toggle {
  align-items: center;
  color: var(--muted);
  display: inline-flex;
  font-size: 12px;
  gap: 4px;
  min-height: 30px;
  white-space: nowrap;
}
#__ROOT_ID__ .viewer-wrap {
  height: 640px;
  position: relative;
  background: linear-gradient(#ffffff, #fbfcfe);
}
#__ROOT_ID__ svg {
  cursor: grab;
  display: block;
  height: 100%;
  width: 100%;
}
#__ROOT_ID__ svg.is-dragging {
  cursor: grabbing;
}
#__ROOT_ID__ .connector {
  fill: none;
  stroke: #a7b1bf;
  stroke-dasharray: 4 5;
  stroke-linecap: round;
  stroke-width: 1.5;
}
#__ROOT_ID__ .beta-link {
  fill: none;
  opacity: 0.34;
  stroke: var(--strand-link);
  stroke-linecap: round;
}
#__ROOT_ID__ .beta-link.link-selected {
  opacity: 0.9;
  stroke: #234f9a;
}
#__ROOT_ID__ .sse {
  cursor: pointer;
}
#__ROOT_ID__ .sse-shape {
  filter: drop-shadow(0 1px 1px rgba(31, 41, 51, 0.16));
  stroke: rgba(31, 41, 51, 0.28);
  stroke-width: 1;
}
#__ROOT_ID__ .sse:hover .sse-shape,
#__ROOT_ID__ .sse.selected .sse-shape {
  stroke: #111827;
  stroke-width: 2;
}
#__ROOT_ID__ .sse.selected .selection-ring {
  opacity: 1;
}
#__ROOT_ID__ .selection-ring {
  fill: none;
  opacity: 0;
  pointer-events: none;
  stroke: #111827;
  stroke-dasharray: 5 4;
  stroke-width: 1.5;
}
#__ROOT_ID__ .element-label {
  fill: #17202a;
  font-size: 12px;
  font-weight: 700;
  pointer-events: none;
  text-anchor: middle;
}
#__ROOT_ID__ .residue-label {
  fill: #687586;
  font-size: 10px;
  pointer-events: none;
  text-anchor: middle;
}
#__ROOT_ID__ .helix-stripe {
  fill: none;
  opacity: 0.45;
  pointer-events: none;
  stroke: #ffffff;
  stroke-linecap: round;
  stroke-width: 2;
}
#__ROOT_ID__ .tooltip {
  background: #111827;
  border-radius: 6px;
  color: white;
  display: none;
  font-size: 12px;
  max-width: 280px;
  padding: 8px 10px;
  pointer-events: none;
  position: absolute;
  z-index: 4;
}
#__ROOT_ID__ .details {
  border-top: 1px solid var(--line);
  color: var(--muted);
  font-size: 13px;
  min-height: 44px;
  padding: 12px 16px;
}
#__ROOT_ID__ .details strong {
  color: var(--ink);
}
#__ROOT_ID__ .details code {
  background: #eef2f7;
  border-radius: 4px;
  color: #344154;
  padding: 1px 4px;
}
@media (max-width: 800px) {
  #__ROOT_ID__ .topology-header {
    grid-template-columns: 1fr;
  }
  #__ROOT_ID__ .controls {
    justify-content: flex-start;
  }
  #__ROOT_ID__ .viewer-wrap {
    height: 540px;
  }
}
""".replace(
        "__ROOT_ID__", root_id
    )

    script = r"""
(function () {
  const root = document.getElementById("__ROOT_ID__");
  const data = JSON.parse(document.getElementById("__DATA_ID__").textContent);
  const svg = root.querySelector("svg");
  const viewport = root.querySelector("[data-role='viewport']");
  const connectorsLayer = root.querySelector("[data-role='connectors']");
  const linksLayer = root.querySelector("[data-role='links']");
  const elementsLayer = root.querySelector("[data-role='elements']");
  const tooltip = root.querySelector("[data-role='tooltip']");
  const details = root.querySelector("[data-role='details']");
  const linkToggle = root.querySelector("[data-role='toggle-links']");
  const labelToggle = root.querySelector("[data-role='toggle-labels']");
  const residueInput = root.querySelector("[data-role='residue-input']");
  const elementById = new Map(data.elements.map((element) => [element.id, element]));
  const elementByResidue = new Map();
  let selectedId = null;
  let transform = { x: 0, y: 0, k: 1 };
  let drag = null;

  data.elements.forEach((element) => {
    element.residue_ids.forEach((residueId) => {
      elementByResidue.set(residueId.toLowerCase(), element.id);
      elementByResidue.set(residueId.split(":").pop().toLowerCase(), element.id);
    });
  });

  function esc(value) {
    return String(value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#039;");
  }

  function makeSvg(tag, attrs) {
    const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
    Object.entries(attrs || {}).forEach(([key, value]) => {
      if (value === null || value === undefined) return;
      if (key === "textContent") node.textContent = value;
      else node.setAttribute(key, value);
    });
    return node;
  }

  function layoutElements() {
    const diagramWidth = 1180;
    const left = 42;
    const top = 54;
    const gap = 54;
    const rowGap = 118;
    let x = left;
    let y = top;

    data.elements.forEach((element, index) => {
      const width = Math.max(58, Math.min(150, 28 + element.length * 4.6));
      const height = element.type === "helix" ? 27 : 29;
      if (index > 0 && x + width + left > diagramWidth) {
        x = left;
        y += rowGap;
      }
      element.x = x;
      element.y = y;
      element.w = width;
      element.h = height;
      element.cx = x + width / 2;
      element.cy = y + height / 2;
      x += width + gap;
    });

    const diagramHeight = y + rowGap;
    svg.setAttribute("viewBox", `0 0 ${diagramWidth} ${diagramHeight}`);
  }

  function helixColor(code) {
    if (code === "G") return "var(--helix-alt)";
    if (code === "I") return "var(--helix-pi)";
    return "var(--helix)";
  }

  function strandColor(code) {
    if (code === "B") return "var(--strand-alt)";
    return "var(--strand)";
  }

  function connectorPath(a, b) {
    const ax = a.x + a.w;
    const ay = a.cy;
    const bx = b.x;
    const by = b.cy;
    if (Math.abs(ay - by) < 2) {
      return `M ${ax} ${ay} C ${ax + 22} ${ay}, ${bx - 22} ${by}, ${bx} ${by}`;
    }
    return `M ${ax} ${ay} C ${ax + 80} ${ay}, ${bx - 80} ${by}, ${bx} ${by}`;
  }

  function linkPath(a, b) {
    const lift = Math.abs(a.cy - b.cy) < 2 ? -58 : 0;
    const midY = (a.cy + b.cy) / 2 + lift;
    return `M ${a.cx} ${a.cy} C ${a.cx} ${midY}, ${b.cx} ${midY}, ${b.cx} ${b.cy}`;
  }

  function drawConnectors() {
    connectorsLayer.replaceChildren();
    for (let i = 0; i < data.elements.length - 1; i += 1) {
      const a = data.elements[i];
      const b = data.elements[i + 1];
      if (a.chain !== b.chain) continue;
      connectorsLayer.appendChild(makeSvg("path", {
        class: "connector",
        d: connectorPath(a, b)
      }));
    }
  }

  function drawLinks() {
    linksLayer.replaceChildren();
    data.links.forEach((link) => {
      const a = elementById.get(link.source);
      const b = elementById.get(link.target);
      if (!a || !b) return;
      const path = makeSvg("path", {
        class: "beta-link",
        d: linkPath(a, b),
        "data-link": `${link.source}:${link.target}`,
        "data-source": link.source,
        "data-target": link.target,
        "stroke-width": Math.min(7, 1.2 + Math.sqrt(link.count))
      });
      path.appendChild(makeSvg("title", {
        textContent: `${link.source} to ${link.target}: ${link.count} DSSP bridge contacts`
      }));
      linksLayer.appendChild(path);
    });
  }

  function drawHelix(group, element) {
    const shape = makeSvg("rect", {
      class: "sse-shape",
      x: element.x,
      y: element.y,
      width: element.w,
      height: element.h,
      rx: element.h / 2,
      fill: helixColor(element.ss_code)
    });
    group.appendChild(shape);

    for (let x = element.x + 11; x < element.x + element.w - 5; x += 18) {
      group.appendChild(makeSvg("path", {
        class: "helix-stripe",
        d: `M ${x} ${element.y + element.h - 4} C ${x + 7} ${element.y + 3}, ${x + 13} ${element.y + 3}, ${x + 18} ${element.y + element.h - 4}`
      }));
    }

    group.appendChild(makeSvg("rect", {
      class: "selection-ring",
      x: element.x - 5,
      y: element.y - 5,
      width: element.w + 10,
      height: element.h + 10,
      rx: element.h / 2 + 5
    }));
  }

  function drawStrand(group, element) {
    const head = Math.min(26, element.w * 0.35);
    const points = [
      `${element.x},${element.y}`,
      `${element.x + element.w - head},${element.y}`,
      `${element.x + element.w},${element.y + element.h / 2}`,
      `${element.x + element.w - head},${element.y + element.h}`,
      `${element.x},${element.y + element.h}`
    ].join(" ");
    group.appendChild(makeSvg("polygon", {
      class: "sse-shape",
      points: points,
      fill: strandColor(element.ss_code)
    }));
    group.appendChild(makeSvg("rect", {
      class: "selection-ring",
      x: element.x - 5,
      y: element.y - 5,
      width: element.w + 10,
      height: element.h + 10,
      rx: 6
    }));
  }

  function drawElements() {
    elementsLayer.replaceChildren();
    data.elements.forEach((element) => {
      const group = makeSvg("g", {
        class: "sse",
        "data-element-id": element.id,
        tabindex: 0
      });
      if (element.type === "helix") drawHelix(group, element);
      else drawStrand(group, element);

      group.appendChild(makeSvg("text", {
        class: "element-label",
        x: element.cx,
        y: element.y - 11,
        textContent: `${element.id} (${element.ss_code})`
      }));
      group.appendChild(makeSvg("text", {
        class: "residue-label",
        x: element.cx,
        y: element.y + element.h + 17,
        textContent: `${element.chain}:${element.start_residue}-${element.end_residue}`
      }));

      group.addEventListener("mouseenter", (event) => showTooltip(event, element));
      group.addEventListener("mousemove", (event) => moveTooltip(event));
      group.addEventListener("mouseleave", hideTooltip);
      group.addEventListener("click", (event) => {
        event.stopPropagation();
        selectElement(element.id);
      });
      group.addEventListener("keydown", (event) => {
        if (event.key === "Enter" || event.key === " ") selectElement(element.id);
      });
      elementsLayer.appendChild(group);
    });
  }

  function showTooltip(event, element) {
    tooltip.innerHTML = `
      <strong>${esc(element.id)} ${esc(element.ss_name)}</strong><br>
      Chain ${esc(element.chain)}, residues ${esc(element.start_residue)}-${esc(element.end_residue)}<br>
      Length ${element.length}; DSSP ${element.start_dssp}-${element.end_dssp}<br>
      ${esc(element.sequence.slice(0, 80))}${element.sequence.length > 80 ? "..." : ""}
    `;
    tooltip.style.display = "block";
    moveTooltip(event);
  }

  function moveTooltip(event) {
    const box = root.getBoundingClientRect();
    tooltip.style.left = `${event.clientX - box.left + 14}px`;
    tooltip.style.top = `${event.clientY - box.top + 14}px`;
  }

  function hideTooltip() {
    tooltip.style.display = "none";
  }

  function selectElement(elementId) {
    selectedId = elementId;
    const element = elementById.get(elementId);
    root.querySelectorAll("[data-element-id]").forEach((node) => {
      node.classList.toggle("selected", node.getAttribute("data-element-id") === elementId);
    });
    root.querySelectorAll("[data-link]").forEach((node) => {
      const linked = node.getAttribute("data-source") === elementId || node.getAttribute("data-target") === elementId;
      node.classList.toggle("link-selected", linked);
    });
    if (!element) return;
    const acc = element.accessibility_mean === null || element.accessibility_mean === undefined
      ? "not available"
      : element.accessibility_mean;
    details.innerHTML = `
      <strong>${esc(element.id)} ${esc(element.ss_name)}</strong>
      <span> chain <code>${esc(element.chain)}</code>, residues <code>${esc(element.start_residue)}-${esc(element.end_residue)}</code>,
      DSSP rows <code>${element.start_dssp}-${element.end_dssp}</code>, length <code>${element.length}</code>,
      mean ACC <code>${esc(acc)}</code>.</span>
      <br><span>Sequence: <code>${esc(element.sequence)}</code></span>
    `;
  }

  function jumpToResidue() {
    const query = residueInput.value.trim().toLowerCase();
    if (!query) return;
    const direct = elementByResidue.get(query);
    const withChain = query.includes(":") ? query : null;
    let elementId = direct;
    if (!elementId && withChain) elementId = elementByResidue.get(withChain);
    if (!elementId) {
      details.innerHTML = `No secondary-structure element contains residue <code>${esc(residueInput.value)}</code>.`;
      return;
    }
    selectElement(elementId);
  }

  function applyTransform() {
    viewport.setAttribute("transform", `translate(${transform.x} ${transform.y}) scale(${transform.k})`);
  }

  function zoomBy(factor) {
    transform.k = Math.max(0.25, Math.min(5, transform.k * factor));
    applyTransform();
  }

  function resetView() {
    transform = { x: 0, y: 0, k: 1 };
    applyTransform();
  }

  function downloadJson() {
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `${(data.name || "dssp-topology").replace(/\W+/g, "_")}.json`;
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(url);
  }

  root.querySelector("[data-action='zoom-in']").addEventListener("click", () => zoomBy(1.2));
  root.querySelector("[data-action='zoom-out']").addEventListener("click", () => zoomBy(1 / 1.2));
  root.querySelector("[data-action='reset']").addEventListener("click", resetView);
  root.querySelector("[data-action='jump']").addEventListener("click", jumpToResidue);
  root.querySelector("[data-action='download-json']").addEventListener("click", downloadJson);
  residueInput.addEventListener("keydown", (event) => {
    if (event.key === "Enter") jumpToResidue();
  });
  linkToggle.addEventListener("change", () => {
    linksLayer.style.display = linkToggle.checked ? "" : "none";
  });
  labelToggle.addEventListener("change", () => {
    const display = labelToggle.checked ? "" : "none";
    root.querySelectorAll(".element-label,.residue-label").forEach((node) => {
      node.style.display = display;
    });
  });

  svg.addEventListener("wheel", (event) => {
    event.preventDefault();
    zoomBy(event.deltaY < 0 ? 1.08 : 1 / 1.08);
  }, { passive: false });
  svg.addEventListener("mousedown", (event) => {
    if (event.button !== 0) return;
    drag = { x: event.clientX, y: event.clientY, moved: false };
    svg.classList.add("is-dragging");
  });
  window.addEventListener("mousemove", (event) => {
    if (!drag) return;
    const dx = event.clientX - drag.x;
    const dy = event.clientY - drag.y;
    if (Math.abs(dx) + Math.abs(dy) > 2) drag.moved = true;
    transform.x += dx;
    transform.y += dy;
    drag.x = event.clientX;
    drag.y = event.clientY;
    applyTransform();
  });
  window.addEventListener("mouseup", () => {
    drag = null;
    svg.classList.remove("is-dragging");
  });
  svg.addEventListener("click", () => {
    if (selectedId) return;
    root.querySelectorAll(".selected").forEach((node) => node.classList.remove("selected"));
  });

  layoutElements();
  drawConnectors();
  drawLinks();
  drawElements();
  applyTransform();
})();
""".replace(
        "__ROOT_ID__", root_id
    ).replace(
        "__DATA_ID__", data_id
    )

    return f"""
<div id="{root_id}">
  <style>{css}</style>
  <div class="topology-shell">
    <div class="topology-header">
      <div>
        <div class="title">{_html_escape(title)}</div>
        <div class="subtitle">{_html_escape(subtitle)}</div>
        <div class="stats">
          {stats["residue_count"]} residues;
          {stats["helix_count"]} helices;
          {stats["strand_count"]} strands;
          {stats["beta_link_count"]} beta links;
          chains: {_html_escape(", ".join(stats["chains"]))}
        </div>
      </div>
      <div class="controls" aria-label="Topology controls">
        <button type="button" data-action="zoom-in" title="Zoom in">+</button>
        <button type="button" data-action="zoom-out" title="Zoom out">-</button>
        <button type="button" data-action="reset" title="Reset pan and zoom">Reset</button>
        <label class="toggle" title="Show DSSP beta bridge links">
          <input type="checkbox" data-role="toggle-links" checked> Links
        </label>
        <label class="toggle" title="Show element and residue labels">
          <input type="checkbox" data-role="toggle-labels" checked> Labels
        </label>
        <input type="text" data-role="residue-input" placeholder="Residue A:120">
        <button type="button" data-action="jump" title="Select a residue">Find</button>
        <button type="button" data-action="download-json" title="Download parsed topology JSON">JSON</button>
      </div>
    </div>
    <div class="viewer-wrap">
      <svg role="img" aria-label="Interactive DSSP topology diagram">
        <g data-role="viewport">
          <g data-role="connectors"></g>
          <g data-role="links"></g>
          <g data-role="elements"></g>
        </g>
      </svg>
      <div class="tooltip" data-role="tooltip"></div>
    </div>
    <div class="details" data-role="details">
      Hover over an element for DSSP details. Click an element to keep its sequence and bridge links selected.
    </div>
  </div>
  <script type="application/json" id="{data_id}">{json_blob}</script>
  <script>{script}</script>
</div>
"""


def _clean_alphafold_topology_html(topology: Dict[str, Any]) -> str:
    root_id = "clean-topology-" + uuid.uuid4().hex
    data_id = root_id + "-data"
    title = topology.get("name") or topology.get("pdbe_entry_id") or "Topology"
    stats = topology["stats"]
    metadata = topology.get("metadata", {})
    subtitle = metadata.get("source") or "Generated topology"
    json_blob = json.dumps(topology, separators=(",", ":")).replace("</", "<\\/")
    molstar_panel = ""
    if topology.get("cif_url") or topology.get("afdb_accession"):
        molstar_panel = """
        <section class="clean-molstar-panel">
          <div class="panel-bar">Mol* 3D view</div>
          <div class="molstar-stage" data-role="molstar"></div>
          <div class="molstar-status" data-role="molstar-status">Loading Mol*...</div>
        </section>
        """

    css = """
#__ROOT_ID__ {
  --ink: #18202a;
  --muted: #627080;
  --line: #d5dde7;
  --panel: #f4f7fa;
  --helix: #d95f72;
  --strand: #2f7fbf;
  --strand-dark: #1f5e94;
  --sheet: #e7f2fb;
  --select: #173f8a;
  color: var(--ink);
  font-family: Inter, "Segoe UI", Arial, sans-serif;
  line-height: 1.35;
}
#__ROOT_ID__ .clean-shell {
  background: #ffffff;
  border: 1px solid #c9d4c1;
}
#__ROOT_ID__ .clean-header {
  background: #dcefd4;
  border-bottom: 1px solid #b9d9ad;
  padding: 12px 14px;
}
#__ROOT_ID__ .title {
  color: #2f6f39;
  font-size: 21px;
  font-weight: 700;
}
#__ROOT_ID__ .subtitle,
#__ROOT_ID__ .stats {
  color: #3e5945;
  font-size: 12px;
  margin-top: 3px;
}
#__ROOT_ID__ .clean-grid {
  display: grid;
  gap: 12px;
  grid-template-columns: minmax(520px, 1.15fr) minmax(360px, 0.85fr);
  padding: 12px;
}
#__ROOT_ID__ .topology-panel {
  border: 1px solid #bfc8d3;
  min-width: 0;
}
#__ROOT_ID__ .toolbar {
  align-items: center;
  background: #68717a;
  color: #ffffff;
  display: flex;
  flex-wrap: wrap;
  gap: 7px;
  min-height: 42px;
  padding: 6px 10px;
}
#__ROOT_ID__ .toolbar .label {
  font-size: 15px;
  font-weight: 600;
  margin-right: auto;
  min-width: 190px;
}
#__ROOT_ID__ button,
#__ROOT_ID__ select,
#__ROOT_ID__ input[type="text"] {
  background: #ffffff;
  border: 1px solid #c8d0d8;
  border-radius: 4px;
  color: #111827;
  font: inherit;
  font-size: 12px;
  min-height: 28px;
}
#__ROOT_ID__ button {
  cursor: pointer;
  padding: 3px 8px;
}
#__ROOT_ID__ input[type="text"] {
  width: 120px;
  padding: 3px 7px;
}
#__ROOT_ID__ select {
  padding: 3px 7px;
}
#__ROOT_ID__ .toggle {
  align-items: center;
  display: inline-flex;
  font-size: 12px;
  gap: 4px;
  min-height: 28px;
  white-space: nowrap;
}
#__ROOT_ID__ .viewer-wrap {
  background: linear-gradient(#ffffff, #fbfcfe);
  height: 690px;
  position: relative;
}
#__ROOT_ID__ svg {
  cursor: grab;
  display: block;
  height: 100%;
  width: 100%;
}
#__ROOT_ID__ svg.dragging {
  cursor: grabbing;
}
#__ROOT_ID__ .domain-band {
  fill: #f8fafc;
  stroke: #edf1f5;
  stroke-width: 1;
}
#__ROOT_ID__ .domain-label {
  fill: #7b8794;
  font-size: 11px;
  font-weight: 700;
  letter-spacing: 0;
}
#__ROOT_ID__ .seq-connector {
  fill: none;
  opacity: 0.78;
  stroke: #1f2933;
  stroke-linecap: round;
  stroke-linejoin: round;
  stroke-width: 1.7;
}
#__ROOT_ID__ .contact-link {
  fill: none;
  opacity: 0.22;
  stroke: #2f7fbf;
  stroke-linecap: round;
}
#__ROOT_ID__ .contact-link.hot {
  opacity: 0.9;
  stroke: #173f8a;
}
#__ROOT_ID__ .sse {
  cursor: pointer;
  outline: none;
}
#__ROOT_ID__ .sse-shape {
  filter: drop-shadow(0 1px 1px rgba(17, 24, 39, 0.12));
  stroke: rgba(17, 24, 39, 0.35);
  stroke-width: 1.2;
}
#__ROOT_ID__ .helix .sse-shape {
  fill: var(--helix);
}
#__ROOT_ID__ .strand .sse-shape {
  fill: var(--strand);
}
#__ROOT_ID__ .sse:hover .sse-shape,
#__ROOT_ID__ .sse.selected .sse-shape {
  stroke: var(--select);
  stroke-width: 2.4;
}
#__ROOT_ID__ .sse.molstar-hot .sse-shape {
  stroke: #f59e0b;
  stroke-width: 2.4;
}
#__ROOT_ID__ .sheet-chip {
  fill: #ffffff;
  opacity: 0.78;
  stroke: rgba(47, 127, 191, 0.38);
  stroke-width: 1;
}
#__ROOT_ID__ .sheet-text {
  fill: #265f8f;
  font-size: 8px;
  font-weight: 700;
  pointer-events: none;
  text-anchor: middle;
}
#__ROOT_ID__ .sse-label {
  fill: #111827;
  font-size: 11px;
  font-weight: 800;
  pointer-events: none;
  text-anchor: middle;
}
#__ROOT_ID__ .range-label {
  fill: #667280;
  font-size: 9px;
  pointer-events: none;
  text-anchor: middle;
}
#__ROOT_ID__ .residue-tick {
  display: none;
  pointer-events: none;
  stroke: #ffe66d;
  stroke-width: 2.4;
}
#__ROOT_ID__ .sse.has-residue .residue-tick {
  display: block;
}
#__ROOT_ID__ .terminus {
  font-size: 15px;
  font-weight: 800;
  pointer-events: none;
}
#__ROOT_ID__ .terminus.n {
  fill: #0b39ff;
}
#__ROOT_ID__ .terminus.c {
  fill: #e11919;
}
#__ROOT_ID__ .tooltip {
  background: #111827;
  border-radius: 5px;
  color: #ffffff;
  display: none;
  font-size: 12px;
  max-width: 320px;
  padding: 8px 10px;
  pointer-events: none;
  position: absolute;
  z-index: 5;
}
#__ROOT_ID__ .details {
  border-top: 1px solid #d7dfe8;
  color: var(--muted);
  font-size: 13px;
  min-height: 46px;
  padding: 10px 12px;
}
#__ROOT_ID__ .details strong {
  color: var(--ink);
}
#__ROOT_ID__ .details code {
  background: #eef2f7;
  border-radius: 3px;
  color: #344154;
  padding: 1px 4px;
}
#__ROOT_ID__ .clean-molstar-panel {
  border: 1px solid #b9c5d1;
  min-width: 0;
  position: relative;
}
#__ROOT_ID__ .panel-bar {
  background: #f0f3f6;
  border-bottom: 1px solid #c9d1da;
  font-size: 13px;
  min-height: 34px;
  padding: 8px 10px;
}
#__ROOT_ID__ .molstar-stage {
  height: 690px;
  position: relative;
}
#__ROOT_ID__ .molstar-stage iframe {
  border: 0;
  display: block;
  height: 100%;
  width: 100%;
}
#__ROOT_ID__ .molstar-status {
  background: rgba(255, 255, 255, 0.92);
  bottom: 8px;
  color: #344154;
  font-size: 12px;
  left: 8px;
  padding: 4px 6px;
  position: absolute;
}
@media (max-width: 1050px) {
  #__ROOT_ID__ .clean-grid {
    grid-template-columns: 1fr;
  }
}
@media (max-width: 760px) {
  #__ROOT_ID__ .clean-grid {
    padding: 8px;
  }
  #__ROOT_ID__ .viewer-wrap,
  #__ROOT_ID__ .molstar-stage {
    height: 560px;
  }
}
""".replace(
        "__ROOT_ID__", root_id
    )

    script = r"""
(function () {
  const root = document.getElementById("__ROOT_ID__");
  const data = JSON.parse(document.getElementById("__DATA_ID__").textContent);
  const svg = root.querySelector("svg");
  const viewport = root.querySelector("[data-role='viewport']");
  const bandsLayer = root.querySelector("[data-role='bands']");
  const connectorsLayer = root.querySelector("[data-role='connectors']");
  const linksLayer = root.querySelector("[data-role='links']");
  const elementsLayer = root.querySelector("[data-role='elements']");
  const tooltip = root.querySelector("[data-role='tooltip']");
  const details = root.querySelector("[data-role='details']");
  const linkToggle = root.querySelector("[data-role='toggle-links']");
  const labelToggle = root.querySelector("[data-role='toggle-labels']");
  const clickMode = root.querySelector("[data-role='click-mode']");
  const residueInput = root.querySelector("[data-role='residue-input']");
  const molstarNode = root.querySelector("[data-role='molstar']");
  const molstarStatus = root.querySelector("[data-role='molstar-status']");
  const residues = data.residues || [];
  const residueBySeq = new Map(residues.map((residue) => [Number(residue.seq), residue]));
  const elementById = new Map(data.elements.map((element) => [element.id, element]));
  const elementByResidue = new Map();
  let selectedId = null;
  let transform = { x: 0, y: 0, k: 1 };
  let drag = null;
  let molstarViewer = null;
  let lastMolstarHit = "";
  const MOLSTAR_BASE_GREY = 0xb8bec7;
  const MOLSTAR_ACTIVE_RED = 0xe11919;
  const SHEET_PALETTE = [
    "#f2a541", "#84cce6", "#f2644a", "#80b86c", "#9a86d4",
    "#d99a5f", "#4aa3a2", "#c879b2", "#7da0d8", "#d4b94f"
  ];

  data.elements.forEach((element) => {
    element.display_id = element.type === "strand"
      ? element.id.replace(/^S/, "E")
      : element.id;
    element.residue_indices.forEach((seq) => {
      elementByResidue.set(String(seq).toLowerCase(), element.id);
      elementByResidue.set(`${element.chain}:${seq}`.toLowerCase(), element.id);
    });
  });

  function esc(value) {
    return String(value ?? "")
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#039;");
  }

  function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
  }

  function makeSvg(tag, attrs) {
    const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
    Object.entries(attrs || {}).forEach(([key, value]) => {
      if (value === null || value === undefined) return;
      if (key === "textContent") node.textContent = value;
      else node.setAttribute(key, value);
    });
    return node;
  }

  function sheetGroups() {
    const strandIds = new Set(data.elements.filter((item) => item.type === "strand").map((item) => item.id));
    const adjacency = new Map([...strandIds].map((id) => [id, []]));
    data.links.forEach((link) => {
      if (!strandIds.has(link.source) || !strandIds.has(link.target)) return;
      adjacency.get(link.source).push(link.target);
      adjacency.get(link.target).push(link.source);
    });
    const groups = [];
    const seen = new Set();
    data.elements.forEach((element) => {
      if (element.type !== "strand" || seen.has(element.id)) return;
      const stack = [element.id];
      const group = [];
      seen.add(element.id);
      while (stack.length) {
        const current = stack.pop();
        group.push(current);
        adjacency.get(current).forEach((next) => {
          if (!seen.has(next)) {
            seen.add(next);
            stack.push(next);
          }
        });
      }
      groups.push(group.sort((a, b) => elementById.get(a).start - elementById.get(b).start));
    });
    groups.sort((a, b) => elementById.get(a[0]).start - elementById.get(b[0]).start);
    return groups;
  }

  const sheetIndexById = new Map();
  sheetGroups().forEach((group, index) => {
    group.forEach((id) => sheetIndexById.set(id, index + 1));
  });

  function strandFill(element) {
    const sheetNo = sheetIndexById.get(element.id);
    if (!sheetNo) return "#7da9c7";
    return SHEET_PALETTE[(sheetNo - 1) % SHEET_PALETTE.length];
  }

  function layout() {
    const columns = data.elements.length > 85 ? 12 : data.elements.length > 55 ? 10 : 9;
    const cellW = 98;
    const rowH = 154;
    const left = 76;
    const top = 78;
    const rows = Math.ceil(data.elements.length / columns);

    data.elements.forEach((element, index) => {
      const row = Math.floor(index / columns);
      const rawCol = index % columns;
      const col = row % 2 === 0 ? rawCol : columns - rawCol - 1;
      const local = row % 2 === 0 ? rawCol : columns - rawCol - 1;
      const len = element.stop - element.start + 1;
      element.row = row;
      element.col = col;
      element.w = element.type === "helix" ? 36 : 42;
      element.h = clamp(44 + len * 3.8, 56, 104);
      element.x = left + col * cellW;
      element.y = top + row * rowH + (element.type === "helix" ? 6 : 0) + ((local % 2) - 0.5) * 5;
      element.direction = (index + row) % 2 === 0 ? 1 : -1;
      element.cx = element.x + element.w / 2;
      element.cy = element.y + element.h / 2;
    });

    const width = left * 2 + (columns - 1) * cellW + 104;
    const height = top + rows * rowH + 64;
    svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
    svg.dataset.viewWidth = width;
    svg.dataset.viewHeight = height;
  }

  function startPoint(element) {
    return element.direction > 0
      ? { x: element.cx, y: element.y }
      : { x: element.cx, y: element.y + element.h };
  }

  function endPoint(element) {
    return element.direction > 0
      ? { x: element.cx, y: element.y + element.h }
      : { x: element.cx, y: element.y };
  }

  function outsideTerminal(element, terminal) {
    const point = terminal === "start" ? startPoint(element) : endPoint(element);
    const isTop = Math.abs(point.y - element.y) < Math.abs(point.y - (element.y + element.h));
    return {
      x: point.x,
      y: point.y + (isTop ? -16 : 16)
    };
  }

  function pathForConnector(a, b) {
    const start = endPoint(a);
    const end = startPoint(b);
    const startOut = outsideTerminal(a, "end");
    const endOut = outsideTerminal(b, "start");
    const viewWidth = Number(svg.dataset.viewWidth) || 1200;
    const routeX = a.row === b.row
      ? (a.cx < b.cx
        ? (a.x + a.w + b.x) / 2
        : (b.x + b.w + a.x) / 2)
      : (a.cx > viewWidth / 2 ? viewWidth - 42 : 42);
    return [
      `M ${start.x} ${start.y}`,
      `L ${startOut.x} ${startOut.y}`,
      `L ${routeX} ${startOut.y}`,
      `L ${routeX} ${endOut.y}`,
      `L ${endOut.x} ${endOut.y}`,
      `L ${end.x} ${end.y}`
    ].join(" ");
  }

  function pathForContact(a, b) {
    const dx = Math.abs(a.cx - b.cx);
    const lift = Math.max(26, Math.min(88, dx * 0.22));
    const sameRow = Math.abs(a.cy - b.cy) < 80;
    const midY = sameRow ? Math.min(a.y, b.y) - lift : (a.cy + b.cy) / 2;
    return `M ${a.cx} ${a.cy} C ${a.cx} ${midY}, ${b.cx} ${midY}, ${b.cx} ${b.cy}`;
  }

  function strandPoints(element) {
    const x = element.x;
    const y = element.y;
    const w = element.w;
    const h = element.h;
    const head = Math.min(30, Math.max(20, h * 0.26));
    const inset = Math.max(7, w * 0.24);
    if (element.direction > 0) {
      return [
        `${x + inset},${y}`,
        `${x + w - inset},${y}`,
        `${x + w - inset},${y + h - head}`,
        `${x + w},${y + h - head}`,
        `${x + w / 2},${y + h}`,
        `${x},${y + h - head}`,
        `${x + inset},${y + h - head}`
      ].join(" ");
    }
    return [
      `${x + w / 2},${y}`,
      `${x + w},${y + head}`,
      `${x + w - inset},${y + head}`,
      `${x + w - inset},${y + h}`,
      `${x + inset},${y + h}`,
      `${x + inset},${y + head}`,
      `${x},${y + head}`
    ].join(" ");
  }

  function drawBands() {
    bandsLayer.replaceChildren();
    const rows = new Map();
    data.elements.forEach((element) => {
      const rowKey = element.row;
      if (!rows.has(rowKey)) rows.set(rowKey, []);
      rows.get(rowKey).push(element);
    });
    rows.forEach((items, rowKey) => {
      const minX = Math.min(...items.map((item) => item.x)) - 24;
      const maxX = Math.max(...items.map((item) => item.x + item.w)) + 24;
      const y = Math.min(...items.map((item) => item.y)) - 42;
      const h = Math.max(...items.map((item) => item.y + item.h)) - y + 38;
      bandsLayer.appendChild(makeSvg("rect", {
        class: "domain-band",
        x: minX,
        y,
        width: maxX - minX,
        height: h,
        rx: 8
      }));
      bandsLayer.appendChild(makeSvg("text", {
        class: "domain-label",
        x: minX + 10,
        y: y + 18,
        textContent: `segment ${rowKey + 1}`
      }));
    });
  }

  function drawConnectors() {
    connectorsLayer.replaceChildren();
    for (let i = 0; i < data.elements.length - 1; i += 1) {
      connectorsLayer.appendChild(makeSvg("path", {
        class: "seq-connector",
        d: pathForConnector(data.elements[i], data.elements[i + 1])
      }));
    }
    if (!data.elements.length) return;
    const first = startPoint(data.elements[0]);
    const last = endPoint(data.elements[data.elements.length - 1]);
    connectorsLayer.appendChild(makeSvg("text", {
      class: "terminus n",
      x: first.x - 18,
      y: first.y + 5,
      textContent: "N"
    }));
    connectorsLayer.appendChild(makeSvg("text", {
      class: "terminus c",
      x: last.x + 14,
      y: last.y + 5,
      textContent: "C"
    }));
  }

  function drawLinks() {
    linksLayer.replaceChildren();
    data.links.forEach((link) => {
      const a = elementById.get(link.source);
      const b = elementById.get(link.target);
      if (!a || !b) return;
      const path = makeSvg("path", {
        class: "contact-link",
        d: pathForContact(a, b),
        "data-source": link.source,
        "data-target": link.target,
        "stroke-width": Math.min(4.5, 1 + Math.sqrt(link.count) / 2.4)
      });
      path.appendChild(makeSvg("title", {
        textContent: `${a.display_id} to ${b.display_id}: ${link.count} inferred contacts`
      }));
      linksLayer.appendChild(path);
    });
    linksLayer.style.display = linkToggle.checked ? "" : "none";
  }

  function drawElements() {
    elementsLayer.replaceChildren();
    data.elements.forEach((element) => {
      const group = makeSvg("g", {
        class: `sse ${element.type}`,
        "data-element-id": element.id,
        tabindex: 0
      });
      if (element.type === "helix") {
        group.appendChild(makeSvg("rect", {
          class: "sse-shape",
          x: element.x,
          y: element.y,
          width: element.w,
          height: element.h,
          rx: element.w / 2,
          ry: element.w / 2,
          style: "fill: var(--helix)"
        }));
      } else {
        group.appendChild(makeSvg("polygon", {
          class: "sse-shape",
          points: strandPoints(element),
          style: `fill: ${strandFill(element)}`
        }));
        const sheetNo = sheetIndexById.get(element.id);
        if (sheetNo) {
          group.appendChild(makeSvg("rect", {
            class: "sheet-chip",
            x: element.x + 5,
            y: element.y + element.h / 2 - 8,
            width: element.w - 10,
            height: 16,
            rx: 6
          }));
          group.appendChild(makeSvg("text", {
            class: "sheet-text",
            x: element.cx,
            y: element.y + element.h / 2 + 3,
            textContent: `B${sheetNo}`
          }));
        }
      }
      group.appendChild(makeSvg("line", {
        class: "residue-tick",
        x1: element.x - 6,
        x2: element.x + element.w + 6,
        y1: element.cy,
        y2: element.cy
      }));
      group.appendChild(makeSvg("text", {
        class: "sse-label",
        x: element.cx,
        y: element.y - 10,
        textContent: element.display_id
      }));
      group.appendChild(makeSvg("text", {
        class: "range-label",
        x: element.cx,
        y: element.y + element.h + 14,
        textContent: `${element.start}-${element.stop}`
      }));

      group.addEventListener("mouseenter", (event) => showResidue(event, element, group));
      group.addEventListener("mousemove", (event) => showResidue(event, element, group));
      group.addEventListener("mouseleave", () => {
        group.classList.remove("has-residue");
        tooltip.style.display = "none";
        try { molstarViewer?.plugin?.managers?.interactivity?.clearHighlights?.(); } catch (error) {}
      });
      group.addEventListener("click", (event) => {
        event.stopPropagation();
        selectElement(element.id, residueAtEvent(event, element));
      });
      group.addEventListener("keydown", (event) => {
        if (event.key === "Enter" || event.key === " ") selectElement(element.id, residueBySeq.get(element.start));
      });
      elementsLayer.appendChild(group);
    });
  }

  function residueAtEvent(event, element) {
    const count = Math.max(1, element.stop - element.start + 1);
    const point = svg.createSVGPoint();
    point.x = event.clientX;
    point.y = event.clientY;
    const matrix = viewport.getScreenCTM();
    if (!matrix) return residueBySeq.get(element.start);
    const local = point.matrixTransform(matrix.inverse());
    const t = clamp((local.y - element.y) / element.h, 0, 1);
    const oriented = element.direction > 0 ? t : 1 - t;
    const seq = element.start + Math.round(oriented * (count - 1));
    return residueBySeq.get(seq) || residueBySeq.get(element.start);
  }

  function updateTick(group, element, residue) {
    if (!residue) return;
    const count = Math.max(1, element.stop - element.start);
    const t = count === 0 ? 0.5 : (residue.seq - element.start) / count;
    const oriented = element.direction > 0 ? t : 1 - t;
    const y = element.y + clamp(oriented, 0, 1) * element.h;
    const tick = group.querySelector(".residue-tick");
    tick.setAttribute("y1", y);
    tick.setAttribute("y2", y);
    group.classList.add("has-residue");
  }

  function showResidue(event, element, group) {
    const residue = residueAtEvent(event, element);
    updateTick(group, element, residue);
    if (!residue) return;
    const plddt = residue.plddt == null ? "not available" : Number(residue.plddt).toFixed(2);
    const sheetNo = sheetIndexById.get(element.id);
    tooltip.innerHTML = `
      <strong>${esc(residue.residue_id)} ${esc(residue.aa)}</strong><br>
      ${esc(element.display_id)} ${esc(element.ss_name)}${sheetNo ? `, sheet B${sheetNo}` : ""}<br>
      Residues ${esc(element.start)}-${esc(element.stop)}; pLDDT ${esc(plddt)}
    `;
    tooltip.style.display = "block";
    const box = root.getBoundingClientRect();
    tooltip.style.left = `${event.clientX - box.left + 14}px`;
    tooltip.style.top = `${event.clientY - box.top + 14}px`;
    highlightMolstarResidue(residue, "highlight");
  }

  function nodeForElement(elementId) {
    let found = null;
    root.querySelectorAll("[data-element-id]").forEach((node) => {
      if (node.getAttribute("data-element-id") === elementId) found = node;
    });
    return found;
  }

  function elementForSeq(seq) {
    const number = Number(seq);
    return data.elements.find((element) => number >= element.start && number <= element.stop);
  }

  function sheetElementsFor(element) {
    if (!element || element.type !== "strand") return element ? [element] : [];
    const sheetNo = sheetIndexById.get(element.id);
    if (!sheetNo) return [element];
    return data.elements.filter((item) => item.type === "strand" && sheetIndexById.get(item.id) === sheetNo);
  }

  function selectionFor(element, residue, requestedMode) {
    const mode = requestedMode || clickMode?.value || "residue";
    if (mode === "residue" && residue) {
      return { mode: "residue", elements: element ? [element] : [] };
    }
    if (mode === "sheet" && element?.type === "strand") {
      return { mode: "sheet", elements: sheetElementsFor(element) };
    }
    return { mode: "range", elements: element ? [element] : [] };
  }

  function selectElement(elementId, residue, options = {}) {
    selectedId = elementId;
    const element = elementById.get(elementId);
    const selection = selectionFor(element, residue, options.mode);
    const selectedIds = new Set(selection.elements.map((item) => item.id));
    root.querySelectorAll("[data-element-id]").forEach((node) => {
      const id = node.getAttribute("data-element-id");
      node.classList.toggle("selected", selectedIds.has(id));
      node.classList.remove("molstar-hot");
    });
    root.querySelectorAll(".contact-link").forEach((node) => {
      const linked = selectedIds.has(node.getAttribute("data-source")) || selectedIds.has(node.getAttribute("data-target"));
      node.classList.toggle("hot", linked);
    });
    if (!element) return;
    const linked = data.links
      .filter((link) => selectedIds.has(link.source) || selectedIds.has(link.target))
      .map((link) => {
        const other = selectedIds.has(link.source) ? link.target : link.source;
        return elementById.get(other)?.display_id || other;
      })
      .filter((id, index, list) => list.indexOf(id) === index);
    const sheetNo = sheetIndexById.get(element.id);
    const selectionLabel = selection.mode === "residue"
      ? "Mol*: residue"
      : selection.mode === "sheet"
        ? (sheetNo ? `Mol*: sheet B${sheetNo}` : "Mol*: sheet")
        : "Mol*: SSE range";
    const ranges = selection.mode === "residue" && residue
      ? `<code>${esc(residue.seq)}-${esc(residue.seq)}</code>`
      : selection.elements.length === 1
        ? `<code>${esc(element.start)}-${esc(element.stop)}</code>`
        : selection.elements.map((item) => `${esc(item.display_id)} <code>${esc(item.start)}-${esc(item.stop)}</code>`).join(", ");
    const totalLength = selection.mode === "residue" && residue
      ? 1
      : selection.elements.reduce((sum, item) => sum + Number(item.length || 0), 0) || element.length;
    const plddt = residue?.plddt == null ? "not available" : Number(residue.plddt).toFixed(2);
    details.innerHTML = `
      <strong>${esc(element.display_id)} ${esc(element.ss_name)}</strong>
      <span>${esc(selectionLabel)}; residues ${ranges}, length <code>${esc(totalLength)}</code>,
      contacts <code>${esc(linked.join(", ") || "none")}</code>.</span>
      ${residue ? `<br><span>Residue <code>${esc(residue.residue_id)}</code> ${esc(residue.aa)}, pLDDT <code>${esc(plddt)}</code>.</span>` : ""}
    `;
    if (options.skipMolstar) return;
    applyMolstarSelection(selection, residue);
  }

  function jumpToResidue() {
    const query = residueInput.value.trim().toLowerCase();
    if (!query) return;
    const elementId = elementByResidue.get(query);
    if (!elementId) {
      details.innerHTML = `No helix or strand contains residue <code>${esc(residueInput.value)}</code>.`;
      return;
    }
    const seq = Number(query.split(":").pop());
    selectElement(elementId, residueBySeq.get(seq), { mode: "residue" });
  }

  function applyTransform() {
    viewport.setAttribute("transform", `translate(${transform.x} ${transform.y}) scale(${transform.k})`);
  }

  function zoomBy(factor) {
    transform.k = clamp(transform.k * factor, 0.35, 5);
    applyTransform();
  }

  function resetView() {
    transform = { x: 0, y: 0, k: 1 };
    applyTransform();
  }

  function downloadJson() {
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `${(data.name || "alphafold-topology").replace(/\W+/g, "_")}.json`;
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(url);
  }

  function loadScriptOnce(url, id) {
    if (id === "molstar-viewer-js" && window.molstar) return Promise.resolve();
    const existing = document.getElementById(id);
    if (existing?.dataset.ready === "true") return Promise.resolve();
    if (existing) {
      return new Promise((resolve, reject) => {
        existing.addEventListener("load", resolve, { once: true });
        existing.addEventListener("error", reject, { once: true });
      });
    }
    return new Promise((resolve, reject) => {
      const script = document.createElement("script");
      script.id = id;
      script.src = url;
      script.onload = () => {
        script.dataset.ready = "true";
        resolve();
      };
      script.onerror = reject;
      document.head.appendChild(script);
    });
  }

  function loadCssOnce(url, id) {
    if (document.getElementById(id)) return;
    const link = document.createElement("link");
    link.id = id;
    link.rel = "stylesheet";
    link.href = url;
    document.head.appendChild(link);
  }

  function molstarLoadOptions() {
    return {
      representationParams: {
        theme: {
          globalName: "uniform",
          globalColorParams: { value: MOLSTAR_BASE_GREY },
          carbonColor: { name: "uniform", params: { value: MOLSTAR_BASE_GREY } }
        }
      }
    };
  }

  function configureMolstarAppearance() {
    const canvas3d = molstarViewer?.plugin?.canvas3d;
    if (!canvas3d) return;
    try {
      const renderer = canvas3d.props?.renderer || {};
      canvas3d.setProps({
        renderer: {
          ...renderer,
          highlightColor: MOLSTAR_ACTIVE_RED,
          selectColor: MOLSTAR_ACTIVE_RED
        }
      });
    } catch (error) {
      if (molstarStatus) molstarStatus.textContent = `Mol* color setup failed: ${error.message || error}`;
    }
  }

  async function initMolstar() {
    if (!molstarNode) return;
    try {
      if (molstarModelUrl()) {
        loadMolstarIframeMvs([], { mode: "initial", elements: [] }, null);
        if (molstarStatus) molstarStatus.textContent = "Mol* view is driven by MolViewSpec. Click a topology residue to make it red.";
        return;
      }
      loadCssOnce("https://cdn.jsdelivr.net/npm/molstar@5.4.2/build/viewer/molstar.css", "molstar-viewer-css");
      await loadScriptOnce("https://cdn.jsdelivr.net/npm/molstar@5.4.2/build/viewer/molstar.js", "molstar-viewer-js");
      molstarViewer = await window.molstar.Viewer.create(molstarNode, {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: true,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        viewportBackgroundColor: "white"
      });
      configureMolstarAppearance();
      if (mvsAvailable() && molstarModelUrl()) {
        await loadMolstarMvs([], { mode: "initial", elements: [] }, null, { replaceExisting: true, focus: false });
      } else if (data.cif_url) {
        await molstarViewer.loadStructureFromUrl(data.cif_url, "mmcif", false, molstarLoadOptions());
      } else if (data.afdb_accession) {
        const afdbCifUrl = `https://alphafold.ebi.ac.uk/files/AF-${data.afdb_accession}-F1-model_v6.cif`;
        try {
          await molstarViewer.loadStructureFromUrl(afdbCifUrl, "mmcif", false, molstarLoadOptions());
        } catch (loadError) {
          if (!molstarViewer.loadAlphaFoldDb) throw loadError;
          await molstarViewer.loadAlphaFoldDb(`AF-${data.afdb_accession}-F1`);
        }
      }
      configureMolstarAppearance();
      const reverseLinked = bindMolstarToTopology();
      if (molstarStatus) {
        molstarStatus.textContent = reverseLinked
          ? "Click topology residues to rebuild the Mol* view with a red MVS selection."
          : "Click topology residues to highlight them in Mol*.";
      }
    } catch (error) {
      if (molstarStatus) molstarStatus.textContent = `Mol* could not load: ${error.message || error}`;
    }
  }

  function molstarRange(chain, start, stop) {
    const chainId = chain || "A";
    const beg = Number(start);
    const end = Number(stop);
    return {
      label: {
        label_asym_id: chainId,
        beg_label_seq_id: beg,
        end_label_seq_id: end
      },
      auth: {
        auth_asym_id: chainId,
        beg_auth_seq_id: beg,
        end_auth_seq_id: end
      },
      labelAnyChain: {
        beg_label_seq_id: beg,
        end_label_seq_id: end
      },
      authAnyChain: {
        beg_auth_seq_id: beg,
        end_auth_seq_id: end
      }
    };
  }

  function molstarModelUrl() {
    if (data.cif_url) return data.cif_url;
    if (data.afdb_accession) return `https://alphafold.ebi.ac.uk/files/AF-${data.afdb_accession}-F1-model_v6.cif`;
    return "";
  }

  function mvsAvailable() {
    return Boolean(window.molstar?.PluginExtensions?.mvs?.MVSData?.createBuilder && window.molstar?.PluginExtensions?.mvs?.loadMVS);
  }

  function mvsSelectorForRange(range, exactResidue) {
    const chain = range.chain || "A";
    const start = Number(range.start);
    const stop = Number(range.stop);
    const selectors = exactResidue && start === stop
      ? [
        { label_asym_id: chain, label_seq_id: start },
        { auth_asym_id: chain, auth_seq_id: start }
      ]
      : [
        { label_asym_id: chain, beg_label_seq_id: start, end_label_seq_id: stop },
        { auth_asym_id: chain, beg_auth_seq_id: start, end_auth_seq_id: stop }
      ];
    if ((data.stats?.chains || []).length <= 1) {
      if (exactResidue && start === stop) {
        selectors.push({ label_seq_id: start });
        selectors.push({ auth_seq_id: start });
      } else {
        selectors.push({ beg_label_seq_id: start, end_label_seq_id: stop });
        selectors.push({ beg_auth_seq_id: start, end_auth_seq_id: stop });
      }
    }
    return selectors;
  }

  function mvsSelectorForRanges(ranges, exactResidue) {
    const selectors = ranges.flatMap((range) => mvsSelectorForRange(range, exactResidue));
    if (!selectors.length) return undefined;
    return selectors.length === 1 ? selectors[0] : selectors;
  }

  function clampSeq(seq) {
    const values = residues.map((residue) => Number(residue.seq)).filter(Number.isFinite);
    if (!values.length) return Number(seq);
    return clamp(Number(seq), Math.min(...values), Math.max(...values));
  }

  function focusRangesForSelection(ranges, selection) {
    if (selection.mode !== "residue" || !ranges.length) return ranges;
    const range = ranges[0];
    return [{
      chain: range.chain,
      start: clampSeq(Number(range.start) - 20),
      stop: clampSeq(Number(range.stop) + 20)
    }];
  }

  function mvsNode(kind, params = {}, children = []) {
    return { kind, params, children };
  }

  function buildMvsJson(ranges, selection, residue) {
    const url = molstarModelUrl();
    if (!url) return null;
    const structureChildren = [
      mvsNode("component", { selector: "polymer" }, [
        mvsNode("representation", { type: "cartoon" }, [
          mvsNode("color", { color: "#b8bec7" })
        ])
      ])
    ];
    if (ranges.length) {
      const exactResidue = selection.mode === "residue" && residue;
      const focusSelector = mvsSelectorForRanges(focusRangesForSelection(ranges, selection), false);
      if (focusSelector) {
        structureChildren.push(mvsNode("component", { selector: focusSelector }, [
          mvsNode("focus", {})
        ]));
      }
      const activeSelector = mvsSelectorForRanges(ranges, exactResidue);
      structureChildren.push(mvsNode("component", { selector: activeSelector }, [
        mvsNode("representation", { type: exactResidue ? "ball_and_stick" : "cartoon" }, [
          mvsNode("color", { color: "#e11919" })
        ]),
        ...(exactResidue ? [mvsNode("label", { text: `${ranges[0].chain || "A"}:${ranges[0].start}` })] : [])
      ]));
    }
    return {
      metadata: {
        title: data.name || "Topology selection",
        version: "1",
        timestamp: new Date().toISOString()
      },
      root: mvsNode("root", {}, [
        mvsNode("download", { url }, [
          mvsNode("parse", { format: "mmcif" }, [
            mvsNode("structure", { type: "model" }, structureChildren)
          ])
        ]),
        mvsNode("canvas", { background_color: "#ffffff" })
      ])
    };
  }

  function molstarViewerUrlForMvs(ranges, selection, residue) {
    const mvs = buildMvsJson(ranges, selection, residue);
    if (!mvs) return "";
    return `https://molstar.org/viewer?mvs-format=mvsj&mvs-data=${encodeURIComponent(JSON.stringify(mvs))}`;
  }

  function loadMolstarIframeMvs(ranges, selection, residue) {
    if (!molstarNode) return false;
    const url = molstarViewerUrlForMvs(ranges, selection, residue);
    if (!url) return false;
    let frame = molstarNode.querySelector("iframe");
    if (!frame) {
      molstarNode.replaceChildren();
      frame = document.createElement("iframe");
      frame.title = "Mol* MolViewSpec view";
      frame.loading = "eager";
      molstarNode.appendChild(frame);
    }
    frame.src = url;
    if (molstarStatus) {
      const label = ranges.length === 1
        ? `${ranges[0].chain || "A"}:${ranges[0].start}-${ranges[0].stop}`
        : ranges.length ? `${ranges.length} selected ranges` : "structure";
      molstarStatus.textContent = ranges.length
        ? `MolViewSpec iframe loading red selection for ${label}.`
        : "MolViewSpec iframe loaded grey structure.";
    }
    return true;
  }

  function buildMolstarMvs(ranges, selection, residue) {
    const url = molstarModelUrl();
    if (!url || !mvsAvailable()) return null;
    const builder = window.molstar.PluginExtensions.mvs.MVSData.createBuilder();
    const structure = builder
      .download({ url })
      .parse({ format: "mmcif" })
      .modelStructure({});
    structure
      .component({ selector: "polymer" })
      .representation({ type: "cartoon" })
      .color({ color: "#b8bec7" });
    if (ranges.length) {
      const exactResidue = selection.mode === "residue" && residue;
      const selector = mvsSelectorForRanges(ranges, exactResidue);
      const active = structure.component({ selector });
      active.focus({});
      if (exactResidue) {
        active
          .representation({ type: "ball_and_stick" })
          .color({ color: "#e11919" });
        active.label({ text: `${ranges[0].chain || "A"}:${ranges[0].start}` });
      } else {
        active
          .representation({ type: "cartoon" })
          .color({ color: "#e11919" });
      }
    }
    return builder.getState();
  }

  async function loadMolstarMvs(ranges, selection, residue, options = {}) {
    if (!mvsAvailable()) return false;
    const mvsData = buildMolstarMvs(ranges, selection, residue);
    if (!mvsData) return false;
    if (molstarStatus) {
      const label = ranges.length === 1
        ? `${ranges[0].chain || "A"}:${ranges[0].start}-${ranges[0].stop}`
        : ranges.length ? `${ranges.length} selected ranges` : "structure";
      molstarStatus.textContent = `Loading MolViewSpec view for ${label}...`;
    }
    await window.molstar.PluginExtensions.mvs.loadMVS(molstarViewer.plugin, mvsData, {
      replaceExisting: options.replaceExisting ?? true,
      sanityChecks: true
    });
    configureMolstarAppearance();
    return true;
  }

  function molstarSchema(items) {
    return items.length === 1 ? items[0] : { items };
  }

  function molstarSchemasForRanges(ranges) {
    const mapped = ranges.map((range) => molstarRange(range.chain, range.start, range.stop));
    const schemas = [
      molstarSchema(mapped.map((range) => range.label)),
      molstarSchema(mapped.map((range) => range.auth))
    ];
    const chains = data.stats?.chains || [];
    if (chains.length <= 1) {
      schemas.push(molstarSchema(mapped.map((range) => range.labelAnyChain)));
      schemas.push(molstarSchema(mapped.map((range) => range.authAnyChain)));
    }
    return schemas;
  }

  function rangesForSelection(selection, residue) {
    if (selection.mode === "residue" && residue) {
      return [{ chain: residue.chain, start: residue.seq, stop: residue.seq }];
    }
    return selection.elements.map((element) => ({
      chain: element.chain,
      start: element.start,
      stop: element.stop
    }));
  }

  function focusOptionsForSelection(selection) {
    if (selection.mode === "residue") {
      return { minRadius: 24, extraRadius: 34, durationMs: 350 };
    }
    if (selection.mode === "sheet") {
      return { minRadius: 42, extraRadius: 34, durationMs: 350 };
    }
    return { minRadius: 34, extraRadius: 28, durationMs: 350 };
  }

  function clearMolstarSelection() {
    try {
      const plugin = molstarViewer?.plugin;
      molstarViewer?.structureInteractivity?.({ action: "select" });
      plugin?.managers?.interactivity?.lociSelects?.deselectAll?.();
      plugin?.managers?.structure?.focus?.clear?.();
      plugin?.managers?.interactivity?.clearHighlights?.();
    } catch (error) {}
  }

  function currentMolstarStructures() {
    const hierarchy = molstarViewer?.plugin?.managers?.structure?.hierarchy;
    const current = hierarchy?.current?.structures || hierarchy?.selection?.structures || [];
    const structures = [];
    current.forEach((item) => {
      const candidates = [
        item?.cell?.obj?.data,
        item?.components?.[0]?.cell?.obj?.data,
        item?.models?.[0]?.cell?.obj?.data
      ];
      candidates.forEach((candidate) => {
        if (candidate && candidate.units && !structures.includes(candidate)) structures.push(candidate);
      });
    });
    return structures;
  }

  function lociIsEmpty(loci) {
    const checker = window.molstar?.StructureElement?.Loci?.isEmpty;
    if (typeof checker === "function") return checker(loci);
    return !loci || !Array.isArray(loci.elements) || loci.elements.length === 0;
  }

  function lociFromSchema(structure, schema) {
    const structureElement = window.molstar?.StructureElement;
    const fromSchema = structureElement?.Loci?.fromSchema || structureElement?.Schema?.toLoci;
    if (typeof fromSchema !== "function") return null;
    try {
      return fromSchema(structure, schema);
    } catch (error) {
      return null;
    }
  }

  function molstarLociForRanges(ranges) {
    const structures = currentMolstarStructures();
    const loci = [];
    if (!structures.length) return loci;
    molstarSchemasForRanges(ranges).forEach((schema) => {
      structures.forEach((structure) => {
        const item = lociFromSchema(structure, schema);
        if (item && !lociIsEmpty(item)) loci.push(item);
      });
    });
    return loci;
  }

  function applyMolstarLoci(ranges, action, focusOptions) {
    const plugin = molstarViewer?.plugin;
    if (!plugin) return 0;
    const loci = molstarLociForRanges(ranges);
    if (!loci.length) return 0;
    loci.forEach((item, index) => {
      if (action === "highlight") {
        plugin.managers?.interactivity?.lociHighlights?.highlightOnly?.({ loci: item }, true);
      } else if (action === "select") {
        plugin.managers?.interactivity?.lociSelects?.select?.({ loci: item }, true);
        if (index === 0) plugin.managers?.structure?.focus?.setFromLoci?.(item);
      } else if (action === "focus") {
        if (index === 0) plugin.managers?.camera?.focusLoci?.(item, focusOptions);
      }
    });
    return loci.length;
  }

  function molstarExpressionForRanges(ranges, numbering, useChain) {
    return (Q) => {
      const groups = ranges.map((range) => {
        const seqProp = numbering === "auth"
          ? Q.struct.atomProperty.macromolecular.auth_seq_id()
          : Q.struct.atomProperty.macromolecular.label_seq_id();
        const group = {
          "residue-test": Number(range.start) === Number(range.stop)
            ? Q.core.rel.eq([seqProp, Number(range.start)])
            : Q.core.rel.inRange([seqProp, Number(range.start), Number(range.stop)])
        };
        if (useChain && range.chain) {
          const chainProp = numbering === "auth"
            ? Q.struct.atomProperty.macromolecular.auth_asym_id()
            : Q.struct.atomProperty.macromolecular.label_asym_id();
          group["chain-test"] = Q.core.rel.eq([chainProp, String(range.chain)]);
        }
        return Q.struct.generator.atomGroups(group);
      });
      return groups.length === 1 ? groups[0] : Q.struct.combinator.merge(groups);
    };
  }

  function invokeMolstarExpressions(ranges, action, focusOptions) {
    if (!molstarViewer || !molstarViewer.structureInteractivity) return 0;
    const singleChain = (data.stats?.chains || []).length <= 1;
    const attempts = [
      ["label", true],
      ["auth", true],
      ["label", false],
      ["auth", false]
    ];
    let sent = 0;
    attempts.forEach(([numbering, useChain]) => {
      if (!useChain && !singleChain) return;
      try {
        molstarViewer.structureInteractivity({
          expression: molstarExpressionForRanges(ranges, numbering, useChain),
          action,
          focusOptions,
          applyGranularity: true
        });
        sent += 1;
      } catch (error) {}
    });
    return sent;
  }

  function invokeMolstarSchemas(ranges, action, focusOptions) {
    if (!molstarViewer || !molstarViewer.structureInteractivity) return;
    let hits = invokeMolstarExpressions(ranges, action, focusOptions);
    hits += applyMolstarLoci(ranges, action, focusOptions);
    molstarSchemasForRanges(ranges).forEach((elements) => {
      molstarViewer.structureInteractivity({
        elements,
        action,
        focusOptions,
        applyGranularity: true
      });
      hits += 1;
    });
    return hits;
  }

  async function applyMolstarSelection(selection, residue) {
    const ranges = rangesForSelection(selection, residue);
    if (!ranges.length) return;
    try {
      if (loadMolstarIframeMvs(ranges, selection, residue)) return;
      if (!molstarViewer || !molstarViewer.structureInteractivity) return;
      if (await loadMolstarMvs(ranges, selection, residue, { replaceExisting: true })) {
        if (molstarStatus) {
          const target = ranges.length === 1
            ? `${ranges[0].chain || "A"}:${ranges[0].start}-${ranges[0].stop}`
            : `${ranges.length} ranges`;
          molstarStatus.textContent = `MolViewSpec highlighted ${target} in red.`;
        }
        return;
      }
      clearMolstarSelection();
      const selected = invokeMolstarSchemas(ranges, "select");
      invokeMolstarSchemas(ranges, "focus", focusOptionsForSelection(selection));
      lastMolstarHit = `${selection.mode}:${ranges.map((range) => `${range.chain}:${range.start}-${range.stop}`).join(",")}`;
      if (molstarStatus) {
        const target = ranges.length === 1
          ? `${ranges[0].chain || "chain"}:${ranges[0].start}-${ranges[0].stop}`
          : `${ranges.length} ranges`;
        molstarStatus.textContent = selected
          ? `Mol* residue/SSE selection sent for ${target}.`
          : `Mol* is loaded, but no selection command could be sent for ${target}.`;
      }
    } catch (error) {
      if (molstarStatus) molstarStatus.textContent = `Mol* selection failed: ${error.message || error}`;
    }
  }

  function highlightMolstarRange(chain, start, stop, action) {
    if (!molstarViewer || !molstarViewer.structureInteractivity) return;
    try {
      invokeMolstarSchemas([{ chain, start, stop }], action);
    } catch (error) {
      if (molstarStatus) molstarStatus.textContent = `Mol* link failed: ${error.message || error}`;
    }
  }

  function highlightMolstarElement(element, action) {
    if (!element) return;
    highlightMolstarRange(element.chain, element.start, element.stop, action);
  }

  function highlightMolstarSheet(element, action) {
    sheetElementsFor(element).forEach((item) => highlightMolstarElement(item, action));
  }

  function highlightMolstarResidue(residue, action) {
    if (!residue) return;
    highlightMolstarRange(residue.chain, residue.seq, residue.seq, action);
  }

  function firstOrderedSetValue(indices) {
    if (!indices) return null;
    if (typeof window.molstar?.OrderedSet?.getAt === "function") {
      const value = window.molstar.OrderedSet.getAt(indices, 0);
      if (Number.isFinite(Number(value))) return Number(value);
    }
    if (Array.isArray(indices) && Number.isFinite(Number(indices[0]))) return Number(indices[0]);
    if (typeof indices[Symbol.iterator] === "function") {
      for (const value of indices) {
        if (Number.isFinite(Number(value))) return Number(value);
        break;
      }
    }
    for (const key of ["array", "indices", "set"]) {
      const value = indices[key]?.[0];
      if (Number.isFinite(Number(value))) return Number(value);
    }
    return null;
  }

  function seqFromMolstarLoci(loci) {
    const elements = loci?.elements;
    if (!elements?.length) return null;
    for (const item of elements) {
      const unit = item?.unit;
      const atomic = unit?.model?.atomicHierarchy;
      const rawIndex = firstOrderedSetValue(item?.indices);
      if (!atomic || rawIndex == null) continue;
      const atomCandidates = [unit.elements?.[rawIndex], rawIndex].filter((value) => Number.isFinite(Number(value)));
      for (const atomIndex of atomCandidates) {
        const residueIndex = typeof unit?.residueIndex === "function"
          ? unit.residueIndex(atomIndex)
          : atomic.residueAtomSegments?.index?.[atomIndex];
        const authSeq = atomic.residues?.auth_seq_id?.value?.(residueIndex);
        const labelSeq = atomic.residues?.label_seq_id?.value?.(residueIndex);
        const seq = Number(authSeq ?? labelSeq);
        if (Number.isFinite(seq)) return seq;
      }
    }
    return null;
  }

  function seqFromMolstarStructure(structure) {
    const units = structure?.units || [];
    for (const unit of units) {
      const atomic = unit?.model?.atomicHierarchy;
      const elements = unit?.elements || [];
      if (!atomic || !elements.length) continue;
      for (const atomIndex of elements) {
        const residueIndex = atomic.residueAtomSegments?.index?.[atomIndex];
        const authSeq = atomic.residues?.auth_seq_id?.value?.(residueIndex);
        const labelSeq = atomic.residues?.label_seq_id?.value?.(residueIndex);
        const seq = Number(authSeq ?? labelSeq);
        if (Number.isFinite(seq)) return seq;
      }
    }
    return null;
  }

  function seqFromMolstarSelection() {
    const entries = molstarViewer?.plugin?.managers?.structure?.selection?.entries;
    if (!entries?.values) return null;
    for (const entry of entries.values()) {
      const seq = seqFromMolstarStructure(entry?.structure);
      if (seq != null) return seq;
    }
    return null;
  }

  function nestedMolstarValue(value, names) {
    const seen = new WeakSet();
    let checked = 0;
    function walk(item, depth) {
      if (!item || typeof item !== "object" || depth > 8 || checked > 1200) return null;
      if (seen.has(item) || ArrayBuffer.isView(item)) return null;
      seen.add(item);
      checked += 1;
      for (const name of names) {
        if (!Object.prototype.hasOwnProperty.call(item, name)) continue;
        const direct = item[name];
        if (typeof direct === "string" || typeof direct === "number") return direct;
        const nested = walk(direct, depth + 1);
        if (nested != null) return nested;
      }
      if (Array.isArray(item)) {
        for (const child of item.slice(0, 40)) {
          const found = walk(child, depth + 1);
          if (found != null) return found;
        }
        return null;
      }
      for (const key of Object.keys(item).slice(0, 80)) {
        const found = walk(item[key], depth + 1);
        if (found != null) return found;
      }
      return null;
    }
    return walk(value, 0);
  }

  function lociFromMolstarEvent(event) {
    return event?.current?.loci || event?.loci || event?.data?.current?.loci || event?.data?.loci || null;
  }

  function clearTopologyMolstarHover() {
    root.querySelectorAll(".molstar-hot").forEach((node) => node.classList.remove("molstar-hot"));
  }

  function syncTopologyFromMolstar(event, persistent) {
    const loci = lociFromMolstarEvent(event);
    const fallbackSeq = nestedMolstarValue(loci, ["auth_seq_id", "label_seq_id", "seq_id"]);
    const seq = seqFromMolstarLoci(loci)
      ?? (persistent ? seqFromMolstarSelection() : null)
      ?? (fallbackSeq == null ? null : Number(fallbackSeq));
    if (!Number.isFinite(Number(seq))) {
      if (!persistent) clearTopologyMolstarHover();
      return false;
    }
    const residue = residueBySeq.get(Number(seq));
    const element = elementForSeq(Number(seq));
    if (!element) return false;
    const group = nodeForElement(element.id);
    if (persistent) {
      selectElement(element.id, residue, { skipMolstar: true, mode: clickMode?.value || "sheet" });
      return true;
    }
    clearTopologyMolstarHover();
    if (group) {
      group.classList.add("molstar-hot");
      if (residue) updateTick(group, element, residue);
    }
    const plddt = residue?.plddt == null ? "not available" : Number(residue.plddt).toFixed(2);
    details.innerHTML = `
      <strong>Mol* hover: ${esc(element.display_id)} ${esc(element.ss_name)}</strong>
      <span> residue <code>${esc(residue?.residue_id || seq)}</code>, range <code>${esc(element.start)}-${esc(element.stop)}</code>, pLDDT <code>${esc(plddt)}</code>.</span>
    `;
    return true;
  }

  function bindMolstarToTopology() {
    const interaction = molstarViewer?.plugin?.behaviors?.interaction;
    let linked = false;
    try {
      if (interaction?.click?.subscribe) {
        interaction.click.subscribe((event) => syncTopologyFromMolstar(event, true));
        linked = true;
      }
      if (interaction?.hover?.subscribe) {
        interaction.hover.subscribe((event) => syncTopologyFromMolstar(event, false));
        linked = true;
      }
    } catch (error) {
      if (molstarStatus) molstarStatus.textContent = `Mol* reverse link failed: ${error.message || error}`;
    }
    return linked;
  }

  root.querySelector("[data-action='zoom-in']").addEventListener("click", () => zoomBy(1.2));
  root.querySelector("[data-action='zoom-out']").addEventListener("click", () => zoomBy(1 / 1.2));
  root.querySelector("[data-action='reset']").addEventListener("click", resetView);
  root.querySelector("[data-action='jump']").addEventListener("click", jumpToResidue);
  root.querySelector("[data-action='download-json']").addEventListener("click", downloadJson);
  residueInput.addEventListener("keydown", (event) => {
    if (event.key === "Enter") jumpToResidue();
  });
  linkToggle.addEventListener("change", () => {
    linksLayer.style.display = linkToggle.checked ? "" : "none";
  });
  labelToggle.addEventListener("change", () => {
    const display = labelToggle.checked ? "" : "none";
    root.querySelectorAll(".sse-label,.range-label,.sheet-text,.terminus,.domain-label").forEach((node) => {
      node.style.display = display;
    });
  });
  svg.addEventListener("wheel", (event) => {
    event.preventDefault();
    zoomBy(event.deltaY < 0 ? 1.08 : 1 / 1.08);
  }, { passive: false });
  svg.addEventListener("mousedown", (event) => {
    if (event.button !== 0) return;
    drag = { x: event.clientX, y: event.clientY };
    svg.classList.add("dragging");
  });
  window.addEventListener("mousemove", (event) => {
    if (!drag) return;
    const dx = event.clientX - drag.x;
    const dy = event.clientY - drag.y;
    transform.x += dx;
    transform.y += dy;
    drag.x = event.clientX;
    drag.y = event.clientY;
    applyTransform();
  });
  window.addEventListener("mouseup", () => {
    drag = null;
    svg.classList.remove("dragging");
  });

  layout();
  drawBands();
  drawConnectors();
  drawLinks();
  drawElements();
  applyTransform();
  initMolstar();
})();
""".replace(
        "__ROOT_ID__", root_id
    ).replace(
        "__DATA_ID__", data_id
    )

    return f"""
<div id="{root_id}">
  <style>{css}</style>
  <div class="clean-shell">
    <div class="clean-header">
      <div class="title">{_html_escape(title)}</div>
      <div class="subtitle">{_html_escape(subtitle)}</div>
      <div class="stats">
        {stats["residue_count"]} residues;
        {stats["helix_count"]} helices;
        {stats["strand_count"]} strands;
        {stats["beta_link_count"]} inferred strand contacts;
        chain {_html_escape(", ".join(stats["chains"]))}
      </div>
    </div>
    <div class="clean-grid">
      <section class="topology-panel">
        <div class="toolbar">
          <div class="label">{_html_escape(title)} | Clean topology</div>
          <button type="button" data-action="zoom-in" title="Zoom in">+</button>
          <button type="button" data-action="zoom-out" title="Zoom out">-</button>
          <button type="button" data-action="reset" title="Reset pan and zoom">Fit</button>
          <label class="toggle" title="Show inferred strand contacts">
            <input type="checkbox" data-role="toggle-links"> Contacts
          </label>
          <label class="toggle" title="Show labels">
            <input type="checkbox" data-role="toggle-labels" checked> Labels
          </label>
          <select data-role="click-mode" title="Choose what a topology click selects in Mol*">
            <option value="residue" selected>Click: residue</option>
            <option value="range">Click: SSE range</option>
            <option value="sheet">Click: helix/sheet</option>
          </select>
          <input type="text" data-role="residue-input" placeholder="Residue A:120">
          <button type="button" data-action="jump" title="Select a residue">Find</button>
          <button type="button" data-action="download-json" title="Download topology JSON">JSON</button>
        </div>
        <div class="viewer-wrap">
          <svg role="img" aria-label="Clean AlphaFold secondary-structure topology">
            <g data-role="viewport">
              <g data-role="bands"></g>
              <g data-role="connectors"></g>
              <g data-role="links"></g>
              <g data-role="elements"></g>
            </g>
          </svg>
          <div class="tooltip" data-role="tooltip"></div>
        </div>
        <div class="details" data-role="details">
          Helices are labeled H1, H2, ... and strands are labeled E1, E2, ... . Click a position on an element to highlight that residue in Mol* with contextual zoom; change click mode to select an SSE or sheet.
        </div>
      </section>
      {molstar_panel}
    </div>
  </div>
  <script type="application/json" id="{data_id}">{json_blob}</script>
  <script>{script}</script>
</div>
"""


def _pdbe_plugin_topology_html(topology: Dict[str, Any]) -> str:
    root_id = "pdbe-topology-" + uuid.uuid4().hex
    data_id = root_id + "-data"
    title = topology.get("name") or topology.get("pdbe_entry_id") or "Topology"
    stats = topology["stats"]
    metadata = topology.get("metadata", {})
    subtitle = metadata.get("source") or "Generated topology"
    json_blob = json.dumps(topology, separators=(",", ":")).replace("</", "<\\/")
    molstar_panel = ""
    if topology.get("cif_url") or topology.get("afdb_accession"):
        molstar_panel = """
        <section class="af-molstar-panel">
          <div class="af-panel-bar">Mol* 3D view</div>
          <div class="af-molstar-stage" data-role="molstar"></div>
          <div class="af-status" data-role="molstar-status">Loading Mol*...</div>
        </section>
        """

    css = """
#__ROOT_ID__ {
  --line: #c8d5c0;
  --header: #dcefd4;
  --ink: #1f2933;
  color: var(--ink);
  font-family: Inter, "Segoe UI", Arial, sans-serif;
}
#__ROOT_ID__ .af-shell {
  border: 1px solid var(--line);
  background: white;
}
#__ROOT_ID__ .af-header {
  background: var(--header);
  border-bottom: 1px solid var(--line);
  padding: 12px 14px;
}
#__ROOT_ID__ .af-title {
  color: #2f6f39;
  font-size: 21px;
  font-weight: 700;
}
#__ROOT_ID__ .af-subtitle,
#__ROOT_ID__ .af-stats {
  color: #3e5945;
  font-size: 12px;
  margin-top: 3px;
}
#__ROOT_ID__ .af-grid {
  display: grid;
  gap: 12px;
  grid-template-columns: minmax(420px, 1.15fr) minmax(360px, 0.85fr);
  padding: 12px;
}
#__ROOT_ID__ .pdbe-topology-target {
  height: 650px;
  min-width: 0;
}
#__ROOT_ID__ .af-molstar-panel {
  border: 1px solid #b9c5d1;
  min-width: 0;
  position: relative;
}
#__ROOT_ID__ .af-panel-bar {
  background: #f0f3f6;
  border-bottom: 1px solid #c9d1da;
  font-size: 13px;
  min-height: 34px;
  padding: 8px 10px;
}
#__ROOT_ID__ .af-molstar-stage {
  height: 650px;
  position: relative;
}
#__ROOT_ID__ .af-status {
  background: rgba(255, 255, 255, 0.92);
  bottom: 8px;
  color: #344154;
  font-size: 12px;
  left: 8px;
  padding: 4px 6px;
  position: absolute;
}
#__ROOT_ID__ .af-details {
  border-top: 1px solid #d7d7d7;
  color: #5d6978;
  font-size: 13px;
  min-height: 42px;
  padding: 10px 14px;
}
#__ROOT_ID__ .af-details strong {
  color: var(--ink);
}
#__ROOT_ID__ .af-details code {
  background: #eef2f7;
  border-radius: 3px;
  color: #344154;
  padding: 1px 4px;
}
@media (max-width: 1050px) {
  #__ROOT_ID__ .af-grid {
    grid-template-columns: 1fr;
  }
}
""".replace(
        "__ROOT_ID__", root_id
    )

    script = r"""
(function () {
  const root = document.getElementById("__ROOT_ID__");
  const data = JSON.parse(document.getElementById("__DATA_ID__").textContent);
  const target = root.querySelector("[data-role='pdbe-topology']");
  const details = root.querySelector("[data-role='details']");
  const molstarNode = root.querySelector("[data-role='molstar']");
  const molstarStatus = root.querySelector("[data-role='molstar-status']");
  const residues = new Map(data.residues.map((residue) => [String(residue.seq), residue]));
  let molstarViewer = null;

  function esc(value) {
    return String(value ?? "")
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#039;");
  }

  function loadScriptOnce(url, id) {
    if (document.getElementById(id)?.dataset.ready === "true") return Promise.resolve();
    if (id === "pdbe-topology-plugin" && window.PdbTopologyViewerPlugin) return Promise.resolve();
    if (id === "molstar-viewer-js" && window.molstar) return Promise.resolve();
    const existing = document.getElementById(id);
    if (existing) {
      return new Promise((resolve, reject) => {
        existing.addEventListener("load", resolve, { once: true });
        existing.addEventListener("error", reject, { once: true });
      });
    }
    return new Promise((resolve, reject) => {
      const script = document.createElement("script");
      script.id = id;
      script.src = url;
      script.onload = () => {
        script.dataset.ready = "true";
        resolve();
      };
      script.onerror = reject;
      document.head.appendChild(script);
    });
  }

  function loadCssOnce(url, id) {
    if (document.getElementById(id)) return;
    const link = document.createElement("link");
    link.id = id;
    link.rel = "stylesheet";
    link.href = url;
    document.head.appendChild(link);
  }

  function renderDetails(eventData) {
    const residue = residues.get(String(eventData.residueNumber));
    if (!residue) return;
    const plddt = residue.plddt == null ? "not available" : Number(residue.plddt).toFixed(2);
    const element = data.elements.find((item) => residue.seq >= item.start && residue.seq <= item.stop);
    const contactCount = data.links.filter((link) => {
      if (!element) return false;
      return link.source === element.id || link.target === element.id;
    }).length;
    details.innerHTML = `
      <strong>${esc(residue.residue_id)} ${esc(residue.aa)}</strong>
      <span> pLDDT <code>${esc(plddt)}</code>, SSE <code>${esc(element ? element.id + " " + element.ss_name : "coil")}</code>,
      strand contacts <code>${contactCount}</code>.</span>
    `;
  }

  function highlightMolstar(eventData, action) {
    if (!molstarViewer || !molstarViewer.structureInteractivity) return;
    const residue = residues.get(String(eventData.residueNumber));
    if (!residue) return;
    try {
      molstarViewer.structureInteractivity({
        elements: {
          chain_id: residue.chain,
          auth_asym_id: residue.chain,
          label_asym_id: residue.chain,
          beg_auth_seq_id: residue.seq,
          end_auth_seq_id: residue.seq,
          beg_label_seq_id: residue.seq,
          end_label_seq_id: residue.seq
        },
        action
      });
    } catch (error) {
      if (molstarStatus) molstarStatus.textContent = `Mol* residue link failed: ${error.message || error}`;
    }
  }

  function handleTopologyEvent(event, action) {
    const eventData = event.eventData;
    if (!eventData) return;
    if (String(eventData.entryId).toLowerCase() !== data.pdbe_entry_id) return;
    if (String(eventData.entityId) !== data.pdbe_entity_id) return;
    if (String(eventData.chainId) !== data.pdbe_chain_id) return;
    renderDetails(eventData);
    highlightMolstar(eventData, action);
  }

  async function initTopology() {
    try {
      await loadScriptOnce("https://cdn.jsdelivr.net/gh/PDBeurope/pdb-topology-viewer@master/build/pdb-topology-viewer-plugin-3.0.1.js", "pdbe-topology-plugin");
      const Plugin = window.PdbTopologyViewerPlugin;
      if (!Plugin) throw new Error("PdbTopologyViewerPlugin was not exposed by the loaded bundle.");
      const viewer = new Plugin();
      viewer.render(target, {
        entryId: data.pdbe_entry_id,
        entityId: data.pdbe_entity_id,
        chainId: data.pdbe_chain_id,
        subscribeEvents: true,
        autoResize: true,
        displayStyle: "height:650px;width:100%;"
      }, data.pdbe_api_data);
    } catch (error) {
      target.innerHTML = `<div style="padding:16px;color:#b42318;">PDBe topology renderer could not load: ${esc(error.message || error)}</div>`;
    }
  }

  async function initMolstar() {
    if (!molstarNode) return;
    try {
      loadCssOnce("https://cdn.jsdelivr.net/npm/molstar@5.4.2/build/viewer/molstar.css", "molstar-viewer-css");
      await loadScriptOnce("https://cdn.jsdelivr.net/npm/molstar@5.4.2/build/viewer/molstar.js", "molstar-viewer-js");
      molstarViewer = await window.molstar.Viewer.create(molstarNode, {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: true,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        viewportBackgroundColor: "white"
      });
      if (data.cif_url) {
        await molstarViewer.loadStructureFromUrl(data.cif_url, "mmcif", false);
      } else if (data.afdb_accession && molstarViewer.loadAlphaFoldDb) {
        await molstarViewer.loadAlphaFoldDb(`AF-${data.afdb_accession}-F1`);
      }
      if (molstarStatus) molstarStatus.textContent = "Hover or click topology residues to link 2D and 3D.";
    } catch (error) {
      if (molstarStatus) molstarStatus.textContent = `Mol* could not load: ${error.message || error}`;
    }
  }

  document.addEventListener("PDB.topologyViewer.mouseover", (event) => handleTopologyEvent(event, "highlight"));
  document.addEventListener("PDB.topologyViewer.click", (event) => handleTopologyEvent(event, "select"));
  document.addEventListener("PDB.topologyViewer.mouseout", () => {
    try { molstarViewer?.plugin?.managers?.interactivity?.clearHighlights?.(); } catch (error) {}
  });
  initTopology();
  initMolstar();
})();
""".replace(
        "__ROOT_ID__", root_id
    ).replace(
        "__DATA_ID__", data_id
    )

    return f"""
<div id="{root_id}">
  <style>{css}</style>
  <div class="af-shell">
    <div class="af-header">
      <div class="af-title">{_html_escape(title)}</div>
      <div class="af-subtitle">{_html_escape(subtitle)}</div>
      <div class="af-stats">
        {stats["residue_count"]} residues;
        {stats["helix_count"]} helices;
        {stats["strand_count"]} strands;
        {stats["beta_link_count"]} strand contacts;
        chain {_html_escape(topology.get("pdbe_chain_id", ""))}
      </div>
    </div>
    <div class="af-grid">
      <section>
        <div class="pdbe-topology-target" data-role="pdbe-topology"></div>
        <div class="af-details" data-role="details">
          Official PDBe topology renderer, using topology JSON generated from the AlphaFold mmCIF. Hover residues for pLDDT and secondary-structure details.
        </div>
      </section>
      {molstar_panel}
    </div>
  </div>
  <script type="application/json" id="{data_id}">{json_blob}</script>
  <script>{script}</script>
</div>
"""


def topology_html(topology: Dict[str, Any]) -> str:
    if topology.get("pdbe_api_data"):
        return _clean_alphafold_topology_html(topology)

    root_id = "dssp-topology-" + uuid.uuid4().hex
    data_id = root_id + "-data"
    title = topology.get("name") or "DSSP topology"
    metadata = topology.get("metadata", {})
    stats = topology["stats"]
    subtitle = metadata.get("molecule") or metadata.get("header") or "Uploaded DSSP"
    json_blob = json.dumps(topology, separators=(",", ":")).replace("</", "<\\/")
    afdb_accession = str(topology.get("afdb_accession") or "").strip()
    afdb_for_url = re.sub(r"[^A-Za-z0-9_.-]", "", afdb_accession)
    molstar_url = (
        f"https://molstar.org/viewer/?afdb={afdb_for_url}&hide-controls=1"
        if afdb_for_url
        else ""
    )
    molstar_panel = ""
    if afdb_accession:
        molstar_panel = f"""
        <section class="molstar-panel" aria-label="AlphaFold Molstar viewer">
          <div class="molstar-topbar">
            <span>AlphaFold DB: <strong>{_html_escape(afdb_accession)}</strong></span>
            <a data-role="molstar-link" href="{_html_escape(molstar_url)}" target="_blank" rel="noopener">Open Mol*</a>
          </div>
          <div class="molstar-stage" data-role="molstar"></div>
          <div class="molstar-status" data-role="molstar-status">Loading Mol* from AlphaFold DB...</div>
        </section>
        """

    css = """
#__ROOT_ID__ {
  --ink: #111111;
  --muted: #5f6670;
  --line: #222222;
  --panel: #f5f7f9;
  --pdbe-green: #dcefd4;
  --pdbe-green-border: #b9d9ad;
  --bar: #707070;
  --select: #1f5fbf;
  color: var(--ink);
  font-family: Inter, "Segoe UI", Arial, sans-serif;
  line-height: 1.35;
}
#__ROOT_ID__ .topology-shell {
  border: 1px solid #c9d4c1;
  background: white;
}
#__ROOT_ID__ .topology-heading {
  background: var(--pdbe-green);
  border-bottom: 1px solid var(--pdbe-green-border);
  padding: 13px 16px;
}
#__ROOT_ID__ .topology-heading .title {
  color: #2f6f39;
  font-size: 22px;
  font-weight: 700;
}
#__ROOT_ID__ .topology-heading .subtitle,
#__ROOT_ID__ .topology-heading .stats {
  color: #3e5945;
  font-size: 12px;
  margin-top: 3px;
}
#__ROOT_ID__ .viewer-grid {
  display: grid;
  gap: 12px;
  grid-template-columns: minmax(360px, 1fr);
  padding: 14px;
}
#__ROOT_ID__.has-molstar .viewer-grid {
  grid-template-columns: minmax(420px, 1.25fr) minmax(360px, 0.75fr);
}
#__ROOT_ID__ .topology-card {
  border: 2px solid #777777;
  min-width: 0;
}
#__ROOT_ID__ .viewer-wrap {
  height: 640px;
  position: relative;
  background: #ffffff;
}
#__ROOT_ID__ svg {
  cursor: grab;
  display: block;
  height: 100%;
  width: 100%;
}
#__ROOT_ID__ svg.is-dragging {
  cursor: grabbing;
}
#__ROOT_ID__ .connector {
  fill: none;
  stroke: var(--line);
  stroke-linecap: square;
  stroke-linejoin: miter;
  stroke-width: 2;
}
#__ROOT_ID__ .chain-break {
  stroke-dasharray: 5 6;
}
#__ROOT_ID__ .beta-link {
  fill: none;
  opacity: 0.22;
  stroke: #4b79bd;
  stroke-linecap: round;
}
#__ROOT_ID__ .beta-link.link-selected {
  opacity: 0.9;
  stroke: #1f4e9d;
}
#__ROOT_ID__ .sse {
  cursor: pointer;
  outline: none;
}
#__ROOT_ID__ .sse-shape {
  fill: #ffffff;
  stroke: #111111;
  stroke-linejoin: miter;
  stroke-width: 2;
}
#__ROOT_ID__ .sse:hover .sse-shape,
#__ROOT_ID__ .sse.selected .sse-shape {
  stroke: var(--select);
  stroke-width: 3;
}
#__ROOT_ID__ .selection-ring {
  fill: none;
  opacity: 0;
  pointer-events: none;
  stroke: var(--select);
  stroke-dasharray: 4 4;
  stroke-width: 2;
}
#__ROOT_ID__ .sse.selected .selection-ring {
  opacity: 1;
}
#__ROOT_ID__ .residue-tick {
  display: none;
  pointer-events: none;
  stroke: #d01f1f;
  stroke-width: 2;
}
#__ROOT_ID__ .sse.has-residue .residue-tick {
  display: block;
}
#__ROOT_ID__ .element-label {
  fill: #111111;
  font-size: 12px;
  font-weight: 700;
  pointer-events: none;
  text-anchor: middle;
}
#__ROOT_ID__ .residue-label {
  fill: #555555;
  font-size: 10px;
  pointer-events: none;
  text-anchor: middle;
}
#__ROOT_ID__ .terminus {
  font-size: 18px;
  font-weight: 700;
  pointer-events: none;
}
#__ROOT_ID__ .terminus.n {
  fill: #0b39ff;
}
#__ROOT_ID__ .terminus.c {
  fill: #ff0000;
}
#__ROOT_ID__ .tooltip {
  background: #111827;
  border-radius: 4px;
  color: white;
  display: none;
  font-size: 12px;
  max-width: 320px;
  padding: 8px 10px;
  pointer-events: none;
  position: absolute;
  z-index: 4;
}
#__ROOT_ID__ .pdbe-footer {
  align-items: stretch;
  background: var(--bar);
  color: white;
  display: grid;
  grid-template-columns: minmax(220px, 1fr) auto;
  min-height: 48px;
}
#__ROOT_ID__ .entry-label {
  align-items: center;
  display: flex;
  font-size: 20px;
  gap: 8px;
  min-width: 0;
  padding: 8px 16px;
}
#__ROOT_ID__ .entry-label span {
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}
#__ROOT_ID__ .entry-dot {
  border: 2px solid #555555;
  box-shadow: inset 0 0 0 4px #f1f1f1;
  display: inline-block;
  height: 18px;
  width: 18px;
}
#__ROOT_ID__ .controls {
  align-items: center;
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  justify-content: flex-end;
  padding: 6px 8px;
}
#__ROOT_ID__ button,
#__ROOT_ID__ input[type="text"],
#__ROOT_ID__ select {
  border: 1px solid #c6c6c6;
  border-radius: 0;
  background: #ffffff;
  color: #111111;
  font: inherit;
  font-size: 12px;
  min-height: 30px;
}
#__ROOT_ID__ button {
  cursor: pointer;
  min-width: 30px;
  padding: 4px 8px;
}
#__ROOT_ID__ button:hover {
  background: #eeeeee;
}
#__ROOT_ID__ input[type="text"] {
  width: 128px;
  padding: 4px 8px;
}
#__ROOT_ID__ select {
  min-width: 120px;
  padding: 4px 8px;
}
#__ROOT_ID__ label.toggle {
  align-items: center;
  color: white;
  display: inline-flex;
  font-size: 12px;
  gap: 4px;
  min-height: 30px;
  white-space: nowrap;
}
#__ROOT_ID__ .details {
  border-top: 1px solid #d7d7d7;
  color: var(--muted);
  font-size: 13px;
  min-height: 48px;
  padding: 10px 14px;
}
#__ROOT_ID__ .details strong {
  color: var(--ink);
}
#__ROOT_ID__ .details code {
  background: #eef2f7;
  border-radius: 3px;
  color: #344154;
  padding: 1px 4px;
}
#__ROOT_ID__ .molstar-panel {
  border: 1px solid #b9c5d1;
  min-width: 0;
  position: relative;
}
#__ROOT_ID__ .molstar-topbar {
  align-items: center;
  background: #f0f3f6;
  border-bottom: 1px solid #c9d1da;
  color: #1f2933;
  display: flex;
  font-size: 13px;
  justify-content: space-between;
  min-height: 38px;
  padding: 6px 10px;
}
#__ROOT_ID__ .molstar-topbar a {
  color: #164f9f;
  text-decoration: none;
}
#__ROOT_ID__ .molstar-stage {
  height: 640px;
  position: relative;
}
#__ROOT_ID__ .molstar-status {
  background: rgba(255, 255, 255, 0.92);
  bottom: 8px;
  color: #344154;
  font-size: 12px;
  left: 8px;
  padding: 4px 6px;
  position: absolute;
}
@media (max-width: 1050px) {
  #__ROOT_ID__.has-molstar .viewer-grid {
    grid-template-columns: 1fr;
  }
}
@media (max-width: 760px) {
  #__ROOT_ID__ .viewer-grid {
    padding: 10px;
  }
  #__ROOT_ID__ .viewer-wrap,
  #__ROOT_ID__ .molstar-stage {
    height: 520px;
  }
  #__ROOT_ID__ .pdbe-footer {
    grid-template-columns: 1fr;
  }
  #__ROOT_ID__ .controls {
    justify-content: flex-start;
  }
}
""".replace(
        "__ROOT_ID__", root_id
    )

    script = r"""
(function () {
  const root = document.getElementById("__ROOT_ID__");
  const data = JSON.parse(document.getElementById("__DATA_ID__").textContent);
  const svg = root.querySelector("svg");
  const viewport = root.querySelector("[data-role='viewport']");
  const connectorsLayer = root.querySelector("[data-role='connectors']");
  const linksLayer = root.querySelector("[data-role='links']");
  const elementsLayer = root.querySelector("[data-role='elements']");
  const tooltip = root.querySelector("[data-role='tooltip']");
  const details = root.querySelector("[data-role='details']");
  const linkToggle = root.querySelector("[data-role='toggle-links']");
  const labelToggle = root.querySelector("[data-role='toggle-labels']");
  const residueInput = root.querySelector("[data-role='residue-input']");
  const molstarNode = root.querySelector("[data-role='molstar']");
  const molstarStatus = root.querySelector("[data-role='molstar-status']");
  const elementById = new Map(data.elements.map((element) => [element.id, element]));
  const residueByDssp = new Map(data.residues.map((residue) => [residue.dssp_index, residue]));
  const elementByResidue = new Map();
  const residueByLookup = new Map();
  let selectedId = null;
  let transform = { x: 0, y: 0, k: 1 };
  let drag = null;
  let molstarViewer = null;
  let lastMolstarResidue = "";

  data.elements.forEach((element) => {
    element.residue_ids.forEach((residueId, index) => {
      const residue = residueByDssp.get(element.residue_indices[index]);
      const keys = [
        residueId,
        residueId.split(":").pop(),
        residue ? `${residue.chain}:${residue.residue_number}${residue.insertion_code || ""}` : "",
        residue ? `${residue.residue_number}${residue.insertion_code || ""}` : ""
      ].filter(Boolean);
      keys.forEach((key) => {
        elementByResidue.set(key.toLowerCase(), element.id);
        if (residue) residueByLookup.set(key.toLowerCase(), residue);
      });
    });
  });

  function esc(value) {
    return String(value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#039;");
  }

  function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
  }

  function makeSvg(tag, attrs) {
    const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
    Object.entries(attrs || {}).forEach(([key, value]) => {
      if (value === null || value === undefined) return;
      if (key === "textContent") node.textContent = value;
      else node.setAttribute(key, value);
    });
    return node;
  }

  function residuesFor(element) {
    return element.residue_indices
      .map((dsspIndex) => residueByDssp.get(dsspIndex))
      .filter(Boolean);
  }

  function layoutElements() {
    const count = Math.max(1, data.elements.length);
    const columns = Math.min(14, Math.max(6, Math.ceil(Math.sqrt(count * 1.8))));
    const left = 58;
    const top = 52;
    const colGap = 78;
    const rowGap = 208;
    const diagramWidth = left * 2 + (columns - 1) * colGap + 70;
    let maxBottom = top;

    data.elements.forEach((element, index) => {
      const row = Math.floor(index / columns);
      const rawCol = index % columns;
      const col = row % 2 === 0 ? rawCol : columns - rawCol - 1;
      const stagger = ((index + row) % 3) * 12;
      const lengthHeight = 52 + element.length * 3.2;
      element.w = element.type === "helix" ? 38 : 42;
      element.h = Math.max(74, Math.min(168, lengthHeight));
      element.x = left + col * colGap;
      element.y = top + row * rowGap + stagger;
      element.cx = element.x + element.w / 2;
      element.cy = element.y + element.h / 2;
      element.direction = index % 2 === 0 ? 1 : -1;
      maxBottom = Math.max(maxBottom, element.y + element.h + 52);
    });

    svg.setAttribute("viewBox", `0 0 ${diagramWidth} ${maxBottom + 36}`);
  }

  function flowStart(element) {
    return element.direction > 0
      ? { x: element.cx, y: element.y }
      : { x: element.cx, y: element.y + element.h };
  }

  function flowEnd(element) {
    return element.direction > 0
      ? { x: element.cx, y: element.y + element.h }
      : { x: element.cx, y: element.y };
  }

  function connectorPath(a, b) {
    const start = flowEnd(a);
    const end = flowStart(b);
    const midY = Math.abs(start.y - end.y) < 16
      ? start.y + (start.y <= end.y ? 26 : -26)
      : (start.y + end.y) / 2;
    return `M ${start.x} ${start.y} L ${start.x} ${midY} L ${end.x} ${midY} L ${end.x} ${end.y}`;
  }

  function linkPath(a, b) {
    const lift = Math.abs(a.cy - b.cy) < 24 ? -58 : 0;
    const midY = (a.cy + b.cy) / 2 + lift;
    return `M ${a.cx} ${a.cy} C ${a.cx} ${midY}, ${b.cx} ${midY}, ${b.cx} ${b.cy}`;
  }

  function drawConnectors() {
    connectorsLayer.replaceChildren();
    for (let i = 0; i < data.elements.length - 1; i += 1) {
      const a = data.elements[i];
      const b = data.elements[i + 1];
      const path = makeSvg("path", {
        class: a.chain === b.chain ? "connector" : "connector chain-break",
        d: connectorPath(a, b)
      });
      connectorsLayer.appendChild(path);
    }
    drawTermini();
  }

  function drawTermini() {
    if (!data.elements.length) return;
    const first = data.elements[0];
    const last = data.elements[data.elements.length - 1];
    const start = flowStart(first);
    const end = flowEnd(last);
    connectorsLayer.appendChild(makeSvg("text", {
      class: "terminus n",
      x: start.x - 18,
      y: start.y + 6,
      textContent: "N"
    }));
    connectorsLayer.appendChild(makeSvg("text", {
      class: "terminus c",
      x: end.x + 12,
      y: end.y + 6,
      textContent: "C"
    }));
  }

  function drawLinks() {
    linksLayer.replaceChildren();
    data.links.forEach((link) => {
      const a = elementById.get(link.source);
      const b = elementById.get(link.target);
      if (!a || !b) return;
      const path = makeSvg("path", {
        class: "beta-link",
        d: linkPath(a, b),
        "data-link": `${link.source}:${link.target}`,
        "data-source": link.source,
        "data-target": link.target,
        "stroke-width": Math.min(7, 1.2 + Math.sqrt(link.count))
      });
      path.appendChild(makeSvg("title", {
        textContent: `${link.source} to ${link.target}: ${link.count} DSSP bridge contacts`
      }));
      linksLayer.appendChild(path);
    });
  }

  function drawHelix(group, element) {
    group.appendChild(makeSvg("rect", {
      class: "sse-shape",
      x: element.x,
      y: element.y,
      width: element.w,
      height: element.h,
      rx: element.w / 2,
      ry: element.w / 2
    }));
    group.appendChild(makeSvg("rect", {
      class: "selection-ring",
      x: element.x - 6,
      y: element.y - 6,
      width: element.w + 12,
      height: element.h + 12,
      rx: element.w / 2 + 6,
      ry: element.w / 2 + 6
    }));
  }

  function strandPoints(element) {
    const head = Math.min(30, Math.max(22, element.h * 0.22));
    const x = element.x;
    const y = element.y;
    const w = element.w;
    const h = element.h;
    if (element.direction > 0) {
      return [
        `${x},${y}`,
        `${x + w},${y}`,
        `${x + w},${y + h - head}`,
        `${x + w * 0.68},${y + h - head}`,
        `${x + w / 2},${y + h}`,
        `${x + w * 0.32},${y + h - head}`,
        `${x},${y + h - head}`
      ].join(" ");
    }
    return [
      `${x + w / 2},${y}`,
      `${x + w * 0.68},${y + head}`,
      `${x + w},${y + head}`,
      `${x + w},${y + h}`,
      `${x},${y + h}`,
      `${x},${y + head}`,
      `${x + w * 0.32},${y + head}`
    ].join(" ");
  }

  function drawStrand(group, element) {
    group.appendChild(makeSvg("polygon", {
      class: "sse-shape",
      points: strandPoints(element)
    }));
    group.appendChild(makeSvg("rect", {
      class: "selection-ring",
      x: element.x - 6,
      y: element.y - 6,
      width: element.w + 12,
      height: element.h + 12,
      rx: 3
    }));
  }

  function drawElements() {
    elementsLayer.replaceChildren();
    data.elements.forEach((element) => {
      const group = makeSvg("g", {
        class: "sse",
        "data-element-id": element.id,
        tabindex: 0
      });
      if (element.type === "helix") drawHelix(group, element);
      else drawStrand(group, element);

      group.appendChild(makeSvg("line", {
        class: "residue-tick",
        x1: element.x - 7,
        x2: element.x + element.w + 7,
        y1: element.cy,
        y2: element.cy
      }));
      group.appendChild(makeSvg("text", {
        class: "element-label",
        x: element.cx,
        y: element.y - 12,
        textContent: element.id
      }));
      group.appendChild(makeSvg("text", {
        class: "residue-label",
        x: element.cx,
        y: element.y + element.h + 16,
        textContent: `${element.start_residue}-${element.end_residue}`
      }));

      group.addEventListener("mouseenter", (event) => showResidueTooltip(event, element, group));
      group.addEventListener("mousemove", (event) => showResidueTooltip(event, element, group));
      group.addEventListener("mouseleave", () => {
        group.classList.remove("has-residue");
        hideTooltip();
        clearMolstarHover();
      });
      group.addEventListener("click", (event) => {
        event.stopPropagation();
        selectElement(element.id, residueAtEvent(event, element));
      });
      group.addEventListener("keydown", (event) => {
        if (event.key === "Enter" || event.key === " ") selectElement(element.id, residuesFor(element)[0]);
      });
      elementsLayer.appendChild(group);
    });
  }

  function residueAtEvent(event, element) {
    const residues = residuesFor(element);
    if (!residues.length) return null;
    const point = svg.createSVGPoint();
    point.x = event.clientX;
    point.y = event.clientY;
    const matrix = viewport.getScreenCTM();
    if (!matrix) return residues[0];
    const local = point.matrixTransform(matrix.inverse());
    const t = clamp((local.y - element.y) / element.h, 0, 1);
    const raw = element.direction > 0 ? t : 1 - t;
    const index = Math.round(raw * (residues.length - 1));
    return residues[clamp(index, 0, residues.length - 1)];
  }

  function updateResidueTick(group, element, residue) {
    if (!residue) {
      group.classList.remove("has-residue");
      return;
    }
    const residues = residuesFor(element);
    const index = residues.findIndex((item) => item.dssp_index === residue.dssp_index);
    if (index < 0) return;
    const t = residues.length <= 1 ? 0.5 : index / (residues.length - 1);
    const visualT = element.direction > 0 ? t : 1 - t;
    const y = element.y + visualT * element.h;
    const tick = group.querySelector(".residue-tick");
    tick.setAttribute("y1", y);
    tick.setAttribute("y2", y);
    group.classList.add("has-residue");
  }

  function showResidueTooltip(event, element, group) {
    const residue = residueAtEvent(event, element);
    updateResidueTick(group, element, residue);
    if (!residue) return;
    const bridges = [residue.bp1, residue.bp2].filter(Boolean).join(", ") || "none";
    const acc = residue.accessibility === null || residue.accessibility === undefined
      ? "not available"
      : residue.accessibility;
    const coords = residue.x === null || residue.x === undefined
      ? "not available"
      : `${Number(residue.x).toFixed(2)}, ${Number(residue.y).toFixed(2)}, ${Number(residue.z).toFixed(2)}`;
    tooltip.innerHTML = `
      <strong>${esc(residue.residue_id)} ${esc(residue.aa)}</strong><br>
      ${esc(residue.ss_name)} (${esc(residue.ss)}) in ${esc(element.id)}<br>
      DSSP row ${esc(residue.dssp_index)}; ACC ${esc(acc)}; BP ${esc(bridges)}<br>
      XYZ ${esc(coords)}
    `;
    tooltip.style.display = "block";
    moveTooltip(event);
    highlightMolstarResidue(residue, "highlight");
  }

  function moveTooltip(event) {
    const box = root.getBoundingClientRect();
    tooltip.style.left = `${event.clientX - box.left + 14}px`;
    tooltip.style.top = `${event.clientY - box.top + 14}px`;
  }

  function hideTooltip() {
    tooltip.style.display = "none";
  }

  function selectElement(elementId, residue) {
    selectedId = elementId;
    const element = elementById.get(elementId);
    root.querySelectorAll("[data-element-id]").forEach((node) => {
      node.classList.toggle("selected", node.getAttribute("data-element-id") === elementId);
    });
    root.querySelectorAll("[data-link]").forEach((node) => {
      const linked = node.getAttribute("data-source") === elementId || node.getAttribute("data-target") === elementId;
      node.classList.toggle("link-selected", linked);
    });
    if (!element) return;
    const acc = element.accessibility_mean === null || element.accessibility_mean === undefined
      ? "not available"
      : element.accessibility_mean;
    const residueLine = residue
      ? `<br><span>Residue: <code>${esc(residue.residue_id)}</code> ${esc(residue.aa)}, DSSP <code>${esc(residue.dssp_index)}</code>, ACC <code>${esc(residue.accessibility ?? "not available")}</code>, BP <code>${esc([residue.bp1, residue.bp2].filter(Boolean).join(", ") || "none")}</code>.</span>`
      : "";
    details.innerHTML = `
      <strong>${esc(element.id)} ${esc(element.ss_name)}</strong>
      <span> chain <code>${esc(element.chain)}</code>, residues <code>${esc(element.start_residue)}-${esc(element.end_residue)}</code>,
      DSSP rows <code>${element.start_dssp}-${element.end_dssp}</code>, length <code>${element.length}</code>,
      mean ACC <code>${esc(acc)}</code>.</span>
      ${residueLine}
      <br><span>Sequence: <code>${esc(element.sequence)}</code></span>
    `;
    if (residue) highlightMolstarResidue(residue, "select");
  }

  function jumpToResidue() {
    const query = residueInput.value.trim().toLowerCase();
    if (!query) return;
    const elementId = elementByResidue.get(query);
    if (!elementId) {
      details.innerHTML = `No secondary-structure element contains residue <code>${esc(residueInput.value)}</code>.`;
      return;
    }
    selectElement(elementId, residueByLookup.get(query));
  }

  function applyTransform() {
    viewport.setAttribute("transform", `translate(${transform.x} ${transform.y}) scale(${transform.k})`);
  }

  function zoomBy(factor) {
    transform.k = Math.max(0.25, Math.min(5, transform.k * factor));
    applyTransform();
  }

  function resetView() {
    transform = { x: 0, y: 0, k: 1 };
    applyTransform();
  }

  function downloadJson() {
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `${(data.name || "dssp-topology").replace(/\W+/g, "_")}.json`;
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(url);
  }

  function loadCssOnce(url, id) {
    if (document.getElementById(id)) return;
    const link = document.createElement("link");
    link.id = id;
    link.rel = "stylesheet";
    link.href = url;
    document.head.appendChild(link);
  }

  function loadScriptOnce(url, id) {
    if (window.molstar) return Promise.resolve();
    const existing = document.getElementById(id);
    if (existing && existing.dataset.ready === "true") return Promise.resolve();
    if (existing && existing.dataset.loading === "true") {
      return new Promise((resolve, reject) => {
        existing.addEventListener("load", resolve, { once: true });
        existing.addEventListener("error", reject, { once: true });
      });
    }
    return new Promise((resolve, reject) => {
      const script = document.createElement("script");
      script.id = id;
      script.dataset.loading = "true";
      script.src = url;
      script.onload = () => {
        script.dataset.ready = "true";
        resolve();
      };
      script.onerror = reject;
      document.head.appendChild(script);
    });
  }

  function afdbModelId(accession) {
    if (!accession) return "";
    if (accession.toUpperCase().startsWith("AF-")) return accession;
    return `AF-${accession}-F1`;
  }

  async function fetchAfdbRecord(accession) {
    const response = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(accession)}`);
    if (!response.ok) throw new Error(`AlphaFold API returned ${response.status}`);
    const payload = await response.json();
    const record = Array.isArray(payload) ? payload[0] : payload;
    if (!record) return null;
    return {
      modelId: record.modelEntityId || record.entryId || afdbModelId(accession),
      cifUrl: record.cifUrl || record.cif_url || record.cifFileUrl || record.downloadUrl || null,
      pdbUrl: record.pdbUrl || record.pdb_url || null
    };
  }

  async function initMolstar() {
    if (!molstarNode || !data.afdb_accession) return;
    try {
      if (molstarStatus) molstarStatus.textContent = "Fetching AlphaFold metadata...";
      const record = await fetchAfdbRecord(data.afdb_accession).catch(() => null);
      loadCssOnce("https://cdn.jsdelivr.net/npm/molstar@5.4.2/build/viewer/molstar.css", "molstar-viewer-css");
      await loadScriptOnce("https://cdn.jsdelivr.net/npm/molstar@5.4.2/build/viewer/molstar.js", "molstar-viewer-js");
      if (molstarStatus) molstarStatus.textContent = "Loading 3D model...";
      molstarViewer = await window.molstar.Viewer.create(molstarNode, {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: true,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        viewportBackgroundColor: "white"
      });
      const primaryId = record?.modelId || afdbModelId(data.afdb_accession);
      if (molstarViewer.loadAlphaFoldDb) {
        await molstarViewer.loadAlphaFoldDb(primaryId).catch(async () => {
          await molstarViewer.loadAlphaFoldDb(data.afdb_accession);
        });
      } else if (record?.cifUrl) {
        await molstarViewer.loadStructureFromUrl(record.cifUrl, "mmcif", false);
      } else {
        await molstarViewer.loadStructureFromUrl(
          `https://alphafold.ebi.ac.uk/files/${primaryId}-model_v6.cif`,
          "mmcif",
          false
        );
      }
      if (molstarStatus) molstarStatus.textContent = "Hover or click topology residues to link 2D and 3D.";
    } catch (error) {
      if (molstarStatus) molstarStatus.textContent = `Mol* could not load: ${error.message || error}`;
    }
  }

  function highlightMolstarResidue(residue, action) {
    if (!molstarViewer || !molstarViewer.structureInteractivity || !residue) return;
    const seq = residue.resseq_int || Number.parseInt(residue.residue_number, 10);
    if (!Number.isFinite(seq)) return;
    const key = `${residue.chain}:${seq}:${action}`;
    if (action === "highlight" && key === lastMolstarResidue) return;
    lastMolstarResidue = key;
    try {
      molstarViewer.structureInteractivity({
        elements: {
          chain_id: residue.chain,
          auth_asym_id: residue.chain,
          label_asym_id: residue.chain,
          beg_auth_seq_id: seq,
          end_auth_seq_id: seq,
          beg_label_seq_id: seq,
          end_label_seq_id: seq
        },
        action
      });
    } catch (error) {
      if (molstarStatus) molstarStatus.textContent = `Mol* residue link failed: ${error.message || error}`;
    }
  }

  function clearMolstarHover() {
    lastMolstarResidue = "";
    try {
      molstarViewer?.plugin?.managers?.interactivity?.clearHighlights?.();
    } catch (error) {
      return;
    }
  }

  root.querySelector("[data-action='zoom-in']").addEventListener("click", () => zoomBy(1.2));
  root.querySelector("[data-action='zoom-out']").addEventListener("click", () => zoomBy(1 / 1.2));
  root.querySelector("[data-action='reset']").addEventListener("click", resetView);
  root.querySelector("[data-action='jump']").addEventListener("click", jumpToResidue);
  root.querySelector("[data-action='download-json']").addEventListener("click", downloadJson);
  residueInput.addEventListener("keydown", (event) => {
    if (event.key === "Enter") jumpToResidue();
  });
  linkToggle.addEventListener("change", () => {
    linksLayer.style.display = linkToggle.checked ? "" : "none";
  });
  labelToggle.addEventListener("change", () => {
    const display = labelToggle.checked ? "" : "none";
    root.querySelectorAll(".element-label,.residue-label,.terminus").forEach((node) => {
      node.style.display = display;
    });
  });

  svg.addEventListener("wheel", (event) => {
    event.preventDefault();
    zoomBy(event.deltaY < 0 ? 1.08 : 1 / 1.08);
  }, { passive: false });
  svg.addEventListener("mousedown", (event) => {
    if (event.button !== 0) return;
    drag = { x: event.clientX, y: event.clientY };
    svg.classList.add("is-dragging");
  });
  window.addEventListener("mousemove", (event) => {
    if (!drag) return;
    const dx = event.clientX - drag.x;
    const dy = event.clientY - drag.y;
    transform.x += dx;
    transform.y += dy;
    drag.x = event.clientX;
    drag.y = event.clientY;
    applyTransform();
  });
  window.addEventListener("mouseup", () => {
    drag = null;
    svg.classList.remove("is-dragging");
  });

  layoutElements();
  drawConnectors();
  drawLinks();
  drawElements();
  applyTransform();
  initMolstar();
})();
""".replace(
        "__ROOT_ID__", root_id
    ).replace(
        "__DATA_ID__", data_id
    )

    root_class = "has-molstar" if afdb_accession else ""
    return f"""
<div id="{root_id}" class="{root_class}">
  <style>{css}</style>
  <div class="topology-shell">
    <div class="topology-heading">
      <div class="title">{_html_escape(title)}</div>
      <div class="subtitle">{_html_escape(subtitle)}</div>
      <div class="stats">
        {stats["residue_count"]} residues;
        {stats["helix_count"]} helices;
        {stats["strand_count"]} strands;
        {stats["beta_link_count"]} beta links;
        chains: {_html_escape(", ".join(stats["chains"]))}
      </div>
    </div>
    <div class="viewer-grid">
      <section class="topology-card" aria-label="PDBe-style DSSP topology map">
        <div class="viewer-wrap">
          <svg role="img" aria-label="Interactive PDBe-style DSSP topology diagram">
            <g data-role="viewport">
              <g data-role="connectors"></g>
              <g data-role="links"></g>
              <g data-role="elements"></g>
            </g>
          </svg>
          <div class="tooltip" data-role="tooltip"></div>
        </div>
        <div class="pdbe-footer">
          <div class="entry-label">
            <i class="entry-dot" aria-hidden="true"></i>
            <span>{_html_escape(title)} | Entity 1 | Chain {_html_escape(", ".join(stats["chains"]))}</span>
          </div>
          <div class="controls" aria-label="Topology controls">
            <button type="button" data-action="zoom-in" title="Zoom in">+</button>
            <button type="button" data-action="zoom-out" title="Zoom out">-</button>
            <button type="button" data-action="reset" title="Reset pan and zoom">Reset</button>
            <select title="Annotation layer">
              <option>DSSP attributes</option>
              <option>Beta links</option>
            </select>
            <label class="toggle" title="Show DSSP beta bridge links">
              <input type="checkbox" data-role="toggle-links" checked> Links
            </label>
            <label class="toggle" title="Show element and residue labels">
              <input type="checkbox" data-role="toggle-labels" checked> Labels
            </label>
            <input type="text" data-role="residue-input" placeholder="Residue A:120">
            <button type="button" data-action="jump" title="Select a residue">Find</button>
            <button type="button" data-action="download-json" title="Download parsed topology JSON">JSON</button>
          </div>
        </div>
        <div class="details" data-role="details">
          Hover over a helix or strand to inspect residue-level DSSP attributes. Click a residue to pin the element and select it in Mol* when available.
        </div>
      </section>
      {molstar_panel}
    </div>
  </div>
  <script type="application/json" id="{data_id}">{json_blob}</script>
  <script>{script}</script>
</div>
"""


def _html_escape(value: Any) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _decode_upload(value: Any) -> Optional[Tuple[str, str]]:
    if not value:
        return None

    if isinstance(value, dict):
        name, item = next(iter(value.items()))
    elif isinstance(value, (tuple, list)):
        item = value[0]
        name = item.get("name", "uploaded.dssp")
    else:
        return None

    content = item.get("content")
    if isinstance(content, memoryview):
        content = content.tobytes()
    if isinstance(content, str):
        text = content
    else:
        text = bytes(content).decode("utf-8", errors="replace")
    return str(name), text


def make_app(default_path: str = DEFAULT_DSSP_PATH) -> Any:
    import ipywidgets as widgets
    from IPython.display import HTML, clear_output, display

    state: Dict[str, Any] = {"name": None, "text": None, "kind": None, "metadata": None}

    uploader = widgets.FileUpload(
        accept=".cif,.mmcif,.dssp",
        multiple=False,
        description="Upload file",
        layout=widgets.Layout(width="180px"),
    )
    afdb_input = widgets.Text(
        value="",
        placeholder="P07949",
        description="AFDB AC",
        tooltip="UniProt/AlphaFold DB accession used to fetch AlphaFold v6 mmCIF and Mol*",
        layout=widgets.Layout(width="220px"),
    )
    fetch_afdb_button = widgets.Button(
        description="Fetch AFDB topology",
        button_style="success",
        tooltip="Fetch the AlphaFold DB mmCIF and derive a PDBe-style topology",
        layout=widgets.Layout(width="180px"),
    )
    visualize_button = widgets.Button(
        description="Visualize topology",
        button_style="primary",
        tooltip="Parse the loaded CIF/DSSP file and draw the topology",
        layout=widgets.Layout(width="180px"),
    )
    example_button = widgets.Button(
        description="Load P07949 example",
        tooltip="Load the DSSP file path supplied with this project",
        disabled=not Path(default_path).exists(),
        layout=widgets.Layout(width="180px"),
    )
    status = widgets.HTML(
        value="<span style='color:#5d6978'>Enter an AFDB accession and fetch topology, or upload a .cif/.mmcif/.dssp file.</span>"
    )
    output = widgets.Output(
        layout=widgets.Layout(border="1px solid #d4dce6", min_height="680px")
    )

    def _kind_from_name(name: str) -> str:
        suffix = Path(name).suffix.lower()
        if suffix in {".cif", ".mmcif"}:
            return "cif"
        return "dssp"

    def set_loaded(
        name: str,
        text: str,
        kind: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        state["name"] = name
        state["text"] = text
        state["kind"] = kind or _kind_from_name(name)
        state["metadata"] = metadata
        if not afdb_input.value.strip():
            afdb_input.value = _guess_afdb_accession(name)
        status.value = (
            f"<span style='color:#1f6f43'>Loaded <b>{_html_escape(name)}</b> "
            f"({len(text.splitlines())} lines, {state['kind']}).</span>"
        )

    def on_upload(change: Dict[str, Any]) -> None:
        decoded = _decode_upload(change.get("new"))
        if decoded is None:
            return
        set_loaded(*decoded)

    def visualize(_: Any = None) -> None:
        if not state["text"]:
            status.value = (
                "<span style='color:#a15c00'>No CIF or DSSP file is loaded yet.</span>"
            )
            return
        afdb_accession = afdb_input.value.strip() or _guess_afdb_accession(
            state["name"] or ""
        )
        try:
            if state.get("kind") == "cif":
                topology = topology_from_alphafold_cif(
                    state["text"],
                    state["name"] or "uploaded.cif",
                    afdb_metadata=state.get("metadata"),
                )
            else:
                topology = topology_from_dssp(
                    state["text"],
                    state["name"] or "uploaded.dssp",
                    afdb_accession=afdb_accession,
                )
            html = topology_html(topology)
        except Exception as exc:  # noqa: BLE001 - surface parser problems in the app.
            status.value = f"<span style='color:#b42318'>Could not parse topology source: {_html_escape(exc)}</span>"
            return

        source_label = "PDBe renderer topology JSON" if topology.get("pdbe_api_data") else "DSSP fallback"
        status.value = (
            f"<span style='color:#1f6f43'>Built {source_label}: {topology['stats']['residue_count']} residues, "
            f"{topology['stats']['helix_count']} helices, "
            f"{topology['stats']['strand_count']} strands"
            f"{'; Mol* linked to ' + _html_escape(afdb_accession) if afdb_accession else ''}.</span>"
        )
        with output:
            clear_output(wait=True)
            display(HTML(html))

    def load_example(_: Any = None) -> None:
        path = Path(default_path)
        if not path.exists():
            status.value = f"<span style='color:#b42318'>Example file not found: {_html_escape(path)}</span>"
            return
        afdb_input.value = _guess_afdb_accession(path.name)
        set_loaded(path.name, path.read_text(encoding="utf-8", errors="replace"))
        visualize()

    def fetch_afdb(_: Any = None) -> None:
        accession = afdb_input.value.strip()
        if not accession:
            status.value = "<span style='color:#a15c00'>Enter an AlphaFold DB accession first.</span>"
            return
        status.value = f"<span style='color:#5d6978'>Fetching AlphaFold DB mmCIF for {_html_escape(accession)}...</span>"
        try:
            metadata, cif_text = fetch_alphafold_cif(accession)
            name = f"{metadata.get('modelEntityId') or metadata.get('entryId') or accession}.cif"
            set_loaded(name, cif_text, kind="cif", metadata=metadata)
            visualize()
        except Exception as exc:  # noqa: BLE001 - show network/parser failures in widget.
            status.value = f"<span style='color:#b42318'>Could not fetch AlphaFold DB model: {_html_escape(exc)}</span>"

    uploader.observe(on_upload, names="value")
    visualize_button.on_click(visualize)
    fetch_afdb_button.on_click(fetch_afdb)
    example_button.on_click(load_example)

    app_style = widgets.HTML(
        """
        <style>
        .dssp-app-title {
          color: #1f2933;
          font-family: Inter, "Segoe UI", Arial, sans-serif;
          margin: 0 0 10px 0;
        }
        .dssp-app-title h2 {
          font-size: 22px;
          margin: 0 0 4px 0;
        }
        .dssp-app-title p {
          color: #5d6978;
          font-size: 13px;
          margin: 0;
        }
        </style>
        """
    )
    header = widgets.HTML(
        """
        <div class="dssp-app-title">
          <h2>DSSP Topology Viewer</h2>
          <p>Fetch an AlphaFold mmCIF, derive secondary-structure topology, and render it with the PDBe topology viewer architecture.</p>
        </div>
        """
    )
    controls = widgets.HBox(
        [afdb_input, fetch_afdb_button, uploader, visualize_button, example_button],
        layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
    )
    return widgets.VBox(
        [app_style, header, controls, status, output],
        layout=widgets.Layout(gap="10px", width="100%"),
    )


if __name__ == "__main__":
    from IPython.display import display

    app = make_app()
    display(app)
