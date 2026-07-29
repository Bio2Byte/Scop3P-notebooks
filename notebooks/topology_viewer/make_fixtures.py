"""Build synthetic structures with known secondary structure.

Ideal helix and strand geometry, so the SS assignment can be scored against a
ground truth we control. Emits PDB and mmCIF, each with and without secondary
structure records, covering the four input shapes the loader has to handle.
"""

import math
from pathlib import Path

HELIX_RADIUS = 2.3
HELIX_RISE = 1.5
HELIX_TWIST = math.radians(100.0)
STRAND_RISE = 3.3
STRAND_ZIGZAG = 0.9


def helix(count, origin, axis_index=2):
    points = []
    for i in range(count):
        angle = i * HELIX_TWIST
        local = [HELIX_RADIUS * math.cos(angle), HELIX_RADIUS * math.sin(angle), i * HELIX_RISE]
        local[axis_index], local[2] = local[2], local[axis_index]
        points.append(tuple(origin[k] + local[k] for k in range(3)))
    return points


def strand(count, origin, direction=1):
    points = []
    for i in range(count):
        offset = STRAND_ZIGZAG if i % 2 else -STRAND_ZIGZAG
        points.append((origin[0] + offset, origin[1], origin[2] + direction * i * STRAND_RISE))
    return points


def coil(start, end, count):
    points = []
    for i in range(1, count + 1):
        t = i / (count + 1.0)
        points.append(tuple(start[k] + (end[k] - start[k]) * t + (2.0 if i % 2 else -2.0) for k in range(3)))
    return points


def build():
    """Beta sheet of four antiparallel strands, with two helices packed above."""
    residues = []
    truth = []

    def add(points, code):
        for point in points:
            residues.append(point)
            truth.append(code)

    h1 = helix(14, (-14.0, 12.0, 0.0), axis_index=0)
    add(h1, "H")
    add(coil(h1[-1], (0.0, 0.0, 0.0), 4), "C")

    s1 = strand(8, (0.0, 0.0, 0.0), 1)
    add(s1, "E")
    add(coil(s1[-1], (4.8, 0.0, 23.1), 5), "C")

    s2 = strand(8, (4.8, 0.0, 23.1), -1)
    add(s2, "E")
    add(coil(s2[-1], (9.6, 0.0, 0.0), 5), "C")

    h2 = helix(12, (2.0, 14.0, -6.0), axis_index=0)
    add(h2, "H")
    add(coil(h2[-1], (9.6, 0.0, 0.0), 4), "C")

    s3 = strand(8, (9.6, 0.0, 0.0), 1)
    add(s3, "E")
    add(coil(s3[-1], (14.4, 0.0, 23.1), 5), "C")

    s4 = strand(8, (14.4, 0.0, 23.1), -1)
    add(s4, "E")

    return residues, truth


def to_ranges(truth):
    ranges = []
    current, start = None, 0
    for index, code in enumerate(truth + ["C"]):
        if code != current:
            if current in {"H", "E"}:
                ranges.append((current, start + 1, index))
            current, start = code, index
    return ranges


def write_pdb(path, residues, truth, with_ss):
    lines = ["HEADER    SYNTHETIC TEST                          01-JAN-26   TEST"]
    if with_ss:
        helix_n = sheet_n = 0
        for kind, start, stop in to_ranges(truth):
            if kind == "H":
                helix_n += 1
                lines.append(
                    f"HELIX  {helix_n:3d} {helix_n:3d} ALA A {start:4d}  ALA A {stop:4d}  1"
                    f"{'':30}{stop - start + 1:5d}"
                )
            else:
                sheet_n += 1
                # Columns per the PDB spec: sheetID 12-14, numStrands 15-16,
                # initChainID 22, initSeqNum 23-26, endChainID 33, endSeqNum 34-37.
                lines.append(
                    f"SHEET  {sheet_n:3d} {'A':>3s}{4:2d} ALA A{start:4d}  ALA A{stop:4d}  0"
                )
    for index, (x, y, z) in enumerate(residues, start=1):
        lines.append(
            f"ATOM  {index:5d}  CA  ALA A{index:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 85.00           C"
        )
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")


def write_cif(path, residues, truth, with_ss):
    lines = ["data_TEST", "#", "_entry.id TEST", "#"]
    if with_ss:
        conf = [entry for entry in to_ranges(truth) if entry[0] == "H"]
        if conf:
            lines += [
                "loop_",
                "_struct_conf.conf_type_id",
                "_struct_conf.id",
                "_struct_conf.beg_label_asym_id",
                "_struct_conf.beg_label_seq_id",
                "_struct_conf.beg_auth_asym_id",
                "_struct_conf.beg_auth_seq_id",
                "_struct_conf.end_label_asym_id",
                "_struct_conf.end_label_seq_id",
                "_struct_conf.end_auth_asym_id",
                "_struct_conf.end_auth_seq_id",
            ]
            for i, (_, start, stop) in enumerate(conf, start=1):
                lines.append(
                    f"HELX_P HELX_P{i} A {start} A {start} A {stop} A {stop}"
                )
            lines.append("#")
        sheets = [entry for entry in to_ranges(truth) if entry[0] == "E"]
        if sheets:
            lines += [
                "loop_",
                "_struct_sheet_range.sheet_id",
                "_struct_sheet_range.id",
                "_struct_sheet_range.beg_auth_asym_id",
                "_struct_sheet_range.beg_auth_seq_id",
                "_struct_sheet_range.end_auth_asym_id",
                "_struct_sheet_range.end_auth_seq_id",
            ]
            for i, (_, start, stop) in enumerate(sheets, start=1):
                lines.append(f"A {i} A {start} A {stop}")
            lines.append("#")

    lines += [
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_seq_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_PDB_model_num",
    ]
    for index, (x, y, z) in enumerate(residues, start=1):
        lines.append(
            f"ATOM {index} CA ALA A {index} A {index} "
            f"{x:.3f} {y:.3f} {z:.3f} 85.00 1"
        )
    lines.append("#")
    Path(path).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    residues, truth = build()
    Path("fixtures").mkdir(exist_ok=True)
    write_pdb("fixtures/annotated.pdb", residues, truth, True)
    write_pdb("fixtures/bare.pdb", residues, truth, False)
    write_cif("fixtures/annotated.cif", residues, truth, True)
    write_cif("fixtures/bare.cif", residues, truth, False)
    Path("fixtures/truth.txt").write_text("".join(truth) + "\n")
    print(f"{len(residues)} residues")
    print("".join(truth))
    print(to_ranges(truth))
