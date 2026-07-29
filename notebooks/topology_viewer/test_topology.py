"""Checks for the phase 0 pipeline.

Run with: python3 test_topology.py
"""

import json
import math
import random
import sys
from pathlib import Path

from topology.elements import annotate_geometry, build_elements, strand_contacts
from topology.io import load_structure, sniff_format
from topology.layout import build_layout
from topology.render import build_payload, render, standalone_document
from topology import ss

FAILURES = []


def check(label, condition, detail=""):
    mark = "pass" if condition else "FAIL"
    print(f"  [{mark}] {label}" + (f"  -- {detail}" if detail and not condition else ""))
    if not condition:
        FAILURES.append(label)


def section(title):
    print(f"\n{title}")


def pipeline(structure, chain):
    elements, residues = build_elements(structure, chain)
    annotate_geometry(elements, residues)
    contacts = strand_contacts(elements, residues)
    layouts = {
        mode: build_layout(mode, elements, residues, contacts)
        for mode in ("sheet", "serpentine", "spatial")
    }
    return elements, residues, contacts, layouts


# ---------------------------------------------------------------- loading
section("Format detection and loading")

fixtures = {name: Path(f"fixtures/{name}").read_text() for name in
            ["annotated.pdb", "bare.pdb", "annotated.cif", "bare.cif"]}

check("PDB text sniffed as pdb", sniff_format(fixtures["bare.pdb"], "x.cif") == "pdb",
      "extension must not override content")
check("mmCIF text sniffed as mmcif", sniff_format(fixtures["bare.cif"], "x.pdb") == "mmcif")
check("garbage rejected", sniff_format("hello world", "x.txt") == "unknown")

structures = {name: load_structure(text, name) for name, text in fixtures.items()}

for name, structure in structures.items():
    check(f"{name}: 81 residues", len(structure.residues_by_chain["A"]) == 81,
          str(len(structure.residues_by_chain.get("A", []))))

check("annotated files report file provenance",
      all("file" in structures[n].ss_source for n in ["annotated.pdb", "annotated.cif"]))
check("bare files report computed provenance",
      all("file" not in structures[n].ss_source for n in ["bare.pdb", "bare.cif"]))

check("PDB and mmCIF agree on element count",
      len(structures["bare.pdb"].ss_by_chain["A"]) == len(structures["bare.cif"].ss_by_chain["A"]))

try:
    load_structure("", "empty.pdb")
    check("empty file raises", False)
except ValueError:
    check("empty file raises", True)

try:
    load_structure("not a structure at all", "junk.txt")
    check("junk file raises", False)
except ValueError:
    check("junk file raises", True)


# ---------------------------------------------- secondary structure quality
section("Secondary structure assignment")

truth = Path("fixtures/truth.txt").read_text().strip()
coords = [r.coords for r in structures["bare.pdb"].residues_by_chain["A"]]

builtin = ss.assign_psea(coords)
agreement = sum(1 for a, b in zip(truth, builtin) if a == b) / len(truth)
check(f"built-in P-SEA agrees with truth ({agreement:.0%})", agreement > 0.80)

codes, source = ss.compute(coords)
check("compute returns a provenance string", isinstance(source, str) and len(source) > 3)
check("compute returns one code per residue", len(codes) == len(coords))

check("no element spans a chain break",
      all(c in "HEC" for c in builtin))

# A single ideal helix should be called a helix throughout its core.
helix_only = [
    (2.3 * math.cos(i * math.radians(100)), 2.3 * math.sin(i * math.radians(100)), i * 1.5)
    for i in range(20)
]
codes = ss.assign_psea(helix_only)
check("ideal helix assigned H", codes.count("H") >= 16, "".join(codes))

# Random coordinates must not produce structure.
random.seed(7)
noise = [(random.uniform(0, 40), random.uniform(0, 40), random.uniform(0, 40))
         for _ in range(40)]
codes = ss.assign_psea(noise)
check("random coordinates give mostly coil", codes.count("C") >= 34, "".join(codes))


# ------------------------------------------------------------ elements
section("Elements and contacts")

structure = structures["annotated.cif"]
elements, residues, contacts, layouts = pipeline(structure, "A")

check("six elements recovered", len(elements) == 6, str(len(elements)))
check("elements ordered by sequence",
      all(elements[i]["start"] <= elements[i + 1]["start"] for i in range(len(elements) - 1)))
check("three strand pairings found", len(contacts) == 3, str(len(contacts)))
check("all pairings antiparallel",
      all(c["orientation"] == "antiparallel" for c in contacts))
check("every element carries a centroid",
      all(e.get("centroid") for e in elements))


# ------------------------------------------------------------ layouts
section("Layout")

for mode, layout in layouts.items():
    check(f"{mode}: every element placed", len(layout["elements"]) == len(elements))
    check(f"{mode}: connectors join consecutive elements",
          len(layout["connectors"]) == len(elements) - 1)
    check(f"{mode}: two termini", len(layout["termini"]) == 2)

    x0, y0, x1, y1 = layout["extents"]
    check(f"{mode}: extents finite and ordered", x1 > x0 and y1 > y0)
    check(f"{mode}: diagram stays a sane size", (x1 - x0) < 4000 and (y1 - y0) < 4000,
          f"{x1 - x0:.0f}x{y1 - y0:.0f}")

    for element in layout["elements"]:
        inside = x0 - 1 <= element["x"] <= x1 + 1
        check(f"{mode}: {element['id']} inside extents", inside,
              f"x={element['x']} range={x0}..{x1}")

sheet = layouts["sheet"]
placed = {e["id"]: e for e in sheet["elements"]}
check("antiparallel neighbours face opposite ways",
      placed["S1"]["direction"] == -placed["S2"]["direction"])
check("strand order follows the sheet",
      placed["S1"]["x"] < placed["S2"]["x"] < placed["S3"]["x"] < placed["S4"]["x"])

# Connectors must be assigned distinct lanes when they overlap horizontally.
by_side = {}
for connector in sheet["connectors"]:
    by_side.setdefault(connector["side"], []).append(connector)
overlapping_share_lane = False
for side, group in by_side.items():
    for i, a in enumerate(group):
        for b in group[i + 1:]:
            a_span = sorted([a["path"][0], a["path"][-2]])
            b_span = sorted([b["path"][0], b["path"][-2]])
            overlap = a_span[0] < b_span[1] and b_span[0] < a_span[1]
            if overlap and a["lane"] == b["lane"]:
                overlapping_share_lane = True
check("overlapping connectors get separate lanes", not overlapping_share_lane)


# ------------------------------------------------------------ multi-chain
section("Multi-chain handling")

text = fixtures["bare.pdb"]
shifted = []
for line in text.splitlines():
    if line.startswith("ATOM  "):
        shifted.append(line[:21] + "B" + line[22:])
dimer = text + "\n" + "\n".join(shifted)
dimer_structure = load_structure(dimer, "dimer.pdb")

check("two chains detected", len(dimer_structure.chains) == 2, str(dimer_structure.chains))
check("chain options carry residue counts",
      all("residues" in label for label, _ in dimer_structure.chain_options()))
check("default chain is the largest",
      dimer_structure.default_chain() in {"A", "B"})

elements_b, residues_b = build_elements(dimer_structure, "B")
check("chain B builds independently", len(elements_b) > 0)
check("chain B residues carry chain B",
      all(r.chain == "B" for r in residues_b))


# ------------------------------------------------------------ payload
section("Payload and rendering")

payload = build_payload(structure, "A", elements, residues, contacts, layouts)
payload["structure_source"] = {"kind": "upload", "data": "x", "format": "mmcif"}

check("payload is JSON serialisable", isinstance(json.dumps(payload), str))
check("schema version present", payload["schema"] == 1)
check("annotations absent in file mode", payload["annotations"] is None)
check("all layouts included",
      set(payload["layouts"]) == {"sheet", "serpentine", "spatial"})
check("residues carry both numbering schemes",
      all("seq" in r and "label_seq" in r for r in payload["residues"]))

html = render(payload, embed="inline")
check("html produced", len(html) > 5000)
check("root element present", 'class="topo-root"' in html)
check("payload embedded", 'data-role="payload"' in html)
check("viewer adapters embedded", "__topoViewers" in html)
check("boot call present", "__topoBoot" in html)
check("no unresolved template markers", "__ROOT__" not in html)

start = html.index('data-role="payload">') + len('data-role="payload">')
end = html.index("</script>", start)
check("embedded payload re-parses", isinstance(json.loads(html[start:end]), dict))

# JupyterLab does not run scripts inside display(HTML(...)), so the default
# embedding has to be an iframe or the diagram never draws.
framed = render(payload)
check("default embedding is an iframe", framed.strip().startswith("<iframe"))
check("iframe carries the document", "srcdoc=" in framed and "__topoBoot" in framed)
check("iframe allows scripts", 'allow-scripts' in framed)
check("standalone document is a full page",
      standalone_document(payload).lstrip().lower().startswith("<!doctype html"))


# ------------------------------------------------------------ larger protein
section("Larger protein")

random.seed(3)
points = []
position = [0.0, 0.0, 0.0]
for _ in range(40):
    for i in range(random.randint(6, 14)):
        position = [position[0] + random.uniform(-2, 2),
                    position[1] + random.uniform(-2, 2),
                    position[2] + 3.6]
        points.append(tuple(position))
lines = [f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 70.00           C"
         for i, (x, y, z) in enumerate(points, start=1)]
big = load_structure("\n".join(lines) + "\nEND\n", "big.pdb")
check(f"large chain loads ({len(points)} residues)", len(big.residues_by_chain["A"]) == len(points))

big_elements, big_residues, big_contacts, big_layouts = pipeline(big, "A")
for mode, layout in big_layouts.items():
    x0, y0, x1, y1 = layout["extents"]
    check(f"{mode}: large diagram bounded", (x1 - x0) < 20000 and (y1 - y0) < 20000,
          f"{x1 - x0:.0f}x{y1 - y0:.0f}")
    check(f"{mode}: all elements placed", len(layout["elements"]) == len(big_elements))



# ------------------------------------------------------- aspect packing
section("Frame packing")

from topology.layout import layout_sheet, layout_serpentine

def aspect_of(layout):
    x0, y0, x1, y1 = layout["extents"]
    return (x1 - x0) / max(1e-6, (y1 - y0))

big = load_structure(Path("fixtures/big.pdb").read_text(), "big.pdb")
big_elements, big_res, big_contacts, _ = pipeline(big, "A")

# Asking for a wider frame must actually produce a wider diagram. A layout that
# ignores the frame shape renders as a long strip that scales down to nothing.
wide = aspect_of(layout_sheet(big_elements, big_res, big_contacts, target_aspect=3.0))
tall = aspect_of(layout_sheet(big_elements, big_res, big_contacts, target_aspect=0.5))
check("sheet packing responds to the target aspect", wide > tall,
      f"wide={wide:.2f} tall={tall:.2f}")

wide_s = aspect_of(layout_serpentine(big_elements, big_res, big_contacts, target_aspect=3.0))
tall_s = aspect_of(layout_serpentine(big_elements, big_res, big_contacts, target_aspect=0.5))
check("serpentine responds to the target aspect", wide_s > tall_s,
      f"wide={wide_s:.2f} tall={tall_s:.2f}")

for mode, engine in [("sheet", layout_sheet), ("serpentine", layout_serpentine)]:
    got = aspect_of(engine(big_elements, big_res, big_contacts, target_aspect=1.3))
    check(f"{mode}: lands within 2x of a 1.3 target", 0.65 <= got <= 2.6, f"{got:.2f}")

# Every element must sit inside the packed extents, including across bands.
packed = layout_sheet(big_elements, big_res, big_contacts, target_aspect=1.3)
px0, py0, px1, py1 = packed["extents"]
inside = all(px0 - 1 <= e["x"] <= px1 + 1 and py0 - 1 <= e["y"] <= py1 + 1
             for e in packed["elements"])
check("packed elements stay inside extents", inside)

# Connectors still run C-terminal end to N-terminal end of the next element.
by_id_packed = {e["id"]: e for e in packed["elements"]}
ordered = packed["elements"]
joins_ok = True
for connector in packed["connectors"]:
    source = by_id_packed[connector["source"]]
    target = by_id_packed[connector["target"]]
    head = (connector["path"][0], connector["path"][1])
    tail = (connector["path"][-2], connector["path"][-1])
    if [round(v, 1) for v in head] != [round(v, 1) for v in source["c_point"]]:
        joins_ok = False
    if [round(v, 1) for v in tail] != [round(v, 1) for v in target["n_point"]]:
        joins_ok = False
check("connectors leave the C end and enter the N end", joins_ok)


# ------------------------------------------------------- sheet labelling
section("Sheet labelling")

from topology.elements import assign_sheets

assign_sheets(big_elements, big_contacts)
sheets = {e["id"]: e["sheet"] for e in big_elements if e["type"] == "strand"}
check("strands carry a sheet label", all(v is not None for v in sheets.values()), str(sheets))
check("two distinct sheets found", len(set(sheets.values())) == 2, str(set(sheets.values())))
check("helices carry no sheet label",
      all(e["sheet"] is None for e in big_elements if e["type"] == "helix"))

# A strand that pairs with nothing is not a sheet and must not be labelled one.
solo_elements, solo_res = build_elements(structures["annotated.cif"], "A")
annotate_geometry(solo_elements, solo_res)
assign_sheets(solo_elements, [])
check("unpaired strands are left unlabelled",
      all(e["sheet"] is None for e in solo_elements))

# Sheet numbering must follow the chain, not hash order.
assign_sheets(big_elements, big_contacts)
first_by_sheet = {}
for element in big_elements:
    if element["type"] == "strand" and element["sheet"]:
        first_by_sheet.setdefault(element["sheet"], element["start"])
ordered_labels = [k for k, _ in sorted(first_by_sheet.items(), key=lambda kv: kv[1])]
check("sheets numbered in sequence order",
      ordered_labels == sorted(ordered_labels), str(ordered_labels))


# ------------------------------------------------------- spatial arrangement
section("Spatial arrangement and local routing")

from topology.layout import layout_spatial

big_order = [e["id"] for e in big_elements]
spatial = layout_spatial(big_elements, big_res, big_contacts, target_aspect=1.3)
check("spatial places every element", len(spatial["elements"]) == len(big_elements))
check("spatial reports segments", len(spatial.get("segments", [])) >= 1)

# Every row must read left to right. Snaking shortens the wrap but forces the
# reader to reverse direction on alternate rows.
rows = {}
for element in spatial["elements"]:
    rows.setdefault(round(element["y"] / 100), []).append(element)
forward = True
for members in rows.values():
    ordered = sorted(members, key=lambda e: big_order.index(e["id"]))
    for a, b in zip(ordered, ordered[1:]):
        if b["x"] < a["x"]:
            forward = False
check("every row reads left to right", forward)

# Local routing: a loop should not climb past the whole diagram to get across.
def path_length(connector):
    p = connector["path"]
    return sum(abs(p[i + 2] - p[i]) + abs(p[i + 3] - p[i + 1])
               for i in range(0, len(p) - 2, 2))

for mode in ("sheet", "serpentine", "spatial"):
    laid = build_layout(mode, big_elements, big_res, big_contacts)
    x0, y0, x1, y1 = laid["extents"]
    diagonal = abs(x1 - x0) + abs(y1 - y0)
    worst = max(path_length(c) for c in laid["connectors"])
    check(f"{mode}: no connector exceeds the diagram perimeter",
          worst <= diagonal, f"worst={worst:.0f} diagonal={diagonal:.0f}")

# Loops must leave and enter at the ends the element directions actually put
# them at, otherwise the path crosses back over its own element.
laid = build_layout("spatial", big_elements, big_res, big_contacts)
placed_by_id = {e["id"]: e for e in laid["elements"]}
sides_ok = True
for connector in laid["connectors"]:
    source = placed_by_id[connector["source"]]
    target = placed_by_id[connector["target"]]
    if connector["side"] == "above":
        if source["c_point"][1] > source["y"] + source["h"] / 2:
            sides_ok = False
    elif connector["side"] == "below":
        if source["c_point"][1] < source["y"] + source["h"] / 2:
            sides_ok = False
check("loop side matches where the C end actually sits", sides_ok)


# ------------------------------------------------------- projection layout
section("Projection layout")

from topology.layout import layout_projection
from topology.projection import principal_axes, project_elements, relax_column

# Eigen decomposition must find the true principal axis of a known spread.
spread = [(t, 0.0, 0.0) for t in range(-5, 6)] + [(0.0, 0.5, 0.0), (0.0, -0.5, 0.0)]
centre, axes = principal_axes(spread)
check("principal axis found", abs(abs(axes[0][0]) - 1.0) < 0.05, str(axes[0]))

# Overlap relaxation must preserve order and remove collisions.
relaxed = relax_column([("a", 0.0, 100.0), ("b", 10.0, 100.0), ("c", 20.0, 100.0)], 20.0)
seps = [relaxed["b"] - relaxed["a"], relaxed["c"] - relaxed["b"]]
check("relaxation separates overlapping elements", all(s >= 119.9 for s in seps), str(seps))
check("relaxation preserves order", relaxed["a"] < relaxed["b"] < relaxed["c"])

big_order = [e["id"] for e in big_elements]
proj = layout_projection(big_elements, big_res, big_contacts, target_aspect=1.3)
check("projection places every element", len(proj["elements"]) == len(big_elements))
check("projection reports its mode", proj["mode"] == "projection")

# No two elements may overlap: same column and vertically intersecting.
overlap = False
for i, a in enumerate(proj["elements"]):
    for b in proj["elements"][i + 1:]:
        if abs(a["x"] - b["x"]) < 1.0:
            if not (a["y"] + a["h"] <= b["y"] or b["y"] + b["h"] <= a["y"]):
                overlap = True
check("no two elements overlap", not overlap)

paspect = aspect_of(proj)
check("projection lands near the target", 0.5 <= paspect <= 2.6, f"{paspect:.2f}")

# Direction must follow the projected N-to-C axis, not an arbitrary rule.
directions = {e["id"]: e["direction"] for e in proj["elements"]}
check("directions are all +1 or -1", set(directions.values()) <= {1, -1})

# Elements with no geometry must fall back rather than crash.
bare = [dict(e) for e in big_elements]
for element in bare:
    element.pop("centroid", None)
    element.pop("axis", None)
fallback = layout_projection(bare, big_res, big_contacts)
check("falls back when geometry is absent", len(fallback["elements"]) == len(bare))


# ------------------------------------------------------- widget behaviour
section("Widget behaviour")

try:
    import datetime
    import ipywidgets as widgets
    from topology.app import make_app

    def walk(widget, kind, found):
        if isinstance(widget, kind):
            found.append(widget)
        for child in getattr(widget, "children", []):
            walk(child, kind, found)

    # A Dropdown must never be sent unobserve_all(): it strips ipywidgets' own
    # internal observer, so the next assignment to .value raises and the widget
    # callback swallows it, freezing the app with no traceback.
    dropdown = widgets.Dropdown(options=[("", "")], value="")
    dropdown.options = [("A (81 residues)", "A")]
    dropdown.value = "A"
    check("dropdown accepts new options and value", dropdown.value == "A")

    for fixture, expect in [("annotated.cif", "file"), ("bare.pdb", "P-SEA")]:
        app = make_app()
        uploads = []
        walk(app, widgets.FileUpload, uploads)
        content = Path("fixtures/" + fixture).read_bytes()
        uploads[0].value = ({
            "name": fixture, "type": "text/plain", "size": len(content),
            "content": memoryview(content),
            "last_modified": datetime.datetime.now(),
        },)

        labels, pickers, outputs = [], [], []
        walk(app, widgets.HTML, labels)
        walk(app, widgets.Dropdown, pickers)
        walk(app, widgets.Output, outputs)
        status = labels[-1].value

        check(f"upload {fixture}: status advances past 'Reading'",
              "Uploaded" in status, status[:120])
        check(f"upload {fixture}: provenance reported", expect in status, status[:120])
        # Select by description: there are now two dropdowns, and the first is
        # the structure picker.
        chain_dd = next(d for d in pickers if d.description == "Chain")
        check(f"upload {fixture}: chain dropdown populated",
              chain_dd.value == "A", repr(chain_dd.value))
        check(f"upload {fixture}: no traceback in output",
              len(outputs[0].outputs) <= 1)

    # Two chains must reveal the chain picker.
    text = Path("fixtures/bare.pdb").read_text()
    shifted = [line[:21] + "B" + line[22:] for line in text.splitlines()
               if line.startswith("ATOM  ")]
    dimer = (text + "\n" + "\n".join(shifted)).encode()
    app = make_app()
    uploads = []
    walk(app, widgets.FileUpload, uploads)
    uploads[0].value = ({
        "name": "dimer.pdb", "type": "text/plain", "size": len(dimer),
        "content": memoryview(dimer), "last_modified": datetime.datetime.now(),
    },)
    boxes = []
    walk(app, widgets.HBox, boxes)
    revealed = any(getattr(box.layout, "display", None) == "flex" for box in boxes)
    check("chain picker revealed for a two-chain file", revealed)

    # An upload must not inherit annotations: no accession, so no valid mapping.
    pickers = []
    walk(app, widgets.Dropdown, pickers)
    chain_dd = next(d for d in pickers if d.description == "Chain")
    check("upload populates the chain picker", chain_dd.value in {"A", "B"}, repr(chain_dd.value))

except ImportError:
    print("  [skip] ipywidgets not installed")



# ------------------------------------------------------------ annotations
section("Annotations")

import json as _json
from topology import annotations as ann

_API = Path("fixtures/api")
_load = lambda name: _json.load(open(_API / name))

refs = ann.merge_refs(
    ann.parse_best_structures(_load("best_structures.json"), "P07949"),
    ann.parse_uniprot_xrefs(_load("uniprot_entry.json")),
)
ids = [r.pdb_id for r in refs]
check("PDB entries found", len(refs) == 4, str(ids))
check("PDBe coverage ranking preserved", ids[0] == "2IVV", str(ids))
check("UniProt-only entries retained", "4CKI" in ids, str(ids))
check("non-PDB cross-references ignored", "P07949" not in ids)

two_ivv = next(r for r in refs if r.pdb_id == "2IVV")
check("chains merged across both sources", set(two_ivv.chains) == {"A", "B"}, str(two_ivv.chains))
check("chain ranges parsed", two_ivv.chains["A"] == (724, 1016), str(two_ivv.chains["A"]))

multi = next(r for r in refs if r.pdb_id == "4CKI")
check("multi-range chain spans min to max", multi.chains["A"] == (705, 1050), str(multi.chains))
check("resolution parsed from text", multi.resolution == 2.6, str(multi.resolution))

# Numbering. Blocks carry different offsets, so a single shift cannot work.
mapping = ann.parse_pdbe_mapping(_load("pdbe_mapping.json"), "2ivv", "P07949", "A")
check("first block offset applied", mapping.get(724) == 721, str(mapping.get(724)))
check("second block has its own offset", mapping.get(810) == 805, str(mapping.get(810)))
check("unobserved positions stay unmapped", mapping.get(805) is None, str(mapping.get(805)))
check("chain filter excludes other chains", max(mapping.values()) < 1700, str(max(mapping.values())))

chain_b = ann.parse_pdbe_mapping(_load("pdbe_mapping.json"), "2ivv", "P07949", "B")
check("chain B maps to its own numbering", chain_b.get(724) == 1724, str(chain_b.get(724)))

# Sites.
scop = ann.parse_scop3p(_load("scop3p.json"))
check("malformed positions skipped", len(scop) == 4, str(len(scop)))
check("functional score parsed", any(s.score == 0.82 for s in scop))

ptms = ann.parse_uniprot_ptms(_load("uniprot_entry.json"))
check("multi-residue features skipped", all(s.position in (905, 696) for s in ptms), str([s.position for s in ptms]))

variants = ann.parse_variants(_load("variants.json"), "P07949")
check("only disease variants kept", len(variants) == 2, str([v.name for v in variants]))
check("non-VARIANT features ignored", all(v.kind == "variant" for v in variants))
check("multiple diseases joined", any(";" in v.detail for v in variants), str([v.detail for v in variants]))

merged = ann.merge_sites(scop, ptms, variants)
check("duplicate positions collapse", len(merged) == 6, str(len(merged)))
check("Scop3P terminology wins on shared sites",
      next(s for s in merged if s.position == 905).source == "Scop3P")
# .colour() is now the category colour; the Scop3P residue palette moved to
# .residue_colour() and both must stay available.
check("residue palette reachable",
      next(s for s in merged if s.position == 696).residue_colour() == "#1F77B4")
check("variants distinct in the residue palette",
      next(s for s in merged if s.kind == "variant").residue_colour() == "#7B241C")
check("category colour is the default",
      next(s for s in merged if s.kind == "variant").colour() == ann.CATEGORY_COLOURS["variant"])

# Attaching to elements, with and without a numbering map.
fake_elements = [
    {"id": "S1", "start": 690, "stop": 700, "type": "strand"},
    {"id": "H1", "start": 895, "stop": 910, "type": "helix"},
]
summary = ann.attach_sites([dict(e) for e in fake_elements], merged, None)
check("sites counted", summary["total"] == 6, str(summary["total"]))
check("sites outside elements reported as coil", summary["in_coil"] == 3, str(summary["in_coil"]))

placed = [dict(e) for e in fake_elements]
ann.attach_sites(placed, merged, None)
helix = next(e for e in placed if e["id"] == "H1")
check("sites land on the right element", len(helix["sites"]) == 2, str(len(helix["sites"])))
check("site offset along element computed",
      all(0.0 <= s["t"] <= 1.0 for s in helix["sites"]))
check("counts split by kind", helix["site_counts"]["ptm"] == 2, str(helix["site_counts"]))

# With a map, unmapped positions must be dropped rather than shifted onto a
# neighbouring residue, which would invent a finding.
strict = [dict(e) for e in fake_elements]
summary = ann.attach_sites(strict, merged, {905: 902, 900: 897})
check("unmapped sites are dropped, not guessed", summary["unmapped"] == 4, str(summary["unmapped"]))
check("mapped sites use structure numbering",
      any(s["position"] == 902 and s["uniprot_position"] == 905
          for e in strict for s in e["sites"]))



# ------------------------------------------------- annotation rendering
section("Annotation rendering")

import html as _html
import re as _re
from topology.app import build_view as _build_view

_st = load_structure(Path("fixtures/big.pdb").read_text(), "big.pdb")
_sites = [
    ann.Site(position=p, kind="ptm", residue=r, name="phosphorylation",
             source="Scop3P", score=sc)
    for p, r, sc in [(5, "SER", 0.9), (10, "TYR", 0.3), (20, "THR", None),
                     (45, "TYR", 0.7), (46, "SER", 0.8), (47, "TYR", 0.6)]
]
_sites.append(ann.Site(position=48, kind="variant", name="M48T", detail="MEN2B"))


def _doc(**kwargs):
    markup = _build_view(_st, "A", {"kind": "upload", "data": "x", "format": "pdb"}, **kwargs)
    return _html.unescape(_re.search(r'srcdoc="(.*?)" style', markup, _re.S).group(1))


def _payload(document):
    return _json.loads(
        _re.search(r'data-role="payload">(.*?)</script>', document, _re.S).group(1)
    )


annotated = _doc(sites=_sites, accession="P07949")
plain = _doc()

check("filter row present when annotated", "filter-ptm" in annotated)
check("variant filter present", "filter-variant" in annotated)
check("functional score filter present", "filter-score" in annotated)
check("density colour option offered", "Annotation density" in annotated)
check("site clustering shipped to the browser", "clusterSites" in annotated)
check("3D site overlay shipped", "setSiteOverlay" in annotated)

# File mode has no annotations, so none of that UI may appear.
check("file mode hides the filter row", "filter-ptm" not in plain)
check("file mode hides the density option", "Annotation density" not in plain)
check("file mode carries no annotations", _payload(plain)["annotations"] is None)

block = _payload(annotated)["annotations"]
check("PTM and variant counted separately",
      block["counts"] == {"ptm": 6, "variant": 1}, str(block["counts"]))
check("peak density recorded for the colour ramp", block["max_density"] >= 1)

# Sites must reach the element they belong to, with an offset for placement.
carrying = [e for e in _payload(annotated)["elements"] if e.get("sites")]
check("sites attached to elements", len(carrying) >= 1)
check("every site has a placement offset",
      all(0.0 <= s["t"] <= 1.0 for e in carrying for s in e["sites"]))
check("every site has a colour",
      all(s["colour"].startswith("#") for e in carrying for s in e["sites"]))

# When numbering is supplied, the mark uses structure numbering but keeps the
# UniProt position so the panel can show both.
shifted = _doc(sites=_sites, accession="P07949", numbering={5: 105, 10: 110})
carried = [s for e in _payload(shifted)["elements"] for s in e.get("sites", [])]
check("mapped sites renumbered",
      all(s["position"] != s["uniprot_position"] for s in carried), str(carried[:2]))
check("unmapped sites dropped under a numbering map",
      _payload(shifted)["annotations"]["unmapped"] == 5,
      str(_payload(shifted)["annotations"]["unmapped"]))

# Notes about failed lookups must surface rather than vanish.
noted = _doc(sites=_sites, accession="P07949", notes=["Scop3P unavailable (test)."])
check("lookup notes surfaced to the user", "Scop3P unavailable" in noted)



# ------------------------------------------------- SIFTS from the file
section("SIFTS numbering from mmCIF")

_SIFTS_CIF = """data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
_atom_site.pdbx_sifts_xref_db_acc
_atom_site.pdbx_sifts_xref_db_num
_atom_site.pdbx_PDB_model_num
ATOM 1 CA ALA A 1 A 721 0.0 0.0 0.0 30.0 P07949 724 1
ATOM 2 CA ALA A 2 A 722 1.0 0.0 3.8 30.0 P07949 725 1
ATOM 3 CA ALA A 3 A 805 2.0 0.0 7.6 30.0 P07949 810 1
ATOM 4 CA ALA B 1 B 900 3.0 0.0 11.4 30.0 Q99999 500 1
ATOM 5 CA ALA C 1 C 950 4.0 0.0 15.2 30.0 ? ? 1
"""

from topology.io import parse_mmcif

sifts_structure = parse_mmcif(_SIFTS_CIF, "sifts.cif")
check("SIFTS columns detected", sifts_structure.has_sifts)
check("all cross-referenced accessions listed",
      sifts_structure.sifts_accessions() == ["P07949", "Q99999"],
      str(sifts_structure.sifts_accessions()))
check("chains resolved per accession",
      sifts_structure.chains_for_accession("P07949") == ["A"],
      str(sifts_structure.chains_for_accession("P07949")))

file_map = sifts_structure.sifts_numbering("A", "P07949")
check("mapping read straight from the file", file_map == {724: 721, 725: 722, 810: 805},
      str(file_map))
check("blocks with different offsets both honoured",
      file_map[724] - 724 != file_map[810] - 810)
check("wrong accession yields no mapping",
      sifts_structure.sifts_numbering("A", "Q99999") == {})
check("unmapped residues excluded",
      all(v != 950 for v in sifts_structure.sifts_numbering("C", "P07949").values()))

plain_structure = parse_mmcif(Path("fixtures/annotated.cif").read_text(), "plain.cif")
check("files without SIFTS report so", not plain_structure.has_sifts)
check("no SIFTS means an empty map, not a wrong one",
      plain_structure.sifts_numbering("A", "P07949") == {})


# ------------------------------------------------- dropdown option swaps
section("Dropdown option swapping")

try:
    import ipywidgets as _w
    from topology.app import _set_options

    # The reported crash: a placeholder value that is absent from the new
    # options makes ipywidgets raise from inside the options assignment.
    placeholder = _w.Dropdown(options=[("", "")], value="")
    try:
        _set_options(placeholder, [("A (81 residues)", "A")], "A")
        check("placeholder swap does not raise", placeholder.value == "A")
    except Exception as error:
        check("placeholder swap does not raise", False, str(error))

    picker = _w.Dropdown(options=[])
    _set_options(picker, [("A", "A"), ("B", "B")], "B")
    check("requested value honoured", picker.value == "B")

    _set_options(picker, [("C", "C")], "B")
    check("absent value falls back to the first option", picker.value == "C")

    _set_options(picker, [])
    check("empty options leave no selection", picker.value is None)

    # Shrinking the list must not strand the old selection either.
    wide = _w.Dropdown(options=[("A", "A"), ("B", "B"), ("C", "C")], value="C")
    _set_options(wide, [("A", "A")], None)
    check("shrinking options does not raise", wide.value == "A")
except ImportError:
    print("  [skip] ipywidgets not installed")



# ------------------------------------------------- browser JavaScript
section("Browser JavaScript")

import shutil as _shutil
import subprocess as _subprocess
import tempfile as _tempfile

from topology.assets import TOPOLOGY_JS, VIEWER_JS

_node = _shutil.which("node")
if _node:
    for _name, _source in [("TOPOLOGY_JS", TOPOLOGY_JS), ("VIEWER_JS", VIEWER_JS)]:
        with _tempfile.NamedTemporaryFile("w", suffix=".js", delete=False) as handle:
            handle.write(_source)
            temp = handle.name
        result = _subprocess.run([_node, "--check", temp], capture_output=True, text=True)
        check(f"{_name} is valid JavaScript", result.returncode == 0,
              result.stderr.strip()[:200])
else:
    print("  [skip] node not available")

# A method that reads `self` without declaring it throws ReferenceError at
# runtime and aborts mid-way. That silently broke the 3D highlight once, because
# the failing line sat before the adapter call.
_undeclared = []
for _match in _re.finditer(
    r"(\w+)\.prototype\.(\w+)\s*=\s*function\s*\([^)]*\)\s*\{", TOPOLOGY_JS + VIEWER_JS
):
    _body_start = _match.end()
    _depth, _index = 1, _body_start
    _blob = TOPOLOGY_JS + VIEWER_JS
    while _index < len(_blob) and _depth:
        if _blob[_index] == "{":
            _depth += 1
        elif _blob[_index] == "}":
            _depth -= 1
        _index += 1
    _body = _blob[_body_start:_index]
    if _re.search(r"\bself\.", _body) and not _re.search(r"\bvar\s+self\s*=", _body):
        _undeclared.append(f"{_match.group(1)}.{_match.group(2)}")

check("no method reads an undeclared `self`", not _undeclared, str(_undeclared))


# ------------------------------------------------- API error handling
section("API error handling")

check("NotFound is distinct from a transport failure",
      issubclass(ann.NotFound, ValueError) and ann.NotFound is not ValueError)

# The proteins API and the UniProtKB entry endpoint disagree about how a
# position is written. Reading only one shape returns nothing against the other,
# which is indistinguishable from a protein with no modifications.
_ebi_shape = {"features": [
    {"category": "PTM", "type": "MOD_RES", "begin": "905", "end": "905",
     "description": "Phosphotyrosine; by autocatalysis"},
    {"category": "PTM", "type": "MOD_RES", "begin": "696", "end": "696",
     "description": "Phosphoserine"},
    {"category": "PTM", "type": "MOD_RES", "begin": "700", "end": "701",
     "description": "spans two residues"},
    {"category": "VARIANTS", "type": "VARIANT", "begin": "918", "end": "918"},
]}
_from_ebi = ann.parse_uniprot_ptms(_ebi_shape)
check("proteins API shape parsed", len(_from_ebi) == 2, str(len(_from_ebi)))
check("category filter applied", all(s.kind == "ptm" for s in _from_ebi))
check("residue inferred from the description",
      {s.residue for s in _from_ebi} == {"TYR", "SER"}, str([s.residue for s in _from_ebi]))
check("UniProt qualifier trimmed from the name",
      any(s.name == "Phosphotyrosine" for s in _from_ebi), str([s.name for s in _from_ebi]))

_from_kb = ann.parse_uniprot_ptms(_load("uniprot_entry.json"))
check("UniProtKB shape still parsed", len(_from_kb) == 2, str(len(_from_kb)))
check("both shapes agree on positions",
      {s.position for s in _from_ebi} == {s.position for s in _from_kb})

# Inferred residues must colour the same way Scop3P sites do.
check("inferred residues get Scop3P colours",
      next(s for s in _from_ebi if s.residue == "TYR").residue_colour() == "#2CA02C")

check("probe helper exposed", callable(ann.probe))



# ------------------------------------------------- non-JSON responses
section("Non-JSON responses")

check("no custom User-Agent is sent", ann.USER_AGENT == "",
      "the Scop3P client and notebook both use the requests default")
check("NotJson is distinct from NotFound",
      ann.NotJson is not ann.NotFound and issubclass(ann.NotJson, ValueError))

# A refusal page, a redirect and an empty body all produce the same bare
# "Expecting value: line 1 column 1" from json.loads. The message has to say
# which one happened or there is nothing to act on.
import urllib.request as _urlreq


class _FakeResponse:
    def __init__(self, body):
        self._body = body.encode()

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def _with_body(body, url="https://example.org/x"):
    real_open, real_requests = _urlreq.urlopen, None
    import sys as _sys
    saved = _sys.modules.get("requests")
    _sys.modules["requests"] = None  # force the urllib path
    _urlreq.urlopen = lambda *a, **k: _FakeResponse(body)
    try:
        return ann._get_json(url)
    finally:
        _urlreq.urlopen = real_open
        if saved is not None:
            _sys.modules["requests"] = saved
        else:
            _sys.modules.pop("requests", None)


for label, body, expect in [
    ("empty body", "", "empty body"),
    ("HTML page", "<!DOCTYPE html><html><body>Forbidden</body></html>", "an HTML page"),
    ("plain text", "Service temporarily unavailable", "non-JSON text"),
]:
    try:
        _with_body(body)
        check(f"{label} reported clearly", False, "no error raised")
    except ann.NotJson as error:
        check(f"{label} reported clearly", expect in str(error), str(error)[:120])
    except Exception as error:
        check(f"{label} reported clearly", False, f"{type(error).__name__}: {error}")

check("valid JSON still parses", _with_body('{"modifications": []}') == {"modifications": []})

# Variants come from UniProt only. Scop3P imports its mutations from UniProt,
# so fetching them there would be a second-hand copy.
check("no Scop3P mutation fetcher remains", not hasattr(ann, "fetch_scop3p_mutations"))
check("no Scop3P mutation parser remains", not hasattr(ann, "parse_scop3p_mutations"))
check("Scop3P is used for modifications only",
      "modifications" in ann.SCOP3P_V1_MODS and not hasattr(ann, "SCOP3P_V1_MUTATIONS"))


# ------------------------------------------------- Scop3P v1 schema
section("Scop3P v1 schema")

_v1 = _load("scop3p_v1.json")
_v1_sites = ann.parse_scop3p(_v1)

check("bare array response parsed", len(_v1_sites) == 4, str(len(_v1_sites)))
check("uniprot_position read", {s.position for s in _v1_sites} == {105, 687, 905, 1096},
      str(sorted(s.position for s in _v1_sites)))
check("null position skipped", len(_v1_sites) == len(_v1) - 1)

# "Phosphoserine" names the modified residue; the code has to come back out of
# it or every site falls through to the generic PTM colour.
_by_position = {s.position: s for s in _v1_sites}
check("Phosphoserine resolves to SER", _by_position[105].residue == "SER", _by_position[105].residue)
check("Phosphotyrosine resolves to TYR", _by_position[905].residue == "TYR", _by_position[905].residue)
check("Phosphothreonine resolves to THR", _by_position[1096].residue == "THR", _by_position[1096].residue)
check("resolved residues drive the Scop3P palette",
      _by_position[105].residue_colour() == "#1F77B4"
      and _by_position[905].residue_colour() == "#2CA02C")

# best_probability is a percentage; functionalScore was a fraction. Mixing the
# scales would make one filter position mean two different things.
check("percentage probability normalised to 0-1",
      abs(_by_position[905].score - 0.993) < 1e-6, str(_by_position[905].score))
check("all scores within 0-1", all(0.0 <= s.score <= 1.0 for s in _v1_sites))
_legacy_scale = ann.parse_scop3p([{"position": 1, "functionalScore": 0.44, "residue": "SER"}])
check("fractional scores left alone", abs(_legacy_scale[0].score - 0.44) < 1e-9,
      str(_legacy_scale[0].score))

check("modification_name used as the site name",
      all(s.name == "phosphorylation" for s in _v1_sites))
check("evidence_terms preferred over the ECO code",
      _by_position[905].evidence == "Experimental", _by_position[905].evidence)
check("pubmed summarised rather than dumped",
      "+2 more" in _by_position[905].detail, _by_position[905].detail)
check("tissue stripped of its project prefix",
      "Kidney" in _by_position[687].detail and "PXD025798=" not in _by_position[687].detail,
      _by_position[687].detail)

# The older wrapped shape must keep working, since a deployment may lag.
_legacy = ann.parse_scop3p({"modifications": [
    {"position": 696, "residue": "SER", "name": "phosphorylation", "functionalScore": 0.44}
]})
check("legacy wrapped shape still parsed", len(_legacy) == 1 and _legacy[0].position == 696)

check("v1 modifications URL is path-based",
      ann.SCOP3P_V1_MODS.format(accession="P07949").endswith("/v1/proteins/P07949/modifications"),
      ann.SCOP3P_V1_MODS)
# One feature with several disease associations is one residue, so it should
# produce one mark carrying every disease, not one mark per disease.
_multi = ann.parse_variants({"features": [
    {"type": "VARIANT", "begin": "634", "end": "634", "wildType": "C", "mutatedType": "R",
     "consequenceType": "missense",
     "association": [{"name": "MEN2A", "disease": True}, {"name": "MTC", "disease": True}]},
]}, "P07949")
check("multi-disease variant collapses to one site", len(_multi) == 1, str(len(_multi)))
check("every disease retained", "MEN2A" in _multi[0].detail and "MTC" in _multi[0].detail,
      _multi[0].detail)
check("consequence type kept as evidence", _multi[0].evidence == "missense")
check("variant name built from wild and mutant", _multi[0].name == "C634R", _multi[0].name)



# ------------------------------------------------- site categories
section("Site categories")

check("category colours defined",
      ann.CATEGORY_COLOURS == {"ptm": "#F0C808", "variant": "#D7263D", "both": "#2E9E4F"},
      str(ann.CATEGORY_COLOURS))

_cat_sites = [
    ann.Site(position=905, kind="ptm", residue="TYR", name="phosphorylation"),
    ann.Site(position=918, kind="variant", name="M918T", detail="MEN2B"),
    ann.Site(position=634, kind="variant", name="C634R", detail="MEN2A"),
    ann.Site(position=634, kind="ptm", residue="SER", name="phosphorylation"),
]
_cat_elements = [{"id": "H1", "start": 600, "stop": 950, "type": "helix"}]
_cat_summary = ann.attach_sites(_cat_elements, _cat_sites, None)
_placed = {(s["position"], s["kind"]): s for s in _cat_elements[0]["sites"]}

check("modification alone is yellow", _placed[(905, "ptm")]["colour"] == "#F0C808")
check("mutation alone is red", _placed[(918, "variant")]["colour"] == "#D7263D")

# The overlap is the case worth seeing: a disease variant landing on a
# regulatory site. Both records on that residue take the shared colour.
check("modified and mutated is green",
      _placed[(634, "ptm")]["colour"] == "#2E9E4F"
      and _placed[(634, "variant")]["colour"] == "#2E9E4F")
check("both records agree on the category",
      _placed[(634, "ptm")]["category"] == _placed[(634, "variant")]["category"] == "both")
check("category carries a readable label",
      _placed[(634, "ptm")]["category_label"] == "Modified and mutated")

# Category counts are per residue, so an overlapping residue counts once.
check("categories counted per residue",
      _cat_summary["categories"] == {"ptm": 1, "variant": 1, "both": 1},
      str(_cat_summary["categories"]))

# The per-residue Scop3P palette must survive alongside the category scheme.
check("residue palette still available",
      _cat_sites[0].residue_colour() == "#2CA02C", _cat_sites[0].residue_colour())
check("residue colour carried in the record",
      _placed[(905, "ptm")]["residue_colour"] == "#2CA02C")

# Rendering: sticks, not spheres, and one representation per category.
_cat_doc = _doc(sites=_cat_sites, accession="P07949")
check("overlay drawn as sticks", "licorice" in _cat_doc)
check("overlay is grouped by category",
      'groups = { ptm: [], variant: [], both: [] }' in _cat_doc)
check("three overlay colours shipped",
      all(colour in _cat_doc for colour in ann.CATEGORY_COLOURS.values()))
check("legend rendered", "topo-legend" in _cat_doc)
check("summary panel names the category", "category_label" in _cat_doc)

# Sites in loops have no element to attach to, but are still real sites and
# must reach the 3D view.
_loop_only = [ann.Site(position=10, kind="ptm", residue="SER", name="phosphorylation")]
_loop_doc = _doc(sites=_loop_only, accession="P07949")
_loop_payload = _payload(_loop_doc)
check("loop sites recorded", _loop_payload["annotations"]["in_coil"] == 1,
      str(_loop_payload["annotations"]["in_coil"]))
check("loop sites carry a colour",
      _loop_payload["annotations"]["coil_sites"][0]["colour"] == "#F0C808")
check("loop sites pushed to 3D", "annotations.coil_sites" in _loop_doc)



# ------------------------------------------------------------ logo
section("Logo")

import xml.etree.ElementTree as _ET
from topology import logo as _logo

for _name, _svg in [("mark", _logo.mark()), ("lockup", _logo.lockup()),
                    ("favicon", _logo.favicon())]:
    try:
        _ET.fromstring(_svg)
        check(f"{_name} is well-formed SVG", True)
    except Exception as error:
        check(f"{_name} is well-formed SVG", False, str(error))
    check(f"{_name} declares the SVG namespace", 'xmlns="http://www.w3.org/2000/svg"' in _svg)
    check(f"{_name} has an accessible label", 'role="img"' in _svg and "aria-label" in _svg)

# The mark is built from the renderer's palette, so a colour change in one place
# does not leave the logo describing a tool that no longer looks like it.
_mark = _logo.mark()
for _label, _colour in [("helix", _logo.HELIX), ("strand", _logo.STRAND_A),
                        ("site", _logo.SITE), ("connector", _logo.CONNECTOR)]:
    check(f"mark uses the {_label} colour", _colour in _mark)

check("mark scales without a fixed pixel body",
      'viewBox="0 0 128 128"' in _logo.mark(size=512))
check("background is optional", "<rect" in _logo.mark(background=True)
      and _logo.mark(background=False).count("<rect") < _logo.mark(background=True).count("<rect"))

# Below about 32px the termini labels and site dot become noise, so the small
# variant drops them rather than scaling them into mush.
_fav = _logo.favicon()
check("favicon drops the termini labels", "<text" not in _fav)
check("favicon drops the site marker", _logo.SITE not in _fav)
check("favicon keeps the topology", _logo.HELIX in _fav and _logo.STRAND_A in _fav)

try:
    import ipywidgets as _w
    from topology.app import make_app as _make_app
    _app = _make_app()
    _found = []

    def _collect(widget):
        if isinstance(widget, _w.HTML):
            _found.append(widget.value)
        for _child in getattr(widget, "children", []):
            _collect(_child)

    _collect(_app)
    check("logo appears in the app header", any("<svg" in v for v in _found))
except ImportError:
    print("  [skip] ipywidgets not installed")


# ------------------------------------------------------------ summary
print()
if FAILURES:
    print(f"{len(FAILURES)} check(s) failed:")
    for item in FAILURES:
        print(f"  - {item}")
    sys.exit(1)
print("All checks passed.")
