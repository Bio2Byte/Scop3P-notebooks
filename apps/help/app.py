"""The Help page: what each protocol is for, and which one to reach for.

Static content, so there are no reactive outputs at all. That is deliberate -- a page
with no outputs cannot hit the hidden-output suspension trap described in
apps/README.md, and it renders identically whether it is opened directly or through the
portal navbar.

Kept in the toolkit rather than in the docs because the question it answers -- "which of
these five tools do I want?" -- arrives while someone is looking at the navbar, not while
they are reading the repository.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

from shiny import App, ui

_APPS_DIR = Path(__file__).resolve().parents[1]
if str(_APPS_DIR) not in sys.path:
    sys.path.append(str(_APPS_DIR))

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from common.vendor import enable_compression, static_assets  # noqa: E402
from common.logging_utils import get_logger, new_trail  # noqa: E402
from common.ui_shell import scop3p_card, scop3p_footer, scop3p_shell  # noqa: E402


LOGGER = get_logger("scop3p.help")


@dataclass(frozen=True, slots=True)
class Protocol:
    key: str            # the ?app= value, so the card can link to the tool
    name: str
    question: str       # the question a user arrives with
    mission: str
    scope: tuple[str, ...]
    use_cases: tuple[str, ...]
    inputs: str
    spec: str           # path to the parity spec, for the curious


PROTOCOLS: tuple[Protocol, ...] = (
    Protocol(
        key="structure-viz",
        name="Structure Visualisation",
        question="What is happening around this residue, structurally and biophysically?",
        mission=(
            "Bring every layer of evidence about one protein into a single "
            "structure-centred workspace: modifications, disease variants, predicted "
            "biophysical properties, contact networks and structural superposition."
        ),
        scope=(
            "Six tabs, each a separate question, sharing one accession set at the top.",
            "Structures can be experimental PDB entries, AlphaFold models, or your own "
            "uploaded coordinates.",
            "Bio2Byte predictions (backbone dynamics, disorder, early folding) are "
            "computed on the sequence and painted onto the structure.",
            "TM-align runs server-side to superpose two structures and report the "
            "aligned region.",
        ),
        use_cases=(
            "Ask whether the PTMs of a protein cluster in ordered or disordered regions.",
            "Check whether a disease variant sits near a modification site in 3D, rather "
            "than merely nearby in sequence.",
            "Colour a structure by predicted flexibility or early-folding propensity to "
            "see whether a site sits in a rigid core or a mobile loop.",
            "Reduce a structure to a residue interaction network and inspect which "
            "residues are the most connected.",
            "Superpose two structures and read off how much of the fold is shared.",
        ),
        inputs="A UniProtKB accession; optionally PDB IDs or uploaded coordinate files.",
        spec="docs/use-cases/structure_viz.md",
    ),
    Protocol(
        key="rinalign",
        name="RIN Alignment",
        question="How do the contact networks of two structures differ?",
        mission=(
            "Treat a structure as a graph of residue contacts, then compare two such "
            "graphs -- either two models of the same protein, or two different proteins "
            "-- so structural differences become countable rather than impressionistic."
        ),
        scope=(
            "Same protein (diff) matches residue positions directly and splits contacts "
            "into conserved, lost and gained, with a per-residue impact table. Only "
            "positions present in both structures are compared, so an unresolved region "
            "never masquerades as a lost contact.",
            "Different proteins (align) solves a graph alignment: Weisfeiler-Lehman "
            "neighbourhood signatures combined with residue type and degree, resolved to "
            "a one-to-one mapping by Hungarian assignment.",
            "Contacts are CB-CB (CA for glycine) within a cutoff you choose, 4 to 14 A, "
            "excluding sequence neighbours. Phospho-residues are treated as real "
            "residues, not heteroatoms.",
            "PTM and disease-variant positions can be overlaid on the network views.",
        ),
        use_cases=(
            "Compare an AlphaFold model against an experimental structure of the same "
            "protein and see which contacts the model gets wrong.",
            "Compare apo and holo forms, or two conformational states, and quantify what "
            "rearranged.",
            "Find the residues whose contact count changes most between two structures.",
            "Ask whether two distantly related proteins share a contact topology even "
            "where their sequences do not align.",
            "Click a residue in the network and see it highlighted in 3D, and the reverse.",
        ),
        inputs="A UniProtKB accession, then two structures chosen from its AlphaFold model and PDB entries.",
        spec="docs/use-cases/rinalign.md",
    ),
    Protocol(
        key="mutation-effect",
        name="Mutation Effect",
        question="What does this mutation do to the protein's biophysical behaviour?",
        mission=(
            "Predict biophysical properties for a wild-type sequence and for a mutated "
            "version of it, then report where the two disagree -- so a substitution can "
            "be read as a change in behaviour rather than just a change in letter."
        ),
        scope=(
            "Bio2Byte predictors: DynaMine backbone dynamics, DisoMine disorder and "
            "EFoldMine early folding.",
            "Sequence-based only. Nothing here uses a structure, which is what lets it "
            "work on any protein regardless of structural coverage.",
            "The inference step converts continuous predictions into categorical labels "
            "and reports shifts at the mutated position and within a five-residue window.",
            "Scop3P PTMs are kept on the plots throughout, so a shift can be read "
            "relative to known modification sites.",
        ),
        use_cases=(
            "Ask whether a disease mutation is predicted to rigidify or loosen its "
            "surroundings.",
            "Check whether a substitution is predicted to push a region across the "
            "order/disorder boundary.",
            "Compare several candidate mutations at the same position in one session.",
            "See whether a predicted change lands on or near a known phosphosite.",
        ),
        inputs="A UniProtKB accession, plus one or more 1-indexed positions and target amino acids.",
        spec="docs/use-cases/mutation_effect.md",
    ),
    Protocol(
        key="peptide-mapper",
        name="Peptide Mapper",
        question="Where do my phospho-peptides sit on the structure?",
        mission=(
            "Turn a list of phospho-peptides into a picture: which parts of the folded "
            "protein your mass-spectrometry evidence actually covers, and where the "
            "modified residues fall within that coverage."
        ),
        scope=(
            "Peptides come either from Scop3P for a given accession, or from your own "
            "search-engine export (TSV/CSV) via the Upload tab.",
            "Structures are AlphaFold models only; there is no experimental-PDB path "
            "here. Use Structure Visualisation for that.",
            "Colour grammar: grey is the whole protein, blue is the union of the "
            "peptides you selected, red is their intersection, magenta marks "
            "modified sites.",
        ),
        use_cases=(
            "Check whether a phosphosite of interest is inside observed peptide coverage "
            "at all, or is only inferred.",
            "See whether several peptides agree on the same region, by looking at the red "
            "intersection.",
            "Bring an in-house search result and inspect coverage per protein without "
            "writing any mapping code.",
            "Export a self-contained HTML session, the PDB, or a TSV of mapped residues "
            "to share with a collaborator.",
        ),
        inputs="A UniProtKB accession, or a peptide table with protein, sequence, start, end and position columns.",
        spec="docs/use-cases/peptide_mapper.md",
    ),
    Protocol(
        key="topology-viewer",
        name="Topology Viewer",
        question="What is the fold's wiring, and where do my sites land on it?",
        mission=(
            "Flatten a three-dimensional fold into a topology diagram -- helices, "
            "strands, sheets and the loops connecting them -- so the architecture is "
            "readable at a glance, with the 3D structure beside it for reference."
        ),
        scope=(
            "Two strictly separate modes. Accession mode fetches an AlphaFold model or a "
            "PDBe entry and can overlay annotations. File mode takes your uploaded "
            "structure, makes no network requests, and deliberately offers no overlay: a "
            "local prediction has no reliable UniProt numbering.",
            "Secondary structure is read from the file when present and derived from "
            "coordinates when it is not. The provenance is always displayed, because "
            "different methods disagree about element boundaries by a residue or two.",
            "Annotations are positioned through the SIFTS numbering map. If no map is "
            "available the sites are hidden rather than drawn at the wrong residues.",
        ),
        use_cases=(
            "Understand the domain architecture of an unfamiliar protein faster than by "
            "rotating a 3D model.",
            "Check whether phosphosites fall on structured elements or in the loops "
            "between them.",
            "Inspect the topology of a locally predicted structure (ColabFold, Boltz, "
            "Chai) fully offline.",
            "Compare which secondary-structure assignment method your structure supports, "
            "using the provenance line.",
        ),
        inputs="A UniProtKB accession, or an uploaded .pdb / .cif / .ent / .mmcif file.",
        spec="docs/use-cases/topology_viewer.md",
    ),
)


def _chooser() -> ui.Tag:
    """The question-first table. Most people arrive with a question, not a tool name."""
    rows = [
        ui.tags.tr(
            ui.tags.td(ui.tags.em(protocol.question)),
            ui.tags.td(
                ui.a(protocol.name, href=f"/?app={protocol.key}", class_="help-tool-link")
            ),
        )
        for protocol in PROTOCOLS
    ]
    return ui.TagList(
        ui.p(
            "Five protocols, each answering a different question about one protein. "
            "Start from the question:",
            class_="scop3p-note",
        ),
        ui.tags.table(
            ui.tags.thead(
                ui.tags.tr(ui.tags.th("If you are asking..."), ui.tags.th("Open"))
            ),
            ui.tags.tbody(*rows),
            class_="help-table",
        ),
    )


def _protocol_card(protocol: Protocol) -> ui.Tag:
    return scop3p_card(
        protocol.name,
        ui.p(protocol.question, class_="help-question"),
        ui.h5("Mission", class_="help-heading"),
        ui.p(protocol.mission),
        ui.h5("Scope", class_="help-heading"),
        ui.tags.ul(*[ui.tags.li(item) for item in protocol.scope]),
        ui.h5("Use cases", class_="help-heading"),
        ui.tags.ul(*[ui.tags.li(item) for item in protocol.use_cases]),
        ui.h5("What you need", class_="help-heading"),
        ui.p(protocol.inputs),
        ui.div(
            ui.a(
                f"Open {protocol.name}",
                href=f"/?app={protocol.key}",
                class_="btn btn-primary help-open-btn",
            ),
            ui.span(
                "Parity spec: ",
                ui.tags.code(protocol.spec),
                class_="scop3p-note help-spec",
            ),
            class_="help-card-footer",
        ),
        extra_class="help-card",
    )


app_ui = scop3p_shell(
    "Help",
    "What each protocol in the toolkit is for, what it covers, and which one answers the "
    "question you arrived with.",
    ui.tags.style(
        """
.help-table {
  width: 100%;
  border-collapse: collapse;
  margin-top: 6px;
}
.help-table th {
  text-align: left;
  padding: 8px 10px;
  border-bottom: 2px solid var(--scop3p-line);
  font-size: 0.85rem;
  text-transform: uppercase;
  letter-spacing: 0.06em;
  color: var(--scop3p-muted);
}
.help-table td {
  padding: 9px 10px;
  border-bottom: 1px solid var(--scop3p-line);
  vertical-align: middle;
}
.help-table td:last-child { width: 210px; }
.help-tool-link {
  font-weight: 700;
  color: var(--scop3p-accent);
  text-decoration: none;
}
.help-tool-link:hover { text-decoration: underline; }
.help-question {
  font-size: 1.02rem;
  font-style: italic;
  color: var(--scop3p-accent);
  margin: -4px 0 14px;
}
.help-heading {
  margin: 16px 0 6px;
  font-size: 0.82rem;
  text-transform: uppercase;
  letter-spacing: 0.08em;
  color: var(--scop3p-muted);
}
.help-card ul { margin: 0; padding-left: 20px; }
.help-card li { margin-bottom: 6px; }
.help-card-footer {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 14px;
  margin-top: 18px;
  padding-top: 14px;
  border-top: 1px solid var(--scop3p-line);
}
.help-open-btn { color: #fff; text-decoration: none; }
.help-open-btn:hover { color: #fff; }
.help-spec { margin: 0; }
.help-stack {
  display: grid;
  gap: 18px;
}
"""
    ),
    scop3p_card("Which tool do I need?", _chooser()),
    ui.div(*[_protocol_card(protocol) for protocol in PROTOCOLS], class_="help-stack"),
    scop3p_card(
        "Common ground",
        ui.tags.ul(
            ui.tags.li(
                "Every protocol identifies a protein the same way, by "
                "UniProtKB accession, and every accession field has a "
                "Load example button if you just want to see the tool work."
            ),
            ui.tags.li(
                "Scop3P is the source of experimentally observed modifications; UniProt "
                "and the EBI Proteins API supply sequences, PTM features and disease "
                "variants; AlphaFold DB, PDBe and RCSB supply structures."
            ),
            ui.tags.li(
                "3D viewers load from public CDNs at page load, so the browser needs "
                "outbound network access even when the server has the data."
            ),
            ui.tags.li(
                "Each session writes a log and a metadata.yml recording the exact "
                "dependency and tool versions used, so a result can be traced back to "
                "the environment that produced it."
            ),
        ),
    ),
)


def server(input, output, session):  # noqa: ARG001 - static page, no reactivity
    # Static page, but consulting Help is part of the session narrative: "read Help,
    # then opened Structure Visualisation" is exactly the kind of thing the record is
    # for. One trail per session, same as every other protocol.
    trail = new_trail()
    trail.opened("Help")


content_ui = ui.div(
    app_ui, scop3p_footer()
)

# static_assets serves the vendored browser libraries; every app mounts the same prefix,
# so /vendor/... resolves whichever app the portal is serving. enable_compression is not
# optional cosmetics: Shiny sends static files raw, and molstar.js is 5 MB uncompressed
# against 1.45 MB gzipped, so without it vendoring would put more bytes on the wire.
app = App(content_ui, server, static_assets=static_assets())
enable_compression(app)
