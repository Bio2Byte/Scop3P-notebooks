"""The renderer.

One renderer, replacing the three near-duplicate variants in the original
module. It draws whatever the layout engines emit, so adding a layout does not
mean adding a renderer, and the annotation layer that lands in phase 3 attaches
here once rather than three times.

The 3D side sits behind an adapter interface so Mol* and NGL are
interchangeable::

    init(node)                       prepare the stage
    loadStructure({url, data, fmt})  fetched entry or uploaded bytes
    highlight(chain, seq)            transient, follows the pointer
    select(chain, seq)               sticky, follows a click
    clearHighlight()
    setScope(mode, chain)            "chain" or "complex"
    dispose()

Residue lookup is keyed on ``chain:seq`` rather than ``seq`` alone, so a
homodimer numbered 1-250 twice highlights the chain the diagram is actually
describing.
"""

from __future__ import annotations

import html as html_module
import json
import uuid
from typing import Any, Dict

from .assets import VIEWER_JS, TOPOLOGY_CSS, TOPOLOGY_JS


def _escape(value: Any) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def build_payload(
    structure: Any,
    chain: str,
    elements: list,
    residues: list,
    contacts: list,
    layouts: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    """Assemble the frozen payload the browser consumes.

    ``annotations`` stays ``None`` for uploaded files, which is what switches
    the renderer out of accession mode: no site marks, no density channel, no
    filter row.
    """
    confidences = [
        residue.plddt for residue in residues if residue.plddt is not None
    ]
    return {
        "schema": 1,
        "name": structure.name,
        "chain": chain,
        "chains": structure.chain_options(),
        "layouts": layouts,
        "elements": elements,
        "contacts": contacts,
        "residues": [
            {
                "seq": residue.seq,
                "label_seq": residue.label_seq,
                "chain": residue.chain,
                "label_chain": residue.label_chain,
                "aa": residue.aa,
                "comp": residue.comp_id,
                "confidence": residue.plddt,
            }
            for residue in residues
        ],
        "source": {
            "format": structure.fmt,
            "entry_id": structure.entry_id,
            "title": structure.title,
            "uniprot": structure.uniprot,
            "ss_source": structure.ss_source,
            "has_confidence": structure.has_plddt,
            "confidence_label": "pLDDT" if structure.has_plddt else "B-factor",
            "confidence_range": (
                [min(confidences), max(confidences)] if confidences else None
            ),
        },
        "structure_source": {},
        "annotations": None,
        "stats": {
            "residues": len(residues),
            "helices": sum(1 for e in elements if e["type"] == "helix"),
            "strands": sum(1 for e in elements if e["type"] == "strand"),
            "contacts": len(contacts),
        },
    }


def render(payload: Dict[str, Any], height: int = 1150, embed: str = "iframe") -> str:
    """Produce HTML for one topology view.

    ``embed="iframe"`` wraps the view in a sandboxed iframe via ``srcdoc``.
    That is not decoration: JupyterLab does not execute ``<script>`` tags inside
    ``display(HTML(...))`` output, so an inline block renders the control bar
    and then sits there with an empty diagram. An iframe document runs its own
    scripts and works the same way under JupyterLab, Voila and nbconvert.
    Pass ``embed="inline"`` to get the bare block, which is what the standalone
    HTML export wants.
    """
    block = _render_block(payload)
    if embed != "iframe":
        return block

    document = (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        "<style>html,body{margin:0;padding:12px;background:#fff;}</style>"
        "</head><body>" + block + "</body></html>"
    )
    srcdoc = html_module.escape(document, quote=True)
    return (
        f'<iframe srcdoc="{srcdoc}" '
        f'style="width:100%;height:{height}px;border:1px solid #d3dce6;border-radius:8px;" '
        'sandbox="allow-scripts allow-downloads allow-popups allow-same-origin" '
        'loading="lazy" title="Protein topology viewer"></iframe>'
    )


def standalone_document(payload: Dict[str, Any]) -> str:
    """A complete, self-contained HTML page for saving or sharing."""
    name = _escape(payload.get("source", {}).get("entry_id") or payload.get("name") or "topology")
    return (
        "<!doctype html><html><head><meta charset='utf-8'>"
        f"<title>Topology {name}</title>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        "<style>body{margin:0;padding:16px;background:#fff;}</style>"
        "</head><body>" + _render_block(payload) + "</body></html>"
    )


def _render_block(payload: Dict[str, Any]) -> str:
    """The view itself, without any page or frame around it."""
    root_id = "topo-" + uuid.uuid4().hex[:12]
    blob = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    source = payload.get("source", {})
    stats = payload.get("stats", {})

    subtitle_bits = [
        _escape(source.get("entry_id") or payload.get("name") or "structure"),
        f"chain {_escape(payload.get('chain'))}",
        f"{stats.get('residues', 0)} residues",
        f"{stats.get('helices', 0)} helices",
        f"{stats.get('strands', 0)} strands",
    ]

    has_confidence = bool(source.get("has_confidence"))
    confidence_label = _escape(source.get("confidence_label") or "Confidence")
    confidence_option = (
        f'<option value="confidence">{confidence_label}</option>' if has_confidence else ""
    )

    annotations = payload.get("annotations") or None
    density_option = (
        '<option value="density">Annotation density</option>' if annotations else ""
    )

    filter_row = ""
    if annotations:
        counts = annotations.get("counts", {})
        notes = annotations.get("notes") or []
        caveats = ""
        if annotations.get("in_coil"):
            caveats += (
                f" &middot; {annotations['in_coil']} in loops "
                "(not on a helix or strand)"
            )
        if annotations.get("unmapped"):
            caveats += (
                f" &middot; {annotations['unmapped']} not present in this structure"
            )
        note_html = "".join(
            f"<div class='muted'>{_escape(note)}</div>" for note in notes
        )

        # Legend. The "both" category only earns a place when it occurs.
        categories = annotations.get("categories", {})
        colours = annotations.get("category_colours", {})
        labels = annotations.get("category_labels", {})
        legend_bits = []
        for key in ("ptm", "variant", "both"):
            if not categories.get(key):
                continue
            legend_bits.append(
                f"<span class='legend-item'>"
                f"<span class='site-chip' style='background:{_escape(colours.get(key, '#999'))}'></span>"
                f"{_escape(labels.get(key, key))} ({categories[key]})</span>"
            )
        legend_html = (
            f"<div class='topo-legend'>{''.join(legend_bits)}</div>" if legend_bits else ""
        )
        filter_row = f"""
  <div class="topo-filters">
    <label title="Modifications from Scop3P, merged with UniProt PTM features">
      <input type="checkbox" data-role="filter-ptm" checked>
      PTMs ({counts.get('ptm', 0)})</label>
    <label title="Disease-associated variants from UniProt">
      <input type="checkbox" data-role="filter-variant" checked>
      Disease variants ({counts.get('variant', 0)})</label>
    <label>Min site probability
      <input type="range" data-role="filter-score" min="0" max="1" step="0.05" value="0">
      <span data-role="score-readout">0%</span></label>
    <span class="muted">{_escape(annotations.get('accession', ''))}{caveats}</span>
    {legend_html}
    {note_html}
  </div>"""

    chain_options = "".join(
        f'<option value="{_escape(value)}"'
        f'{" selected" if value == payload.get("chain") else ""}>{_escape(label)}</option>'
        for label, value in payload.get("chains", [])
    )
    multi_chain = len(payload.get("chains", [])) > 1
    chain_row_style = "" if multi_chain else ' style="display:none"'

    css = TOPOLOGY_CSS.replace("__ROOT__", f"#{root_id}")

    return f"""
<div id="{root_id}" class="topo-root">
  <style>{css}</style>

  <header class="topo-head">
    <div class="topo-title">{" &middot; ".join(subtitle_bits)}</div>
    <div class="topo-provenance">Secondary structure: {_escape(source.get('ss_source') or 'unknown')}</div>
  </header>

  <div class="topo-controls">
    <label class="topo-field"{chain_row_style}>
      <span>Chain</span>
      <select data-role="chain">{chain_options}</select>
    </label>

    <label class="topo-field">
      <span>Layout</span>
      <select data-role="layout">
        <option value="sheet">Sheet topology</option>
        <option value="serpentine">Sequence order</option>
        <option value="spatial">Spatial arrangement</option>
      </select>
    </label>

    <label class="topo-field">
      <span>Colour by</span>
      <select data-role="colour">
        <option value="ss">Structure type</option>
        {confidence_option}
        {density_option}
      </select>
    </label>

    <label class="topo-field"{chain_row_style}>
      <span>3D scope</span>
      <select data-role="scope">
        <option value="complex">Complex, chain emphasised</option>
        <option value="chain">Selected chain only</option>
      </select>
    </label>

    <label class="topo-field">
      <span>3D engine</span>
      <select data-role="engine">
        <option value="ngl">NGL</option>
        <option value="molstar">Mol*</option>
      </select>
    </label>

    <div class="topo-actions">
      <button type="button" data-action="zoom-in" title="Zoom in">+</button>
      <button type="button" data-action="zoom-out" title="Zoom out">&minus;</button>
      <button type="button" data-action="reset" title="Fit to view">Fit</button>
      <input type="text" data-role="jump" placeholder="Residue" size="7">
      <button type="button" data-action="jump">Go</button>
      <button type="button" data-action="download">Save JSON</button>
    </div>
  </div>

  {filter_row}

  <div class="topo-grid">
    <section class="topo-panel topo-2d">
      <svg data-role="svg" xmlns="http://www.w3.org/2000/svg" role="img"
           aria-label="Secondary structure topology diagram">
        <g data-role="viewport"></g>
      </svg>
      <div class="topo-tooltip" data-role="tooltip" hidden></div>
    </section>

    <section class="topo-panel topo-3d">
      <div class="topo-stage" data-role="stage"></div>
      <div class="topo-status" data-role="status">Preparing the 3D view.</div>
    </section>
  </div>

  <div class="topo-details" data-role="details">
    Select an element in the diagram to see its residues.
  </div>

  <script type="application/json" data-role="payload">{blob}</script>
  <script>{VIEWER_JS}</script>
  <script>{TOPOLOGY_JS}</script>
  <script>window.__topoBoot && window.__topoBoot("{root_id}");</script>
</div>
"""
