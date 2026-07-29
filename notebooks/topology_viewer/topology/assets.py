"""Static CSS and JavaScript for the renderer.

Kept out of render.py so the Python there stays readable, and so the browser
code can be linted or extracted to real .js files later without unpicking
f-strings.
"""

TOPOLOGY_CSS = """
__ROOT__ {
  --ink: #16202b;
  --muted: #5f6b7a;
  --line: #d3dce6;
  --panel: #ffffff;
  --wash: #f6f8fb;
  --helix: #d4693f;
  --helix-dim: #f0c9b8;
  --strand: #3f72b0;
  --strand-dim: #bed2e8;
  --coil: #93a1b3;
  --select: #12203a;
  --font: "Inter", "Segoe UI", system-ui, sans-serif;
  --mono: "SF Mono", "Cascadia Mono", Menlo, Consolas, monospace;
  font-family: var(--font);
  color: var(--ink);
  display: block;
  width: 100%;
}
__ROOT__ * { box-sizing: border-box; }

__ROOT__ .topo-head { margin-bottom: 10px; }
__ROOT__ .topo-title { font-size: 14px; font-weight: 600; letter-spacing: -0.01em; }
__ROOT__ .topo-provenance { font-size: 12px; color: var(--muted); margin-top: 2px; }

__ROOT__ .topo-controls {
  display: flex; flex-wrap: wrap; align-items: flex-end;
  gap: 10px 14px; padding: 10px 12px; margin-bottom: 10px;
  background: var(--wash); border: 1px solid var(--line); border-radius: 8px;
}
__ROOT__ .topo-field { display: flex; flex-direction: column; gap: 3px; }
__ROOT__ .topo-field > span {
  font-size: 10px; font-weight: 600; text-transform: uppercase;
  letter-spacing: 0.06em; color: var(--muted);
}
__ROOT__ select, __ROOT__ input[type="text"] {
  font: inherit; font-size: 13px; padding: 5px 8px;
  border: 1px solid var(--line); border-radius: 6px;
  background: var(--panel); color: var(--ink);
}
__ROOT__ .topo-actions { display: flex; gap: 6px; align-items: center; margin-left: auto; }
__ROOT__ .topo-actions button {
  font: inherit; font-size: 13px; padding: 5px 11px; cursor: pointer;
  border: 1px solid var(--line); border-radius: 6px;
  background: var(--panel); color: var(--ink);
}
__ROOT__ .topo-actions button:hover { background: #eef2f7; }
__ROOT__ select:focus-visible, __ROOT__ button:focus-visible, __ROOT__ input:focus-visible {
  outline: 2px solid var(--strand); outline-offset: 1px;
}

__ROOT__ .topo-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; }
@media (max-width: 900px) { __ROOT__ .topo-grid { grid-template-columns: 1fr; } }

__ROOT__ .topo-panel {
  position: relative; background: var(--panel);
  border: 1px solid var(--line); border-radius: 8px; overflow: hidden;
}
__ROOT__ .topo-2d { min-height: 560px; }
__ROOT__ .topo-2d svg { display: block; width: 100%; height: 560px; cursor: grab; }
__ROOT__ .topo-2d svg.dragging { cursor: grabbing; }
__ROOT__ .topo-stage { width: 100%; height: 528px; background: #fff; }
__ROOT__ .topo-status {
  font-size: 11px; color: var(--muted); padding: 6px 10px;
  border-top: 1px solid var(--line); background: var(--wash);
  white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
}

__ROOT__ .sse { cursor: pointer; }
__ROOT__ .sse-shape { stroke: var(--ink); stroke-width: 1.4; transition: opacity 120ms; }
__ROOT__ .sse-label {
  font-size: 11.5px; font-weight: 700; fill: var(--ink);
  pointer-events: none; text-anchor: middle; dominant-baseline: middle;
}
__ROOT__ .sheet-badge {
  font-size: 10px; font-weight: 700; fill: #1d2733;
  pointer-events: none; text-anchor: middle; dominant-baseline: middle;
}
__ROOT__ .sse-resnum {
  font-size: 9.5px; font-weight: 500; fill: var(--muted);
  pointer-events: none; dominant-baseline: middle;
}
__ROOT__ .sse:hover .sse-shape { opacity: 0.82; }
__ROOT__ .sse.selected .sse-shape { stroke: var(--select); stroke-width: 3; }
__ROOT__ .segment-band { fill: #f4f7fa; stroke: #e4ebf2; stroke-width: 1; }
__ROOT__ .segment-label {
  font-size: 10px; font-weight: 600; fill: #93a1b3;
  letter-spacing: 0.04em; pointer-events: none;
}
__ROOT__ .connector { fill: none; stroke: #6b7887; stroke-width: 1.4; stroke-linejoin: miter; }
__ROOT__ .connector.dim { opacity: 0.35; }
__ROOT__ .connector.active { stroke: var(--select); stroke-width: 2.6; }
__ROOT__ .terminus {
  font-size: 12px; font-weight: 700; fill: var(--muted);
  text-anchor: middle; dominant-baseline: middle;
}
__ROOT__ .terminus-stub { stroke: var(--coil); stroke-width: 1.8; fill: none; }
__ROOT__ .site-mark { stroke: #ffffff; stroke-width: 1.4; pointer-events: none; }
__ROOT__ .site-count {
  font-size: 8.5px; font-weight: 700; fill: #ffffff;
  text-anchor: middle; dominant-baseline: middle; pointer-events: none;
}
__ROOT__ .site-chip {
  display: inline-block; width: 9px; height: 9px; border-radius: 50%;
  margin-right: 5px; vertical-align: -1px;
}
__ROOT__ .topo-filters {
  display: flex; flex-wrap: wrap; align-items: center; gap: 6px 16px;
  padding: 8px 12px; margin-bottom: 10px; font-size: 12.5px;
  background: #fbfcfe; border: 1px solid var(--line); border-radius: 8px;
}
__ROOT__ .topo-filters label { display: flex; align-items: center; gap: 5px; }
__ROOT__ .topo-legend {
  display: flex; flex-wrap: wrap; gap: 4px 14px;
  width: 100%; padding-top: 4px; font-size: 12px; color: var(--ink);
}
__ROOT__ .legend-item { display: flex; align-items: center; }
__ROOT__ .topo-filters .muted { color: var(--muted); }

__ROOT__ .topo-tooltip {
  position: absolute; z-index: 20; pointer-events: none;
  background: rgba(19, 28, 39, 0.95); color: #fff;
  font-size: 11.5px; line-height: 1.45; padding: 6px 9px;
  border-radius: 5px; max-width: 260px;
}
__ROOT__ .topo-tooltip[hidden] { display: none; }
__ROOT__ .topo-tooltip code { font-family: var(--mono); color: #cfe3ff; }

__ROOT__ .topo-details {
  margin-top: 10px; padding: 10px 12px; font-size: 13px; line-height: 1.55;
  background: var(--wash); border: 1px solid var(--line); border-radius: 8px;
  min-height: 44px; color: var(--ink);
}
__ROOT__ .topo-details code { font-family: var(--mono); font-size: 12px; }
__ROOT__ .topo-details .muted { color: var(--muted); }

@media (prefers-reduced-motion: reduce) {
  __ROOT__ .sse-shape { transition: none; }
}
"""


VIEWER_JS = r"""
(function () {
  if (window.__topoViewers) return;

  function loadScript(url, id) {
    return new Promise(function (resolve, reject) {
      var existing = document.getElementById(id);
      if (existing) {
        if (existing.dataset.loaded === "1") return resolve();
        existing.addEventListener("load", function () { resolve(); }, { once: true });
        existing.addEventListener("error", reject, { once: true });
        return;
      }
      var tag = document.createElement("script");
      tag.id = id;
      tag.src = url;
      tag.async = true;
      tag.addEventListener("load", function () { tag.dataset.loaded = "1"; resolve(); }, { once: true });
      tag.addEventListener("error", reject, { once: true });
      document.head.appendChild(tag);
    });
  }

  function loadCss(url, id) {
    if (document.getElementById(id)) return;
    var link = document.createElement("link");
    link.id = id;
    link.rel = "stylesheet";
    link.href = url;
    document.head.appendChild(link);
  }

  /* ----------------------------------------------------------------
     Mol*

     The original code called viewer.structureInteractivity() on every
     pointer move. That is the high-level convenience API: each call
     compiles a MolQL query and runs it over the whole structure, so
     dragging across a helix queued dozens of compilations a second and
     starved the click that followed.

     Here the structure is walked once after load to build a
     chain:seq -> Loci map, turning every later highlight into a map
     lookup. Mol*'s UMD bundle does not formally guarantee the internals
     that walk needs, so the build is wrapped and falls back to the
     public query API if the shape ever changes. Even on the fallback
     path the pointer-move calls are gone, which is the larger win.
     ---------------------------------------------------------------- */
  function MolstarAdapter(options) {
    this.options = options || {};
    this.viewer = null;
    this.plugin = null;
    this.index = null;
    this.lastHighlight = "";
    this.scope = { mode: "complex", chain: "" };
  }

  MolstarAdapter.prototype.name = "Mol*";

  MolstarAdapter.prototype.init = function (node) {
    var self = this;
    var version = this.options.molstarVersion || "4.19.0";
    var base = "https://cdn.jsdelivr.net/npm/molstar@" + version + "/build/viewer/";
    loadCss(base + "molstar.css", "topo-molstar-css");
    return loadScript(base + "molstar.js", "topo-molstar-js").then(function () {
      return window.molstar.Viewer.create(node, {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: false,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        viewportBackgroundColor: { r: 255, g: 255, b: 255 }
      });
    }).then(function (viewer) {
      self.viewer = viewer;
      self.plugin = viewer.plugin || viewer._plugin || null;
      return self;
    });
  };

  MolstarAdapter.prototype.loadStructure = function (spec) {
    var self = this;
    var format = spec.format === "pdb" ? "pdb" : "mmcif";
    var promise;
    if (spec.url) {
      promise = this.viewer.loadStructureFromUrl(spec.url, format, false);
    } else {
      promise = this.viewer.loadStructureFromData(spec.data, format, false);
    }
    return promise.then(function () {
      self.index = null;
      try { self.index = self._buildIndex(); } catch (error) { self.index = null; }
      return { indexed: !!self.index };
    });
  };

  MolstarAdapter.prototype._structure = function () {
    if (!this.plugin) return null;
    var current = this.plugin.managers.structure.hierarchy.current;
    var entry = current && current.structures && current.structures[0];
    return (entry && entry.cell && entry.cell.obj && entry.cell.obj.data) || null;
  };

  /* Walk every unit once and record, per residue, the contiguous span of
     element positions it occupies. Spans are stored as plain numbers; the
     Loci object is assembled on demand from cached bounds. */
  MolstarAdapter.prototype._buildIndex = function () {
    var structure = this._structure();
    if (!structure || !structure.units || !structure.units.length) return null;

    var molstar = window.molstar;
    var OrderedSet =
      (molstar.OrderedSet) ||
      (molstar.StructureElement && molstar.StructureElement.OrderedSet) ||
      null;
    if (!OrderedSet || !OrderedSet.ofBounds) return null;

    var index = new Map();

    for (var u = 0; u < structure.units.length; u++) {
      var unit = structure.units[u];
      var model = unit.model;
      if (!model || !model.atomicHierarchy) continue;

      var hierarchy = model.atomicHierarchy;
      var residueOf = hierarchy.residueAtomSegments && hierarchy.residueAtomSegments.index;
      var chainOf = hierarchy.chainAtomSegments && hierarchy.chainAtomSegments.index;
      if (!residueOf || !chainOf) continue;

      var authSeq = hierarchy.residues.auth_seq_id;
      var authAsym = hierarchy.chains.auth_asym_id;
      if (!authSeq || !authAsym) continue;

      var elements = unit.elements;
      var currentResidue = -1;
      var spanStart = 0;

      for (var i = 0; i <= elements.length; i++) {
        var residueIndex = i < elements.length ? residueOf[elements[i]] : -2;
        if (residueIndex === currentResidue) continue;

        if (currentResidue >= 0) {
          var atom = elements[spanStart];
          var chainIndex = chainOf[atom];
          var key = authAsym.value(chainIndex) + ":" + authSeq.value(currentResidue);
          if (!index.has(key)) {
            index.set(key, { unit: unit, start: spanStart, end: i, structure: structure });
          }
        }
        currentResidue = residueIndex;
        spanStart = i;
      }
    }

    if (!index.size) return null;
    this._OrderedSet = OrderedSet;
    return index;
  };

  MolstarAdapter.prototype._loci = function (chain, seq) {
    if (!this.index) return null;
    var record = this.index.get(chain + ":" + seq);
    if (!record) return null;
    return {
      kind: "element-loci",
      structure: record.structure,
      elements: [
        { unit: record.unit, indices: this._OrderedSet.ofBounds(record.start, record.end) }
      ]
    };
  };

  MolstarAdapter.prototype.highlight = function (chain, seq) {
    var key = chain + ":" + seq;
    if (key === this.lastHighlight) return;
    this.lastHighlight = key;

    var loci = this._loci(chain, seq);
    if (loci && this.plugin) {
      this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci: loci });
      return;
    }
    this._fallback(chain, seq, "highlight");
  };

  MolstarAdapter.prototype.select = function (chain, start, stop) {
    var seq = start;
    var loci = this._loci(chain, seq);
    if (loci && this.plugin) {
      this.plugin.managers.interactivity.lociSelects.selectOnly({ loci: loci });
      if (this.options.followCamera) {
        try { this.plugin.managers.camera.focusLoci(loci); } catch (error) { /* no-op */ }
      }
      return;
    }
    this._fallback(chain, seq, "select");
  };

  /* Public API. Correct but slow, so it only runs when the index is absent. */
  MolstarAdapter.prototype._fallback = function (chain, seq, action) {
    if (!this.viewer || !this.viewer.structureInteractivity) return;
    try {
      this.viewer.structureInteractivity({
        elements: {
          auth_asym_id: chain,
          beg_auth_seq_id: seq,
          end_auth_seq_id: seq
        },
        action: action
      });
    } catch (error) { /* no-op */ }
  };

  MolstarAdapter.prototype.clearHighlight = function () {
    this.lastHighlight = "";
    if (!this.plugin) return;
    try { this.plugin.managers.interactivity.lociHighlights.clearHighlights(); }
    catch (error) { /* no-op */ }
  };

  MolstarAdapter.prototype.setScope = function (mode, chain) {
    this.scope = { mode: mode, chain: chain };
    if (!this.plugin || !this.viewer) return;
    try {
      if (this.viewer.setSubtreeVisibility) return;
      var hidden = mode === "chain";
      this.plugin.managers.structure.component.setOptions(
        Object.assign({}, this.plugin.managers.structure.component.state.options)
      );
      void hidden;
    } catch (error) { /* no-op */ }
  };

  MolstarAdapter.prototype.dispose = function () {
    try { this.viewer && this.viewer.dispose && this.viewer.dispose(); }
    catch (error) { /* no-op */ }
    this.viewer = null; this.plugin = null; this.index = null;
  };

  /* ----------------------------------------------------------------
     NGL

     Smaller bundle, so first paint is quicker. Selection strings are
     cheap here, so no index is needed: a highlight is one representation
     whose selection is reassigned.
     ---------------------------------------------------------------- */
  function NglAdapter(options) {
    this.options = options || {};
    this.stage = null;
    this.component = null;
    this.baseRep = null;
    this.selectCartoon = null;
    this.siteRep = null;
    this.overlayReps = null;
    this.hoverRep = null;
    this.lastHighlight = "";
    this.chain = "";
  }

  NglAdapter.prototype.name = "NGL";

  /* Everything is grey so the one thing that matters can be red. A structure
     coloured by chain or by spectrum has no free visual channel left: a
     highlight has to compete with colour that is already carrying meaning.
     Neutral grey costs nothing and makes the selection unmissable. */
  NglAdapter.prototype.BASE_COLOUR = 0xc8ccd2;
  NglAdapter.prototype.SELECT_COLOUR = 	0x020066;
  NglAdapter.prototype.HOVER_COLOUR = 0xf2a71e;

  NglAdapter.prototype.init = function (node) {
    var self = this;
    var url = "https://cdn.jsdelivr.net/npm/ngl@2.3.1/dist/ngl.js";
    return loadScript(url, "topo-ngl-js").then(function () {
      self.stage = new window.NGL.Stage(node, {
        backgroundColor: "white",
        tooltip: false
      });
      self._onResize = function () { self.stage && self.stage.handleResize(); };
      window.addEventListener("resize", self._onResize);
      return self;
    });
  };

  NglAdapter.prototype.loadStructure = function (spec) {
    var self = this;
    var source = spec.url ? spec.url : new Blob([spec.data], { type: "text/plain" });
    var settings = spec.url ? {} : { ext: spec.format === "pdb" ? "pdb" : "cif" };

    return this.stage.loadFile(source, settings).then(function (component) {
      self.component = component;

      self.baseRep = component.addRepresentation("cartoon", {
        color: self.BASE_COLOUR,
        opacity: 1.0,
        smoothSheet: true
      });

      // Drawn once with an empty selection, then only the selection string is
      // reassigned. Rebuilding representations per click is what makes a linked
      // viewer feel sluggish; setSelection is a cheap update.
      self.selectCartoon = component.addRepresentation("cartoon", {
        sele: "none",
        color: self.SELECT_COLOUR,
        smoothSheet: true
      });
      // Sticks are reserved for a single clicked site. Drawing them for a
      // whole element buries the cartoon in geometry and makes the shape of
      // the selection harder to read, not easier.
      self.siteRep = component.addRepresentation("licorice", {
        sele: "none",
        color: self.SELECT_COLOUR,
        radiusScale: 1.8
      });
      self.hoverRep = component.addRepresentation("spacefill", {
        sele: "none",
        color: self.HOVER_COLOUR,
        radiusScale: 0.32
      });

      component.autoView();
      return { indexed: true };
    });
  };

  NglAdapter.prototype._range = function (chain, start, stop) {
    var span = (start === stop) ? String(start) : (start + "-" + stop);
    return chain ? (span + ":" + chain) : span;
  };

  NglAdapter.prototype.highlight = function (chain, seq) {
    var key = chain + ":" + seq;
    if (key === this.lastHighlight || !this.hoverRep) return;
    this.lastHighlight = key;
    this.hoverRep.setSelection(this._range(chain, seq, seq));
  };

  /* A click selects the whole secondary-structure element, not one residue.
     The diagram's unit of meaning is the element, so that is what lights up --
     by colour alone, which is enough to locate it. */
  NglAdapter.prototype.select = function (chain, start, stop) {
    if (!this.selectCartoon) return;
    var sele = this._range(chain, start, stop);
    this.selectCartoon.setSelection(sele);
    if (this.siteRep) this.siteRep.setSelection("none");
    if (this.options.followCamera && this.component) {
      this.component.autoView(sele, 400);
    }
  };

  /* A single site is the one case where sticks earn their place: the question
     there is which way a sidechain points, which a cartoon cannot answer. */
  NglAdapter.prototype.selectSite = function (chain, seq) {
    if (!this.siteRep) return;
    this.siteRep.setSelection(this._range(chain, seq, seq) + " and sidechainAttached");
  };

  NglAdapter.prototype.clearHighlight = function () {
    this.lastHighlight = "";
    if (this.hoverRep) this.hoverRep.setSelection("none");
  };

  NglAdapter.prototype.clearSelection = function () {
    if (this.selectCartoon) this.selectCartoon.setSelection("none");
    if (this.siteRep) this.siteRep.setSelection("none");
  };

  NglAdapter.prototype.setScope = function (mode, chain) {
    this.chain = chain;
    if (!this.baseRep) return;
    this.baseRep.setSelection(mode === "chain" && chain ? (":" + chain) : "all");
    if (this.component) this.component.autoView(400);
  };

  /* Every annotated site at once, as sticks over the grey cartoon.

     Three representations rather than one, because the categories carry
     different colours and NGL colours a representation as a whole. They are
     built once and only their selection strings are rewritten, so filtering
     stays cheap. */
  NglAdapter.prototype.SITE_COLOURS = {
    ptm: "#F0C808",
    variant: "#D7263D",
    both: "#2E9E4F"
  };

  NglAdapter.prototype.setSiteOverlay = function (chain, groups) {
    if (!this.component) return;
    var self = this;
    groups = groups || {};

    if (!this.overlayReps) {
      this.overlayReps = {};
      ["ptm", "variant", "both"].forEach(function (category) {
        self.overlayReps[category] = self.component.addRepresentation("licorice", {
          sele: "none",
          color: self.SITE_COLOURS[category],
          radiusScale: 1.6,
          /* Sidechains attached to the backbone, so a modified residue reads as
             a residue rather than a floating fragment. */
          multipleBond: "symmetric"
        });
      });
    }

    ["ptm", "variant", "both"].forEach(function (category) {
      var positions = groups[category] || [];
      var unique = positions.filter(function (value, index, all) {
        return all.indexOf(value) === index;
      });
      if (!unique.length) {
        self.overlayReps[category].setSelection("none");
        return;
      }
      var sele = "(" + unique.join(" or ") + ")" +
                 (chain ? " and :" + chain : "") + " and sidechainAttached";
      self.overlayReps[category].setSelection(sele);
    });
  };

  NglAdapter.prototype.dispose = function () {
    if (this._onResize) window.removeEventListener("resize", this._onResize);
    try { this.stage && this.stage.dispose(); } catch (error) { /* no-op */ }
    this.stage = null; this.component = null;
  };

  window.__topoViewers = {
    molstar: MolstarAdapter,
    ngl: NglAdapter
  };
})();
"""


TOPOLOGY_JS = r"""
(function () {
  if (window.__topoBoot) return;

  var SVG_NS = "http://www.w3.org/2000/svg";

  var HELIX_FILL = "#d9606b";

  /* Strands are coloured by the sheet they pair into, so the sheets are
     readable at a glance instead of having to be traced pairing by pairing.
     Hues are well separated and hold up next to the single helix colour. */
  var SHEET_FILL = [
    "#e8912d", "#5fb9e0", "#e0563f", "#6fbb6f",
    "#9186cf", "#cf9a63", "#3fa39b", "#c76fae"
  ];
  var LONE_STRAND_FILL = "#8fa2b8";

  function el(name, attributes) {
    var node = document.createElementNS(SVG_NS, name);
    for (var key in attributes) {
      if (attributes[key] !== null && attributes[key] !== undefined) {
        node.setAttribute(key, attributes[key]);
      }
    }
    return node;
  }

  function escapeHtml(value) {
    return String(value === null || value === undefined ? "" : value)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }

  /* Blue-to-orange ramp for confidence, matching AlphaFold's convention
     closely enough to be read the same way. */
  function confidenceColour(value, range) {
    if (value === null || value === undefined || !range) return "var(--coil)";
    var low = range[0], high = range[1];
    var t = high > low ? (value - low) / (high - low) : 0.5;
    t = Math.max(0, Math.min(1, t));
    var stops = [
      [0.00, [230, 116, 58]],
      [0.50, [244, 213, 96]],
      [0.75, [101, 203, 243]],
      [1.00, [0, 83, 214]]
    ];
    for (var i = 1; i < stops.length; i++) {
      if (t <= stops[i][0]) {
        var span = stops[i][0] - stops[i - 1][0];
        var local = span > 0 ? (t - stops[i - 1][0]) / span : 0;
        var a = stops[i - 1][1], b = stops[i][1];
        return "rgb(" +
          Math.round(a[0] + (b[0] - a[0]) * local) + "," +
          Math.round(a[1] + (b[1] - a[1]) * local) + "," +
          Math.round(a[2] + (b[2] - a[2]) * local) + ")";
      }
    }
    return "rgb(0,83,214)";
  }

  function roundedPolyline(points, radius) {
    if (points.length < 4) return "";
    var parts = ["M " + points[0] + " " + points[1]];
    for (var i = 2; i < points.length - 2; i += 2) {
      var px = points[i - 2], py = points[i - 1];
      var cx = points[i], cy = points[i + 1];
      var nx = points[i + 2], ny = points[i + 3];

      var inLength = Math.hypot(cx - px, cy - py);
      var outLength = Math.hypot(nx - cx, ny - cy);
      var r = Math.min(radius, inLength / 2, outLength / 2);
      if (!isFinite(r) || r <= 0.5) { parts.push("L " + cx + " " + cy); continue; }

      var t1 = r / inLength, t2 = r / outLength;
      parts.push("L " + (cx - (cx - px) * t1) + " " + (cy - (cy - py) * t1));
      parts.push("Q " + cx + " " + cy + " " +
                 (cx + (nx - cx) * t2) + " " + (cy + (ny - cy) * t2));
    }
    parts.push("L " + points[points.length - 2] + " " + points[points.length - 1]);
    return parts.join(" ");
  }

  function Topology(root) {
    this.root = root;
    this.payload = JSON.parse(root.querySelector("[data-role='payload']").textContent);
    this.svg = root.querySelector("[data-role='svg']");
    this.viewport = root.querySelector("[data-role='viewport']");
    this.tooltip = root.querySelector("[data-role='tooltip']");
    this.details = root.querySelector("[data-role='details']");
    this.status = root.querySelector("[data-role='status']");
    this.stage = root.querySelector("[data-role='stage']");

    this.mode = "sheet";
    this.colourBy = "ss";
    this.showSites = true;
    this.filters = { ptm: true, variant: true, minScore: 0 };
    this.selectedId = null;
    this.transform = { x: 0, y: 0, k: 1 };
    this.adapter = null;
    this.pendingHighlight = null;
    this.frameQueued = false;

    this.residueByKey = new Map();
    for (var i = 0; i < this.payload.residues.length; i++) {
      var residue = this.payload.residues[i];
      this.residueByKey.set(residue.chain + ":" + residue.seq, residue);
    }
    this.elementById = new Map();
    for (var j = 0; j < this.payload.elements.length; j++) {
      this.elementById.set(this.payload.elements[j].id, this.payload.elements[j]);
    }

    this.bindControls();
    this.draw();
    this.initViewer();
  }

  Topology.prototype.layout = function () {
    return this.payload.layouts[this.mode] || this.payload.layouts.sheet;
  };

  Topology.prototype.draw = function () {
    var layout = this.layout();
    var extents = layout.extents;
    var pad = 48;
    this.svg.setAttribute("viewBox",
      (extents[0] - pad) + " " + (extents[1] - pad) + " " +
      (extents[2] - extents[0] + pad * 2) + " " + (extents[3] - extents[1] + pad * 2));

    while (this.viewport.firstChild) this.viewport.removeChild(this.viewport.firstChild);

    var segmentLayer = el("g", { "data-layer": "segments" });
    this.viewport.appendChild(segmentLayer);
    (layout.segments || []).forEach(function (segment) {
      segmentLayer.appendChild(el("rect", {
        x: segment.x, y: segment.y, width: segment.width, height: segment.height,
        rx: 10, ry: 10, class: "segment-band"
      }));
      var caption = el("text", {
        x: segment.x + 12, y: segment.y + 15, class: "segment-label"
      });
      caption.textContent = segment.label;
      segmentLayer.appendChild(caption);
    });

    var connectorLayer = el("g", { "data-layer": "connectors" });
    var elementLayer = el("g", { "data-layer": "elements" });
    this.viewport.appendChild(connectorLayer);
    this.viewport.appendChild(elementLayer);

    var self = this;

    layout.connectors.forEach(function (connector) {
      var path = el("path", {
        d: roundedPolyline(connector.path, 4),
        class: "connector",
        "data-source": connector.source,
        "data-target": connector.target
      });
      connectorLayer.appendChild(path);
    });

    layout.termini.forEach(function (terminus) {
      connectorLayer.appendChild(el("path", {
        d: "M " + terminus.anchor[0] + " " + terminus.anchor[1] +
           " L " + terminus.x + " " + terminus.y,
        class: "terminus-stub"
      }));
      var label = el("text", {
        x: terminus.x,
        y: terminus.y + (terminus.type === "N" ? -12 : 12),
        class: "terminus"
      });
      label.textContent = terminus.type;
      connectorLayer.appendChild(label);
    });

    layout.elements.forEach(function (item) {
      elementLayer.appendChild(self.drawElement(item));
    });

    this.applyTransform();
  };

  /* Merge marks that would land within a few pixels of each other. */
  Topology.prototype.clusterSites = function (sites, item) {
    var ordered = sites.slice().sort(function (a, b) { return a.t - b.t; });
    var minimum = item.h > 0 ? (11 / item.h) : 1;
    var clusters = [];
    ordered.forEach(function (site) {
      var last = clusters[clusters.length - 1];
      if (last && Math.abs(site.t - last.t) < minimum) {
        last.members.push(site);
        last.t = (last.t * (last.members.length - 1) + site.t) / last.members.length;
      } else {
        clusters.push({ t: site.t, members: [site] });
      }
    });
    return clusters;
  };

  Topology.prototype.passesFilter = function (site) {
    if (site.kind === "ptm" && !this.filters.ptm) return false;
    if (site.kind === "variant" && !this.filters.variant) return false;
    if (this.filters.minScore > 0) {
      if (site.score === null || site.score === undefined) return false;
      if (site.score < this.filters.minScore) return false;
    }
    return true;
  };

  Topology.prototype.sheetFill = function (meta) {
    if (!meta || !meta.sheet_index) return LONE_STRAND_FILL;
    return SHEET_FILL[(meta.sheet_index - 1) % SHEET_FILL.length];
  };

  Topology.prototype.fillFor = function (item) {
    var meta = this.elementById.get(item.id);
    if (this.colourBy === "confidence") {
      var source = this.payload.source;
      var values = [];
      for (var seq = item.start; seq <= item.stop; seq++) {
        var residue = this.residueByKey.get(this.payload.chain + ":" + seq);
        if (residue && residue.confidence !== null && residue.confidence !== undefined) {
          values.push(residue.confidence);
        }
      }
      if (!values.length) return "var(--coil)";
      var mean = values.reduce(function (a, b) { return a + b; }, 0) / values.length;
      return confidenceColour(mean, source.confidence_range);
    }
    if (this.colourBy === "density") {
      var annotations = this.payload.annotations;
      var peak = (annotations && annotations.max_density) || 0;
      var here = (meta.sites || []).filter(this.passesFilter, this).length;
      if (!peak) return "#e7ecf2";
      var t = Math.min(1, here / peak);
      /* Pale to saturated on a single hue. Shape still separates helix from
         strand, so the fill is free to carry a quantity instead. */
      var r = Math.round(238 + (155 - 238) * t);
      var g = Math.round(242 + (34 - 242) * t);
      var b = Math.round(247 + (94 - 247) * t);
      return "rgb(" + r + "," + g + "," + b + ")";
    }
    if (item.kind === "helix") return HELIX_FILL;
    return this.sheetFill(meta);
  };

  Topology.prototype.drawElement = function (item) {
    var self = this;
    var meta = this.elementById.get(item.id) || {};
    var group = el("g", {
      class: "sse",
      "data-element-id": item.id,
      tabindex: "0",
      role: "button",
      "aria-label": item.kind + " " + item.id + ", residues " + item.start + " to " + item.stop
    });

    var fill = this.fillFor(item);

    if (item.kind === "strand") {
      var points = [];
      for (var i = 0; i < item.path.length; i += 2) {
        points.push(item.path[i] + "," + item.path[i + 1]);
      }
      group.appendChild(el("polygon", {
        points: points.join(" "), fill: fill, class: "sse-shape"
      }));

      /* Sheet badge. Colour alone cannot be relied on -- there may be more
         sheets than distinguishable hues, and some readers cannot separate
         them -- so the label carries the same information redundantly. */
      if (meta.sheet) {
        var badgeY = item.y + item.h / 2;
        group.appendChild(el("rect", {
          x: item.x - 13, y: badgeY - 9, width: 26, height: 18,
          rx: 5, ry: 5, fill: "rgba(255,255,255,0.92)",
          stroke: "rgba(0,0,0,0.18)", "stroke-width": 0.8,
          "pointer-events": "none"
        }));
        var badge = el("text", {
          x: item.x, y: badgeY, class: "sheet-badge"
        });
        badge.textContent = meta.sheet;
        group.appendChild(badge);
      }
    } else {
      /* A plain capsule. The rib hatching that was here read as texture rather
         than structure and fought with the label at small sizes. */
      group.appendChild(el("rect", {
        x: item.x - 16, y: item.y, width: 32, height: item.h,
        rx: 15, ry: 15, fill: fill, class: "sse-shape"
      }));
    }

    /* Site marks. Clustering is not cosmetic: a 20-residue strand is about
       100px wide, so eight sites would overlap into an unreadable smear.
       Nearby marks merge into one glyph carrying the count, and the detail is
       recovered on click. */
    var sites = (meta.sites || []).filter(function (site) { return self.passesFilter(site); });
    if (sites.length && this.showSites) {
      var clusters = self.clusterSites(sites, item);
      clusters.forEach(function (cluster) {
        var cy = item.direction > 0
          ? item.y + cluster.t * item.h
          : item.y + (1 - cluster.t) * item.h;
        var cx = item.x + (item.kind === "strand" ? 20 : 22);
        if (cluster.members.length === 1) {
          group.appendChild(el("circle", {
            cx: cx, cy: cy, r: 5,
            fill: cluster.members[0].colour, class: "site-mark"
          }));
        } else {
          group.appendChild(el("circle", {
            cx: cx, cy: cy, r: 8,
            fill: cluster.members[0].colour, class: "site-mark"
          }));
          var count = el("text", { x: cx, y: cy, class: "site-count" });
          count.textContent = cluster.members.length;
          group.appendChild(count);
        }
      });
    }

    // Identifier above the element, clear of the shape.
    var label = el("text", { x: item.x, y: item.y - 11, class: "sse-label" });
    label.textContent = item.id;
    group.appendChild(label);

    /* Residue numbers at the ends, ordered N first so the numbers themselves
       show which way the element runs even in a black-and-white printout. */
    var topNumber = item.direction > 0 ? item.start : item.stop;
    var bottomNumber = item.direction > 0 ? item.stop : item.start;
    var startText = el("text", {
      x: item.x - 20, y: item.y + 6, class: "sse-resnum", "text-anchor": "end"
    });
    startText.textContent = topNumber;
    var endText = el("text", {
      x: item.x - 20, y: item.y + item.h - 2, class: "sse-resnum", "text-anchor": "end"
    });
    endText.textContent = bottomNumber;
    group.appendChild(startText);
    group.appendChild(endText);

    group.addEventListener("mouseenter", function (event) { self.onHover(event, item); });
    group.addEventListener("mousemove", function (event) { self.onHover(event, item); });
    group.addEventListener("mouseleave", function () { self.onLeave(); });
    group.addEventListener("click", function (event) {
      self.select(item.id, self.residueAt(event, item));
    });
    group.addEventListener("keydown", function (event) {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        self.select(item.id, item.start);
      }
    });

    return group;
  };

  /* Map a pointer position onto a residue by how far along the element it
     falls, respecting the element's N-to-C direction. */
  Topology.prototype.residueAt = function (event, item) {
    var point = this.svg.createSVGPoint();
    point.x = event.clientX;
    point.y = event.clientY;
    var local;
    try { local = point.matrixTransform(this.viewport.getScreenCTM().inverse()); }
    catch (error) { return item.start; }

    var t = (local.y - item.y) / (item.h || 1);
    if (item.direction < 0) t = 1 - t;
    t = Math.max(0, Math.min(1, t));
    var span = item.stop - item.start;
    return item.start + Math.round(t * span);
  };

  /* Pointer moves update the tooltip only. The 3D highlight is deferred to
     the next animation frame and skipped when the residue has not changed,
     so a drag across a long helix issues one 3D call, not one per pixel. */
  Topology.prototype.onHover = function (event, item) {
    var seq = this.residueAt(event, item);
    var residue = this.residueByKey.get(this.payload.chain + ":" + seq);

    var box = this.root.getBoundingClientRect();
    this.tooltip.style.left = (event.clientX - box.left + 14) + "px";
    this.tooltip.style.top = (event.clientY - box.top + 14) + "px";
    this.tooltip.hidden = false;

    var confidence = residue && residue.confidence !== null && residue.confidence !== undefined
      ? Number(residue.confidence).toFixed(1) : "n/a";
    this.tooltip.innerHTML =
      "<strong>" + escapeHtml(item.id) + "</strong> " + escapeHtml(item.kind) +
      "<br>Residue <code>" + escapeHtml(seq) + "</code> " +
      escapeHtml(residue ? residue.comp : "?") +
      "<br>" + escapeHtml(this.payload.source.confidence_label) + " <code>" +
      confidence + "</code>";

    this.queueHighlight(seq);
  };

  Topology.prototype.queueHighlight = function (seq) {
    this.pendingHighlight = seq;
    if (this.frameQueued) return;
    this.frameQueued = true;
    var self = this;
    window.requestAnimationFrame(function () {
      self.frameQueued = false;
      if (self.adapter && self.pendingHighlight !== null) {
        self.adapter.highlight(self.payload.chain, self.pendingHighlight);
      }
    });
  };

  Topology.prototype.onLeave = function () {
    this.tooltip.hidden = true;
    this.pendingHighlight = null;
    if (this.adapter) this.adapter.clearHighlight();
  };

  Topology.prototype.select = function (elementId, seq) {
    var self = this;
    this.selectedId = elementId;
    var item = this.elementById.get(elementId);

    this.root.querySelectorAll("[data-element-id]").forEach(function (node) {
      node.classList.toggle("selected", node.getAttribute("data-element-id") === elementId);
    });
    this.root.querySelectorAll(".connector").forEach(function (node) {
      var linked = node.getAttribute("data-source") === elementId ||
                   node.getAttribute("data-target") === elementId;
      node.classList.toggle("active", linked);
    });

    if (!item) return;
    var residue = this.residueByKey.get(this.payload.chain + ":" + seq);
    var contacts = (this.payload.contacts || []).filter(function (contact) {
      return contact.source === elementId || contact.target === elementId;
    });

    var partners = contacts.map(function (contact) {
      var other = contact.source === elementId ? contact.target : contact.source;
      return escapeHtml(other) + " (" + escapeHtml(contact.orientation) + ")";
    }).join(", ");

    this.details.innerHTML =
      "<strong>" + escapeHtml(item.id) + "</strong> &middot; " + escapeHtml(item.ss_name) +
      " &middot; residues <code>" + escapeHtml(item.start) + "&ndash;" +
      escapeHtml(item.stop) + "</code> &middot; length <code>" +
      escapeHtml(item.length) + "</code>" +
      (residue ? " &middot; clicked <code>" + escapeHtml(residue.comp) + " " +
        escapeHtml(seq) + "</code>" : "") +
      "<br><span class='muted'>Sequence</span> <code>" + escapeHtml(item.sequence) + "</code>" +
      (partners ? "<br><span class='muted'>Pairs with</span> " + partners : "");

    var siteList = (item.sites || []).filter(function (site) { return self.passesFilter(site); });
    if (siteList.length) {
      var rows = siteList.map(function (site) {
        var where = site.uniprot_position && site.uniprot_position !== site.position
          ? escapeHtml(site.position) + " (UniProt " + escapeHtml(site.uniprot_position) + ")"
          : escapeHtml(site.position);
        return "<span class='site-chip' style='background:" + escapeHtml(site.colour) + "'></span>" +
          "<code>" + where + "</code> " + escapeHtml(site.name) +
          (site.category_label ? " <span class='muted'>" +
            escapeHtml(site.category_label) + "</span>" : "") +
          (site.score !== null && site.score !== undefined
            ? " <span class='muted'>" + Math.round(site.score * 100) + "%</span>" : "") +
          (site.detail ? " <span class='muted'>" + escapeHtml(site.detail) + "</span>" : "");
      });
      this.details.innerHTML += "<br><span class='muted'>Sites</span><br>" + rows.join("<br>");
    }

    if (this.adapter) {
      this.adapter.select(this.payload.chain, item.start, item.stop, seq);
      if (seq && this.adapter.selectSite) this.adapter.selectSite(this.payload.chain, seq);
    }
  };

  Topology.prototype.applyTransform = function () {
    this.viewport.setAttribute("transform",
      "translate(" + this.transform.x + " " + this.transform.y + ") scale(" + this.transform.k + ")");
  };

  Topology.prototype.bindControls = function () {
    var self = this;
    var root = this.root;

    root.querySelector("[data-role='layout']").addEventListener("change", function (event) {
      self.mode = event.target.value;
      self.draw();
      if (self.selectedId) self.select(self.selectedId, null);
    });

    root.querySelector("[data-role='colour']").addEventListener("change", function (event) {
      self.colourBy = event.target.value;
      self.draw();
    });

    ["ptm", "variant"].forEach(function (kind) {
      var box = root.querySelector("[data-role='filter-" + kind + "']");
      if (!box) return;
      box.addEventListener("change", function () {
        self.filters[kind] = box.checked;
        self.draw();
        self.pushSiteOverlay();
      });
    });
    var score = root.querySelector("[data-role='filter-score']");
    if (score) {
      score.addEventListener("input", function () {
        self.filters.minScore = parseFloat(score.value) || 0;
        var readout = root.querySelector("[data-role='score-readout']");
        if (readout) readout.textContent = Math.round(self.filters.minScore * 100) + "%";
        self.draw();
        self.pushSiteOverlay();
      });
    }

    var scope = root.querySelector("[data-role='scope']");
    if (scope) {
      scope.addEventListener("change", function (event) {
        if (self.adapter) self.adapter.setScope(event.target.value, self.payload.chain);
      });
    }

    root.querySelector("[data-role='engine']").addEventListener("change", function (event) {
      self.switchEngine(event.target.value);
    });

    root.querySelector("[data-action='zoom-in']").addEventListener("click", function () {
      self.transform.k = Math.min(6, self.transform.k * 1.25); self.applyTransform();
    });
    root.querySelector("[data-action='zoom-out']").addEventListener("click", function () {
      self.transform.k = Math.max(0.2, self.transform.k / 1.25); self.applyTransform();
    });
    root.querySelector("[data-action='reset']").addEventListener("click", function () {
      self.transform = { x: 0, y: 0, k: 1 }; self.applyTransform();
    });
    root.querySelector("[data-action='download']").addEventListener("click", function () {
      var blob = new Blob([JSON.stringify(self.payload, null, 2)], { type: "application/json" });
      var link = document.createElement("a");
      link.href = URL.createObjectURL(blob);
      link.download = (self.payload.source.entry_id || "topology") + "_topology.json";
      link.click();
      URL.revokeObjectURL(link.href);
    });

    var jump = root.querySelector("[data-role='jump']");
    function doJump() {
      var seq = parseInt(jump.value, 10);
      if (!isFinite(seq)) return;
      var match = null;
      self.elementById.forEach(function (item) {
        if (!match && seq >= item.start && seq <= item.stop) match = item;
      });
      if (match) {
        self.select(match.id, seq);
      } else {
        self.details.innerHTML = "<span class='muted'>Residue " + escapeHtml(seq) +
          " is not inside a helix or strand.</span>";
      }
    }
    root.querySelector("[data-action='jump']").addEventListener("click", doJump);
    jump.addEventListener("keydown", function (event) {
      if (event.key === "Enter") doJump();
    });

    this.svg.addEventListener("wheel", function (event) {
      event.preventDefault();
      var factor = event.deltaY < 0 ? 1.1 : 1 / 1.1;
      self.transform.k = Math.max(0.2, Math.min(6, self.transform.k * factor));
      self.applyTransform();
    }, { passive: false });

    var dragging = false, originX = 0, originY = 0;
    this.svg.addEventListener("mousedown", function (event) {
      if (event.target.closest(".sse")) return;
      dragging = true; originX = event.clientX; originY = event.clientY;
      self.svg.classList.add("dragging");
    });
    window.addEventListener("mousemove", function (event) {
      if (!dragging) return;
      self.transform.x += event.clientX - originX;
      self.transform.y += event.clientY - originY;
      originX = event.clientX; originY = event.clientY;
      self.applyTransform();
    });
    window.addEventListener("mouseup", function () {
      dragging = false; self.svg.classList.remove("dragging");
    });
  };

  /* All annotated sites shown at once in 3D, so the spatial clustering of
     modifications is visible without clicking through them one at a time. */
  Topology.prototype.pushSiteOverlay = function () {
    if (!this.adapter || !this.adapter.setSiteOverlay) return;
    var self = this;
    var groups = { ptm: [], variant: [], both: [] };

    function collect(site) {
      if (!self.passesFilter(site)) return;
      var category = site.category || site.kind;
      if (!groups[category]) groups[category] = [];
      groups[category].push(site.position);
    }

    (this.payload.elements || []).forEach(function (element) {
      (element.sites || []).forEach(collect);
    });
    // Sites in loops are still real sites, so they belong in the 3D view even
    // though the diagram has no element to hang them on.
    var annotations = this.payload.annotations;
    if (annotations && annotations.coil_sites) {
      annotations.coil_sites.forEach(collect);
    }

    this.adapter.setSiteOverlay(this.payload.chain, groups);
  };

  Topology.prototype.setStatus = function (message) {
    if (this.status) this.status.textContent = message;
  };

  Topology.prototype.initViewer = function () {
    var engine = this.root.querySelector("[data-role='engine']").value;
    this.switchEngine(engine);
  };

  Topology.prototype.switchEngine = function (engine) {
    var self = this;
    var source = this.payload.structure_source || {};
    if (!source.url && !source.data) {
      this.setStatus("No 3D coordinates are attached to this view.");
      return;
    }

    if (this.adapter) { this.adapter.dispose(); this.adapter = null; }
    this.stage.innerHTML = "";

    var Adapter = window.__topoViewers[engine];
    if (!Adapter) { this.setStatus("Unknown 3D engine: " + engine); return; }

    var adapter = new Adapter({ followCamera: false });
    this.setStatus("Loading " + adapter.name + ".");

    adapter.init(this.stage).then(function () {
      self.setStatus("Loading structure.");
      return adapter.loadStructure(source);
    }).then(function (result) {
      self.adapter = adapter;
      var scopeControl = self.root.querySelector("[data-role='scope']");
      if (scopeControl) adapter.setScope(scopeControl.value, self.payload.chain);
      self.pushSiteOverlay();
      var note = result && result.indexed === false
        ? " Using the fallback selection path."
        : "";
      self.setStatus(adapter.name + " ready. Click the diagram to link 2D and 3D." + note);
    }).catch(function (error) {
      self.setStatus(adapter.name + " could not load: " + (error && error.message ? error.message : error));
    });
  };

  window.__topoBoot = function (rootId) {
    var root = document.getElementById(rootId);
    if (!root || root.dataset.booted === "1") return;
    root.dataset.booted = "1";
    try { new Topology(root); }
    catch (error) {
      var details = root.querySelector("[data-role='details']");
      if (details) details.textContent = "Could not draw the topology: " + error.message;
    }
  };
})();
"""
