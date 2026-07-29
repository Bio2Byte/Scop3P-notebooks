# Protein topology viewer — phase 0

Restructures the original 4,837-line `dssp_topology_app.py` into a package, and
replaces the DSSP-file assumption with a format-agnostic structure loader.

```
topology/
  io.py        format sniffing, PDB + mmCIF parsing, SS extraction
  ss.py        secondary structure from coordinates (3 backends)
  elements.py  helix/strand elements, beta-sheet contacts
  layout.py    sheet + serpentine layouts, shared schema
  render.py    the single renderer
  assets.py    CSS and browser JS
  app.py       ipywidgets layer, accession + file modes
dssp_topology_app.py    shim so existing notebooks keep working
```

Run it — from *inside* this folder, so `topology/` is importable:

```bash
voila topology_viewer.ipynb        # or open it in JupyterLab
python3 test_topology.py           # 84 assertions, all passing
python3 preview.py fixtures/big.pdb
```

## What changed

**Uploads accept structures, not DSSP.** The old uploader listed
`.cif,.mmcif,.dssp` and routed anything that was not `.cif` to the DSSP parser,
so `.pdb` could not be selected at all and would have been misparsed if it had
been. `.pdb`, `.ent`, `.cif` and `.mmcif` all load now, and the format is
decided by sniffing content rather than trusting the extension — several
prediction tools emit PDB-format text under a `.cif` name.

**Predicted structures work.** `topology_from_alphafold_cif` read secondary
structure only from `_struct_conf` and raised if it was missing, which is
exactly the case for raw ColabFold / Boltz / Chai output. Secondary structure is
now derived from CA coordinates when the file does not carry it, via `mkdssp`
if the binary is on PATH, else biotite, else a built-in P-SEA implementation
with no dependencies at all. The provenance is shown in the UI, because P-SEA
and DSSP disagree about element boundaries by a residue or two.

**Three renderers became one.** `_legacy_topology_html`,
`_clean_alphafold_topology_html` and `_pdbe_plugin_topology_html` were ~2,500
lines of near-duplicate JavaScript in f-strings. The layout engines now emit a
shared coordinate schema and one renderer consumes it.

**Layout defects fixed.** Helix direction came from the last character of the
element id (`element["id"].endswith(("1","3","5","7","9"))`), which is arbitrary
and made connectors flip and cross; it now follows N-to-C continuity. Helices
were positioned beside whichever strand was nearest *in sequence*, which has no
spatial meaning; they are now placed by projecting their 3D centroid onto a
frame built from the sheet itself. Connectors all collapsed onto a shared
midline at `mid_y` and drew on top of each other; they now get separate lanes by
interval colouring and route over the nearer edge.

**Accession and file modes are separate.** `set_loaded` used to guess an
accession from the filename and populate the AFDB field, so uploading
`P07949_relaxed_rank_1.pdb` drew your predicted model's topology while showing
AlphaFold's *different* coordinates in 3D beside it, unlabelled. Uploads now
feed their own bytes to the viewer and never touch the network.

**Chain handling.** The old code silently picked the largest chain. There is now
a chain dropdown (hidden when there is only one) plus a 3D scope control for
chain-only vs complex-with-emphasis. Residue lookup is keyed on `chain:seq`, not
`seq`, so a homodimer numbered 1–250 twice highlights the right copy.

## Fixed after the first run

**The app froze on "Fetching…" / "Reading…".** `adopt()` called
`chain_picker.unobserve_all()` before reassigning `options` and `value`.
`unobserve_all()` removes ipywidgets' *own internal* observer — the one that
rebuilds `_options_values` when `options` change — so the next assignment to
`.value` raised `TraitError: Invalid selection: value not found`. Widget
callbacks swallow exceptions, so no traceback appeared anywhere and the status
simply stopped updating. A suppression flag replaces it, and every handler now
runs through `guarded()`, which prints the traceback into the output area rather
than discarding it.

**The diagram never drew in JupyterLab.** JupyterLab does not execute `<script>`
tags inside `display(HTML(...))` output, so the control bar rendered and the
diagram stayed blank. The view is now wrapped in an iframe via `srcdoc`, which
runs its own scripts and behaves the same under JupyterLab, Voilà and nbconvert.
`render(payload, embed="inline")` still returns the bare block, and
`standalone_document(payload)` gives a full saveable page. Adjust the frame
height with `make_app(viewer_height=1400)`.

## Test coverage

`test_topology.py` — 84 assertions covering format sniffing against misleading
extensions, all four input shapes (PDB/mmCIF × annotated/bare), SS assignment
against known-geometry ground truth, element and contact recovery, layout
invariants (everything placed, bounded extents, antiparallel strands facing
opposite ways, overlapping connectors on distinct lanes), multi-chain
separation, payload serialisation, iframe embedding, and the widget callbacks
themselves — including a regression test that drives a real `FileUpload` and
asserts the status advances past "Reading".

The built-in P-SEA agrees 90% per-residue with ground truth on the synthetic
fixture and recovers all elements with correct ordering.

## Not yet verified

**The 3D viewer integration has not been executed.** The build sandbox had no
network access to CDNs, so the Mol\* and NGL adapters are written but untested
against a live library. Two specifics to check:

- The Mol\* residue→Loci index in `assets.py::MolstarAdapter._buildIndex` reads
  `unit.model.atomicHierarchy` and needs `OrderedSet.ofBounds` from the UMD
  bundle. That export is not formally guaranteed. The code detects its absence
  and falls back to the public `structureInteractivity` API — correct but slow.
  The status line reports which path is active, so check it says "ready" without
  "Using the fallback selection path".
- Mol\* is pinned to `4.19.0` in `MolstarAdapter.init`. Adjust as needed.

Independent of the index, the pointer-move path to 3D is gone: hover updates the
tooltip only, and 3D highlights are deferred to `requestAnimationFrame` and
skipped when the residue has not changed. That was the dominant cause of the
click delay — `showResidueTooltip` was wired to `mousemove` and called
`structureInteractivity` per pixel, compiling a MolQL query each time, so clicks
queued behind a saturated main thread. Camera focus is off by default
(`followCamera: false`), since its ~500 ms transition reads as lag.

## Not built yet

The **UniProt → PDB entry dropdown is phase 2**, not in this drop. Only the
chain dropdown exists so far, populated from whatever structure is loaded. The
port target is `fetch_uniprot_pdb_xrefs` and `_chain_options_for_uniprot_pdb`
from the Scop3P notebook, plus PDBe `best_structures` for ranking and the SIFTS
numbering map.

| Phase | Work |
|---|---|
| 1 | Verify the Mol\* index, measure click latency, tune the NGL adapter |
| 2 | PDB entry dropdown from UniProt xrefs + SIFTS numbering map |
| 3 | PTMs and variants: marks on SSEs, density channel, filters, 3D overlay |
| 4 | Layout polish informed by real structures |

The payload reserves `annotations: None` for phase 3. When populated, the
renderer switches into accession mode and shows the site marks, density colour
channel and filter row; while it is `None` those controls stay hidden.
