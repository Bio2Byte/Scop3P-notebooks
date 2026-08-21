# Topology Viewer

## Source

- Package: [`notebooks/topology_viewer/topology/`](../../notebooks/topology_viewer/topology/)
- Launcher notebook: [`notebooks/topology_viewer/topology_viewer.ipynb`](../../notebooks/topology_viewer/topology_viewer.ipynb)
- Shiny app: [`apps/topology_viewer/app.py`](../../apps/topology_viewer/app.py)

Unlike the other protocols, the science here was already a package with the notebook
reduced to a three-line launcher. The migration therefore replaces only the
`ipywidgets` control layer in `topology/app.py::make_app`; `build_view()` and every
line of browser code are reused untouched. `make_app` stays in place, because the
notebook and `notebooks/topology_viewer/test_topology.py` still drive it.

## Functional scope

- Draw a 2D secondary-structure topology diagram beside a 3D structure viewer.
- Two sources: an AlphaFold DB model or PDBe entry by UniProtKB accession, or an
  uploaded `.pdb` / `.ent` / `.cif` / `.mmcif` file.
- Format is decided by sniffing the content, not the extension, because several
  prediction tools emit PDB-format text under a `.cif` name.
- Secondary structure comes from the file when present, and is otherwise derived from
  CA coordinates. Three interchangeable backends: `mkdssp` if on PATH, then `biotite`,
  then a built-in P-SEA implementation with no dependencies at all. The provenance is
  always shown, because P-SEA and DSSP disagree about element boundaries by a residue
  or two.
- Three layouts: sheet topology, sequence order, spatial arrangement.
- Colour by structure type, by pLDDT/B-factor, or by annotation density.
- Chain selection, with the chain that actually maps to the accession preferred over
  merely the largest one.
- In accession mode only: overlay Scop3P PTMs (optionally union'd with UniProt PTM
  features) and UniProt disease variants, positioned through the SIFTS UniProt↔author
  numbering map.
- 3D engines: NGL 2.3.1 or Mol\* 4.18.0, behind one adapter interface.

## UI mapping to Shiny

| ipywidgets (`topology/app.py`) | Shiny |
|---|---|
| `accession_input` `Text` | `ui.input_text("accession", ACCESSION_LABEL)` in a `scop3p_field_row`, with a Load example button |
| `fetch_button` `Button` | `ui.input_action_button("fetch_btn")` |
| `uploader` `FileUpload` | `ui.input_file("structure_upload")` |
| `structure_picker` / `chain_picker` `Dropdown` | `ui.input_select` + `ui.update_select` |
| `ptm_button` / `variant_button` / `clear_button` | action buttons |
| `uniprot_ptm_toggle` `Checkbox` | `ui.input_checkbox("include_uniprot_ptms")` |
| `status` `HTML` + `say()` | `ui.output_text_verbatim("status")` + a `reactive.value` |
| `output` `Output` + `redraw()` | `@render.ui def topology_view()` |
| `annotation_bar.layout.display` | `ui.panel_conditional` on the mode radio |
| `guarded()` | `try/except` + `LOGGER.exception` + a status string |
| `state` dict | per-session `reactive.value`s inside `server()` |
| `_set_options()` / `suspend_redraw` | **removed** — see below |
| — (new) | `ui.input_radio_buttons("mode")` making the two sources explicit |

## Behavior preserved

- The two modes stay strictly separate. File mode makes no network requests and
  offers no PTM or variant overlay, because a locally predicted model has no reliable
  UniProt numbering. The original guessed an accession from the uploaded filename and
  drew a prediction's topology beside AlphaFold's *different* coordinates in 3D, with
  nothing to signal the mismatch.
- The annotation gate is unchanged: marks are drawn only when sites exist **and** the
  loaded structure came from an accession.
- The numbering ladder is unchanged: SIFTS columns in the updated mmCIF first (free,
  and cannot disagree with the coordinates being drawn), then the PDBe mapping API,
  and if neither is available the sites are **hidden** with an explanatory note rather
  than drawn at UniProt positions on the wrong residues.
- Residue identity is keyed `chain:seq`, so a homodimer numbered 1–250 twice
  highlights the copy the diagram is describing.
- Secondary-structure backend precedence and the provenance strings.

## Differences from Voilà

- `_set_options()` and the `suspend_redraw` flag are **gone**. They existed because
  assigning `.options` then `.value` re-entered ipywidgets' observer and redrew from
  half-updated state, and because `unobserve_all()` — the obvious workaround —
  removes ipywidgets' own internal options observer, so the next `.value` assignment
  raised `TraitError` and the widget callback swallowed it, freezing the app on
  "Fetching…". Neither failure mode exists in Shiny: `ui.update_select` carries
  choices and selection in one message with no re-validation, and everything a
  handler writes lands in a single reactive flush, so the view re-renders exactly
  once. A defensive `input.chain() in structure.residues_by_chain` fallback is kept.
- The chain selector is always visible, showing e.g. `A (352 residues)` even for a
  single-chain structure — information rather than absence.
- Session state is per-connection.
- The package is imported through `apps/common/topology_bridge.py` rather than by
  `sys.path` luck. See "Where the topology package lives" in `apps/README.md`.

## Fixed during migration

- **Mol\* never loaded.** The adapter pinned `molstar@4.19.0`, which does not exist on
  npm — the pin was written in a build environment with no CDN access and had never
  been resolved, so selecting the Mol\* engine always failed with "Mol\* could not
  load". Now pinned to `4.18.0`, the nearest published release in the same line.
  Its UMD bundle exports only `Viewer` and `ViewerAutoPreset`, so
  `MolstarAdapter._buildIndex`'s `OrderedSet` probes all miss and the adapter takes
  its documented fallback selection path; the in-view status line says so. This was
  the "not yet verified" item in the package README, and it is now verified: both
  engines load and render against the live CDNs.

## Shared UI conventions

This app uses the toolkit-wide vocabulary from `apps/common/ui_shell.py`: the accession
field is labelled **UniProtKB accession** (`ACCESSION_LABEL`), it sits in a
`scop3p_field_row` so its buttons share the input's baseline, and it carries a
**Load example** button wired to this app's `EXAMPLE_ACCESSION`. Result cards stretch to
the height of their controls card. See "Shared UI Vocabulary" in
[`apps/README.md`](../../apps/README.md).

## Validation targets

- File mode with `fixtures/annotated.pdb` reports `SS from file (HELIX/SHEET)`; with
  `fixtures/bare.pdb` it reports `built-in P-SEA` and still finds the helix that
  `fixtures/truth.txt` specifies.
- Accession mode on `P07949` offers the AlphaFold model plus the PDBe entries;
  selecting a PDB entry reports its numbering source.
- Fetching PTMs reports how many sites fall in loops and how many are absent from the
  chosen structure — the latter proves the numbering map is actually being applied.
- Loading an accession with annotations and then uploading a file must clear the marks.
- Unknown accessions and non-coordinate uploads produce a message, never a traceback.
- The 3D panel reaches "ready" on both engines.
