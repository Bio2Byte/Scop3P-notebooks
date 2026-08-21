# Protocol status

Where every protocol in this repository sits along the delivery pipeline, and what is
still missing for the ones that have not made it all the way through.

Looking for how to *run* these? See [`QUICKSTART.md`](../QUICKSTART.md). This document
is about what is and is not migrated.

## The pipeline

Protocols move through four stages. Each stage has a different audience, and a protocol is
only "done" once it reaches the stage its audience needs.

| Stage | Artifact | Audience | Entry point |
|---|---|---|---|
| 1. Notebook | `notebooks/*.ipynb` | Bioinformaticians who want to read and modify the code | JupyterLab, via Binder |
| 2. Voilà app | the same notebook, widget-driven | Bench scientists who want the tool, not the code | Binder `?urlpath=voila/render/...` |
| 3. Shiny app | `apps/<name>/app.py` + logic in `apps/common/` | Anyone with a browser; multi-user, no Jupyter kernel | `shiny run`, or the portal navbar |
| 4. Container | an image from `docker/apps/<app>.Dockerfile` | Deployment — the `scop3p-toolkit` portal image is the Galaxy interactive-tool entrypoint | Docker Hub `bio2byte/scop3p-toolkit` |

Stage 3 is not a cosmetic re-skin. It replaces `ipywidgets` with Shiny reactives, moves the
science out of notebook cells into importable, unit-tested modules under `apps/common/`, and
gives every protocol the shared Scop3P visual shell (`apps/common/ui_shell.py`) plus
structured logging and FAIR session metadata (`apps/common/logging_utils.py`,
`apps/common/session_metadata.py`).

## Status matrix

| Protocol | Notebook | Voilà | Shiny app | Docker target | Portal | Use-case spec | Tests |
|---|---|---|---|---|---|---|---|
| Peptide mapper — Scop3P peptides | [`Peptide_mapper_scop3p_voila.ipynb`](../notebooks/Peptide_mapper_scop3p_voila.ipynb) | yes | [`apps/peptide_mapper`](../apps/peptide_mapper/app.py) | `peptide-mapper` (:8001) | yes | [spec](use-cases/peptide_mapper.md) | unit + smoke |
| Peptide mapper — upload your own | [`Peptide_mapper_fileupload_voila.ipynb`](../notebooks/Peptide_mapper_fileupload_voila.ipynb) | yes | same app, *Upload your own* tab | (same image) | yes | [spec](use-cases/peptide_mapper.md) | unit + smoke |
| PTM structure visualisation | [`Scop3P_PTM_structure_viz_voila_app.ipynb`](../notebooks/Scop3P_PTM_structure_viz_voila_app.ipynb) | yes | [`apps/structure_viz`](../apps/structure_viz/app.py) | `structure-viz` (:8002) | yes | [spec](use-cases/structure_viz.md) | unit + smoke |
| Biophysical prediction & mutation effect | [`Scop3P_b2b_mutation_effect_voila_app.ipynb`](../notebooks/Scop3P_b2b_mutation_effect_voila_app.ipynb) | yes | [`apps/mutation_effect`](../apps/mutation_effect/app.py) | `mutation-effect` (:8003) | yes | [spec](use-cases/mutation_effect.md) | unit + smoke |
| Protein topology viewer | [`notebooks/topology_viewer/`](../notebooks/topology_viewer/) | yes | [`apps/topology_viewer`](../apps/topology_viewer/app.py) | `topology-viewer` (:8004) | yes | [spec](use-cases/topology_viewer.md) | unit + smoke + its own script |
| RINAlign — network alignment & comparison | [`RINAlign_align_and compare_networks.ipynb`](../notebooks/RINAlign_align_and%20compare_networks.ipynb) | yes | [`apps/rinalign`](../apps/rinalign/app.py) | `rinalign` (:8005) | yes | [spec](use-cases/rinalign.md) | unit + smoke |
| Scop3P → py3Dmol peptide mapper | [`scop3p_py3dmol_mapper.ipynb`](../notebooks/scop3p_py3dmol_mapper.ipynb) | yes | not planned | — | — | — | — |
| Modifications endpoint walkthrough | `Scop3P_API.ipynb` — **file no longer in the repo** | broken badge | not planned | — | — | — | — |

Every protocol that has a Shiny app is reachable from the published all-in-one
`scop3p-toolkit` image; the per-app images are build-and-run-locally conveniences.

The portal also serves a **Help** page (`/?app=help`) describing the mission, scope and
use cases of each protocol. It is documentation rather than a protocol, so it has no
notebook, no Voila entry and no image of its own; a smoke test keeps it in step with the
navbar.

Ports are the `docker-compose.yml` host ports. The portal image (`scop3p-toolkit`, :8000)
contains every Shiny app and is the only image published by CI; the per-app images are
build-and-run-locally conveniences.

## Scop3P moved to a v1 API

Scop3P replaced its query-string endpoints with a v1 REST API: the accession moved from
a query parameter into the path, responses became bare JSON lists, and every field was
renamed to snake_case.

| | Retired | Current |
|---|---|---|
| Modifications | `/scop3p/api/modifications?accession=X` | `/scop3p/api/v1/proteins/X/modifications` |
| Peptides | `/scop3p/api/get-peptides-modifications?accession=X` | `/scop3p/api/v1/proteins/X/peptides` |

Spec: <https://iomics.ugent.be/scop3p/api/v1/openapi.json>

**Everything in this repository now targets v1** — the five Shiny apps via
`apps/common/services.py::Scop3PClient`, the five affected notebooks via a
`_scop3p_get` helper of their own, and `notebooks/topology_viewer/`, which already
tried v1 first.

The `scop3p` PyPI client is deliberately unused: 1.1.0, the latest release, still
targets all four retired endpoints.

### Logging: one scheme across five protocols

Every protocol writes a **step-by-step experiment record** -- protocol opened, values
entered, options selected, actions clicked, and what each produced -- to its own file
alongside the combined log. Steps are numbered and carry a session id, so the sequence is
explicit and two concurrent browser sessions stay separable.

Levels carry meaning: `DEBUG` for captured third-party output, `INFO` for the record,
`WARNING` for a blocked-but-continuing action, `ERROR` for a failure (traceback on the app
logger, one ordered step in the record). A failure is never recorded at `INFO`, because
"0 results" and "the lookup failed" must not read the same.

`b2bTools` reports progress with `print` across 215 call sites and its dependencies warn
once per prediction, none of which reached the log file at all (`warnings` go to stderr).
Both are now captured to `DEBUG`: console output over three predictions went from 20 lines
to 1, with nothing discarded. There are no `print` statements left in `apps/`, and tests
walk each module's AST to keep it that way, to keep every b2bTools call inside the capture
wrapper, and to keep each app's trail per-session rather than per-process.

Interactive lookups (PTMs, variants, sequence, cross-references, SIFTS) share one bounded
HTTP policy: a 5s connect timeout, a 20s read timeout and one retry. Both bounds matter and
both were observed in practice -- `rest.uniprot.org` timing out on connect while a direct
request moments earlier answered in 0.04s, and the EBI Proteins API accepting a connection
in 0.06s and then never answering. These calls run in synchronous reactive effects, so an
unbounded read freezes every connected session, not just the one that made it. File
downloads are deliberately excluded and keep their longer timeout.

### Why this failure was hard to spot

Scop3P serves its single-page app from a `GET /scop3p/{catchall}` route, so a request
to a retired endpoint returns **200 OK with `content-type: text/html`**, not a 404.
`raise_for_status()` passes and the only symptom is
`JSONDecodeError: Expecting value: line 1 column 1` from deep inside `requests` — or,
where the caller swallowed exceptions, a silent report of zero modifications. Every
fetch path in the repo now checks the content type and says so instead.

## Per-protocol notes

### Peptide mapper (Scop3P) — stage 4

Maps Scop3P phospho-peptides onto an AlphaFold model. Colour grammar: grey cartoon = whole
protein, blue = union of selected peptide spans, red = their intersection, magenta =
modified sites. Filter DSL supports substring, `70-90`, `>=150`, `<=300`, and a bare `154`.

- APIs: Scop3P `get-peptides-modifications`; AlphaFold DB (`v6` with a `v4` fallback).
- Browser libraries: NGL from unpkg, injected by `apps/common/viewer.py`.
- Both gaps recorded here previously are closed: the app now has an *Upload your own*
  source tab, and session state is per-connection (it was a module-level
  `PeptideMapperController` shared by every browser, which showed one user's peptide
  table to another). The server-side `exports/` write is kept for parity alongside the
  new browser downloads.

### Peptide mapper (upload) — stage 4

Same mapping and colour grammar, but the peptide table comes from the user's own MS search
instead of Scop3P. Columns are auto-detected by keyword across five dropdowns
(protein ID / peptide sequence / start / end / UniProt position); exports include the raw
PDB, a styled HTML session, and a TSV of mapped residues.

- APIs: AlphaFold DB only. No Scop3P call.
- Migrated as a second *data source* inside `apps/peptide_mapper` rather than as a
  separate app: the tabs cover only how the peptide table is obtained, and an uploaded
  table is normalised into the same column schema the Scop3P path produces, so
  `filter_peptides`, `build_options` and `map_selection` are reused unchanged.
- Improvements over the notebook: protein identifiers such as `sp|P07949|RET_HUMAN` and
  `P07949;Q12345` are reduced to the bare accession AlphaFold needs; an undetectable
  column is named instead of silently defaulting to the first one; and the exports are
  browser downloads rather than writes to a shared server-side `exports/` path.

### PTM structure visualisation — stage 4

The flagship six-tab protocol: Scop3P + UniProt PTMs, UniProt disease variants, 3D mapping
onto experimental PDB *and* AlphaFold structures, Bio2Byte biophysical properties painted
onto structure, residue interaction networks, and TM-align superposition.

- Binary: `TM-align` (committed at the repo root as a linux/amd64 ELF, which is why every
  image is `platform: linux/amd64`).
- APIs: Scop3P, UniProt REST, EBI Proteins (features + variation), PDBe (mappings +
  updated mmCIF), RCSB, AlphaFold DB.
- Browser libraries: py3Dmol, NGL, pyvis/vis.js.
- Notebook health: 28 cells, one of which is a 137 kB mega-cell holding ~130 functions,
  with two superseded "tombstone" cells. The Shiny port covers the protocol but the
  notebook itself is the least maintainable file in the repo.
- Features back-ported from the notebook after the first conversion: the **Include UniProt
  PTMs** source, **PDB entry pickers driven by the accession's UniProt cross-references**
  (replacing free-text entry codes, with a cascading chain picker that shows each chain's
  UniProt range), and **colouring the residue interaction network by a predicted Bio2Byte
  property**. The Bio2Byte tab itself stays as the Shiny app has it -- its table management
  is ahead of the notebook's.
- SIFTS numbering is ported: UniProt positions are translated onto each PDB entry's author
  numbering through the PDBe SIFTS API, falling back to the SIFTS-enriched mmCIF, then to a
  chain-range offset, then to direct numbering (which is the correct answer for AlphaFold).
  Whichever tier answered is shown to the user, and a site SIFTS cannot place is not drawn.
  Measured on `1A3N`, where the cleaved initiator Met offsets the numbering by one, this took
  correct mark placement from 1/19 to 19/19.

### Biophysical prediction & mutation effect — stage 4

Runs DynaMine, DisoMine and EFoldMine on a wild-type sequence and a mutant, then reports
categorical label shifts around the mutated position (rigid→flexible, ordered→disordered,
and so on). Thresholds: backbone `>1.0 / >0.8 / >0.69`, disorder `>0.50`, early folding
`>0.169`.

- In-process predictor: `b2bTools` 3.0.9rc3, pinned against CPU-only PyTorch
  (`torch==2.2.2+cpu`). The notebook needs a `torch.load` monkeypatch and a
  scikit-learn version-warning suppression that the Shiny app does not.
  See "Dependencies" in `apps/README.md` for why the torch flavour and version are
  what they are.
- APIs: UniProt FASTA, Scop3P `modifications`.
- Browser libraries: Bokeh, embedded via `CDN.render()` inside an iframe.

### Protein topology viewer — stage 4, and the best-engineered protocol here

Draws a 2D secondary-structure topology diagram beside a 3D viewer, with PTM and disease
variant overlays. Unlike every other protocol it is already a **package**
(`notebooks/topology_viewer/topology/`, 11 modules) with the notebook reduced to a
three-line launcher, a 1,100-assertion test script, and synthetic fixtures carrying
ground-truth secondary structure.

- Dependencies: **none required.** It parses PDB and mmCIF itself, and derives secondary
  structure through three interchangeable backends — `mkdssp` if on PATH, then `biotite`,
  then a from-scratch P-SEA implementation that needs nothing at all. `requests`,
  `numpy`/`biotite` and `ipywidgets` are all imported lazily inside the functions that
  need them.
- APIs (accession mode only): AlphaFold DB, PDBe (`best_structures`, SIFTS mappings,
  updated mmCIF), Scop3P, UniProt/EBI Proteins.
- Browser libraries: Mol\* pinned at `4.19.0` and NGL pinned at `2.3.1`, both from
  jsDelivr, behind one adapter interface. Its own diagram renderer is hand-written SVG.
- Why it ports cleanly: `build_view(...)` already returns a complete, ready-to-embed
  `<iframe srcdoc=...>` block, so the Shiny app only has to replace the `ipywidgets`
  control layer — the whole pipeline and all the browser code are reused untouched.
- **Known gap, quoted from its own README: the 3D viewer integration has never been
  executed.** The build environment had no network access to CDNs, so the Mol\* and NGL
  adapters are written but untested against a live library. Two specifics to watch when it
  first runs for real: `MolstarAdapter._buildIndex` depends on `OrderedSet.ofBounds` from
  the UMD bundle, which is not a formally guaranteed export (there is a correct-but-slow
  fallback, and the in-view status line says which path is live); and Mol\* is version-pinned
  precisely because of that dependency.
- Also note: the package README is stale relative to the code — it describes the UniProt→PDB
  dropdown and annotations as future phases, but `__version__` is `1.6.0-phase3` and
  `annotations.py` is 917 lines.

### RINAlign — stage 4

Builds residue interaction networks from two structures and compares them two ways. *Same
protein (diff)* matches residue positions directly and splits the contact sets into
conserved / lost / gained, plus a mutation list. *Different proteins (align)* is a real
graph alignment: node similarity is `0.25·residue-type + 0.25·degree + 0.5·Weisfeiler-Lehman
signature agreement` at WL depth 3, solved to a one-to-one mapping by Hungarian assignment,
then scored by edge Jaccard.

Networks use Cβ–Cβ contacts (Cα for glycine and where Cβ is absent) within a user-set 4–14 Å
cutoff, excluding sequence-adjacent pairs within a chain. `MSE`, `PTR`, `SEP` and `TPO` are
whitelisted so phospho-residues become real nodes.

- Binaries: **none.** This is the most portable of the unmigrated protocols.
- APIs: UniProt REST, PDBe (`best_structures`, SIFTS mappings), AlphaFold DB, RCSB, Scop3P,
  EBI Proteins.
- Browser libraries: D3 v7 from cdnjs and NGL from unpkg (unpinned `@latest`).
- Notebook health: cleanly layered — cells 2–5 are pure, dependency-light, testable
  functions, and all the risk is concentrated in one 50 kB cell of JavaScript held in Python
  f-strings. Five views: numeric summary, a canvas contact map with hand-rolled zoom/pan/pick,
  a static SVG network overlay, a D3 force layout, and a linked view that bridges the force
  graph to an NGL stage in both directions.
- Fixed during migration: per-click temp-directory leaks and missing HTTP timeouts;
  the three structure-discovery fallbacks now log which one failed instead of making a
  network outage look like "this protein has no PDB entries"; Compare is enabled only
  after both networks build; and changing a structure or the cutoff invalidates the
  built networks instead of silently comparing stale graphs. Contact detection moved to
  a KD-tree, with a test asserting the edge set is identical to the original O(n²) loop.
- Remaining gaps, both documented in [the spec](use-cases/rinalign.md): the contact map
  accepts PTM and variant positions but draws no overlay (the notebook computed the
  payloads and never interpolated them), and selenomethionine has no entry in the
  residue-similarity table, so `MSE` scores 0.0 against `MET` during cross-protein
  alignment.

### Scop3P → py3Dmol peptide mapper — not planned

Scientifically identical to the Scop3P peptide mapper, rendered with py3Dmol instead of
NGL. It is the cleanest notebook in the repo (a single closure, one API, four cells) and is
worth keeping as a teaching example, but a Shiny port would duplicate `apps/peptide_mapper`
rather than add a protocol.

### Modifications endpoint walkthrough — broken

The root `README.md` still links a Binder badge for `Scop3P_API.ipynb`, but that file was
removed when the notebooks were reorganised into `notebooks/`. Either restore the notebook
or drop the badge.

## What "migrated" requires

A protocol reaching stage 3–4 needs all of the following. Missing any one of them leaves it
unreachable or unshippable:

1. `apps/<name>/app.py` exposing a module-level `app = App(content_ui, server)`, with the
   science factored into `apps/common/<name>.py` so it is unit-testable without a browser.
2. Registration in `apps/portal/main.py` — one import plus one `APP_OPTIONS` entry — or the
   app is invisible in the portal image.
3. A `docker/apps/<name>.Dockerfile` (four lines, `FROM` the shared base) plus the name
   added to `APPS` in the `Makefile`. Ordering no longer matters: the Dockerfile you
   build is the app you get.
4. A `docker-compose.yml` service on the next free host port.
5. A `docs/use-cases/<name>.md` parity spec recording which notebook behaviours are
   preserved and which deliberately differ.
6. Unit tests for the pure logic and a construct/label assertion in
   `tests/integration/test_app_smoke.py`. The suite runs offline — network calls must be
   monkeypatched.
7. Entries in `apps/README.md` (apps list, ports, logs) and the root `README.md`.

CI publishes only the `scop3p-toolkit` target, so a protocol that reaches steps 1–2 ships
with the next release even without its own image.
