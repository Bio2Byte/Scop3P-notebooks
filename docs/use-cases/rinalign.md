# RIN Alignment (RINAlign)

## Source

- Notebook: [`notebooks/RINAlign_align_and compare_networks.ipynb`](../../notebooks/RINAlign_align_and%20compare_networks.ipynb)
  (the filename contains a space)
- Logic: [`apps/common/rinalign.py`](../../apps/common/rinalign.py) — notebook cells 2–5
- Views: [`apps/common/rinalign_views.py`](../../apps/common/rinalign_views.py) — cell 6, verbatim
- App: [`apps/rinalign/app.py`](../../apps/rinalign/app.py) — cells 7–9

## Functional scope

- Load a UniProt accession and show protein name, gene, length and organism.
- Discover candidate structures through three fallbacks in order: UniProt PDB
  cross-references (already fetched, so free and most reliable), the PDBe
  `best_structures` endpoint, then PDBe SIFTS `mappings/uniprot`. The AlphaFold model
  is offered first when one exists.
- Pick a Left (Model A) and Right (Model B) structure independently.
- Build a residue interaction network per structure: Cβ–Cβ contacts (Cα for glycine
  and where Cβ is absent) within a 4–14 Å cutoff, default 8.0 Å, excluding
  sequence-adjacent pairs within a chain. `MSE`, `PTR`, `SEP` and `TPO` are treated as
  standard residues so phospho-residues become real network nodes.
- Two comparison modes:
  - **Same protein (diff)** — residue positions are matched directly and the contact
    sets differenced into conserved / lost / gained, plus a list of positions where
    the residue identity differs, plus a per-residue contact-impact table sorted by
    absolute net change, plus the Jaccard overlap.
  - **Different proteins (align)** — node similarity is
    `0.25·residue-type + 0.25·degree + 0.5·Weisfeiler–Lehman agreement` at WL depth 3,
    solved to a one-to-one node mapping by Hungarian assignment
    (`scipy.optimize.linear_sum_assignment`), then scored by edge Jaccard.
- Overlay Scop3P PTM positions (optionally union'd with UniProt PTM features) and
  UniProt disease-variant positions on the network views.
- Five views: numeric summary, interactive contact map, aligned network overlay,
  D3 force-directed network, and a linked network↔3D view.

## UI mapping to Shiny

| ipywidgets | Shiny |
|---|---|
| `mode_toggle` `ToggleButtons` | `ui.input_radio_buttons("mode")` |
| `uniprot_input` `Text` + `fetch_btn` | `scop3p_field_row(ui.input_text("accession", ACCESSION_LABEL), Fetch, Load example)` |
| `left_dd` / `right_dd` `Dropdown` | `ui.input_select` + a cascading `ui.update_select` pair fired from the Fetch handler |
| `cutoff_slider` `FloatSlider` | `ui.input_slider("cutoff", 4.0, 14.0, 8.0, step=0.5)` |
| `gen_btn` | `ui.input_action_button("generate")` |
| `cmp_btn` with `disabled` | `@render.ui def compare_actions()` returning a `disabled=` button |
| `ptm_uniprot_chk` `Checkbox` | `ui.input_checkbox("include_uniprot_ptms")` |
| `left_rin_out` / `right_rin_out` `Output` | `@render.ui` fed by `rin_html()` |
| `status_out` / `ptmvar_out` `Output` | `ui.output_text_verbatim` in `.scop3p-status` cards |
| `progress` `IntProgress` | the status string (no existing app has a progress bar) |
| `diff_out` — four stacked iframes | a radio selector plus one `@render.ui`, so exactly one heavy view is live |

The Shiny layout is two cards rather than the notebook's single column. **Controls** on
the left runs three numbered steps in the order the work happens: 1. structure
selection, 2. annotation sources (optional), 3. contact cutoff plus Generate and
Compare. **Networks** on the right carries the per-structure stats, the summary, the
view selector and the view itself.
The view selector is disabled until Compare / Align has produced a payload, and is
trimmed to the views the current mode actually produces.
| `app_state` module dict | per-session `reactive.value`s inside `server()` |
| `_voila_iframe()` | `ui.tags.iframe(srcdoc=...)` |

Structure dropdowns are keyed on index strings (`"0"`, `"1"`, …) with a
`structure_entries` map alongside, because the notebook stored whole dicts as
`Dropdown` values and Shiny select values must be strings.

## Behavior preserved

- Structure-discovery order and fallback precedence, the `resolution or 999` sort with
  unknown resolutions last, and `{pdb_id}_{chain}` de-duplication.
- `Chains` property parsing (`A/B=1-100, C=10-200`) from UniProt cross-references.
- Contact definition, glycine → Cα fallback, sequence-adjacency exclusion, distances
  rounded to two decimals.
- All diff metrics, and the fact that only positions present in **both** structures
  are compared — a contact is "lost" only when the residue exists on both sides,
  never merely because it is unresolved in one model.
- Alignment scoring weights, WL depth 3, and the top-20 node-mapping table.
- PTM and variant overlays take effect on the next Compare, matching the notebook's
  explicit instruction.
- Every line of view JavaScript, including the two-way
  `window.__RIN_HL` / `window.__RIN_ONSELECT` / `LVgo` bridge.

## Differences from Voilà

- `_voila_iframe` is not ported: Shiny escapes `srcdoc` itself, so the notebook's
  manual `html.escape` would double-escape. `rinalign_views.html_document()` supplies
  the same document shell without it.
- Each script-bearing view gets its **own** iframe, so their duplicated element ids
  and `window['<uid>*']` helpers cannot collide. `aligned_network_html` in particular
  uses entirely unprefixed globals and ids, so inlining it would be unsafe.
- The linked view remains a **single** iframe. Its D3 force graph and its NGL stage
  must share one `window` for the highlight handshake; splitting them breaks the
  bridge silently — nothing errors, clicks just stop crossing.
- Element ids are derived with md5 rather than `hash()`, which is `PYTHONHASHSEED`-salted
  and so changed on every interpreter restart.
- Contacts come from a `scipy.spatial.KDTree` rather than an O(n²) Python loop. The
  edge set is identical — `query_pairs(r)` selects `distance <= r`, the same
  comparison — and there is a test asserting equality against the original loop across
  three structures and three cutoffs. The motivation is that a `@reactive.effect` runs
  on the ASGI event loop, so a 2000-residue model would otherwise freeze every session.
- Structure downloads are cached in one per-session working directory, and every HTTP
  call has a timeout. The notebook called `tempfile.mkdtemp()` per click and never
  cleaned up, and its download and UniProt calls had no timeout at all.
- The three structure-discovery fallbacks log which one failed. The notebook's bare
  `except: pass` made a network outage indistinguishable from "this protein has no PDB
  entries".
- `align_rins` refuses networks above 1500 nodes with an actionable message rather
  than allocating a dense matrix and blocking the worker.
- Compare is enabled only after **both** networks build. The notebook enabled it once
  the left side succeeded and then raised on `GR.number_of_nodes()`.
- Changing a structure selection or the cutoff invalidates the built networks and
  disables Compare. The notebook read the cutoff only inside Generate, so results
  could silently disagree with the slider.
- Session state is per-connection.
- Enter-to-submit on the accession field is not reproduced.
- The BioPython import guard is dropped: it is a hard dependency, and the notebook's
  fallback printed to stdout, which would land in the structured log.

## Known gaps

- **The contact map draws no annotation overlay.** It accepts `ptm_pos` and `var_pos`
  and ignores them. The notebook computed `ptm_js`, `var_js` and a mutation list in
  that function and then never interpolated any of them into its script, which only
  receives positions, edges and residue names. The parameters are kept so all views
  share one signature. The other three views do mark PTMs and variants.
- **Selenomethionine has no similarity group.** `MSE` is in `THREE_TO_ONE` (so it
  becomes a network node) but was never added to `RESIDUE_GROUPS` alongside PTR/SEP/TPO.
  `restype_sim("MSE", "MET")` therefore returns 0.0 instead of 1.0. This only affects
  cross-protein alignment, where residue type is 25% of the node score, but
  selenomethionine is common in X-ray structures. Left unchanged because fixing it
  changes published alignment scores and needs a scientific decision.
- NGL is loaded from `unpkg.com/ngl@latest`, unpinned. This matches the existing
  convention in `apps/common/viewer.py` and `apps/common/structure_viz.py`, so it is
  not a new regression, but all three call sites would benefit from one pinned constant.
- The linked view embeds the left structure's coordinates in the `srcdoc`, roughly
  500 KB escaped for a 400-residue protein, sent over the websocket on each Compare.
  Same as the notebook, but now on a shared server; serving the structure from a
  session-scoped route would be the fix.

## Shared UI conventions

This app uses the toolkit-wide vocabulary from `apps/common/ui_shell.py`: the accession
field is labelled **UniProtKB accession** (`ACCESSION_LABEL`), it sits in a
`scop3p_field_row` so its buttons share the input's baseline, and it carries a
**Load example** button wired to this app's `EXAMPLE_ACCESSION`. Result cards stretch to
the height of their controls card. See "Shared UI Vocabulary" in
[`apps/README.md`](../../apps/README.md).

## Validation targets

- An unknown accession, and a protein with neither PDB entries nor an AlphaFold model,
  both produce a clear message.
- Generate fails cleanly — Compare stays disabled — if either download or either
  network build fails.
- Moving the cutoff slider after Generate disables Compare and says why.
- In align mode the four diff views explain that they need diff mode, rather than
  rendering nothing.
- PTM and variant fetches degrade to an empty set with an explanation; Scop3P mainly
  covers human phosphoproteins.
- In the linked view, clicking a network node highlights the residue in 3D **and**
  clicking an atom highlights the node — that round trip is the bridge.
- The linked view falls back to a visible message when NGL cannot load, and to a note
  when the left structure's coordinates were unavailable.
