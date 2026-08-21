# Peptide Mapper: Use-Case Analysis and Functional Mapping

## Source Notebook
- `notebooks/Peptide_mapper_scop3p_voila.ipynb`

## Core User Goal
Map phosphopeptide coverage from Scop3P onto AlphaFold structures for a UniProtKB accession and inspect overlap and modification context.

## Use Cases

### 1. Load peptide/modification data by accession
- Input: UniProtKB accession (example: `O00571`, available via Load example).
- Action: call Scop3P peptide-modification endpoint.
- Output: peptide rows with start/end, modified residue position, and score.
- Failure paths:
  - Empty accession.
  - Scop3P request error.
  - No peptide rows returned.

### 2. Filter peptide candidates
- Search supported:
  - Sequence substring (case-insensitive).
  - Coordinate range (`70-90`).
  - Lower bound (`>=150`).
  - Upper bound (`<=300`).
  - Single residue (`154`).
- Behavior: filtered set updates selectable peptide entries.

### 3. Select mapping mode and peptide entries
- Modes:
  - `Unique peptide spans`: group by sequence/start/end.
  - `All rows`: one selectable entry per API row.
- Multi-select peptide entries used for mapping.
- `Map all (filtered)` selects every visible option.

### 4. Render mapped coverage on AlphaFold structure
- AlphaFold PDB download fallback: `v6`, then `v4`.
- Mapping behavior:
  - Protein cartoon context in grey.
  - Union peptide span cartoon in blue.
  - Intersection residues across selected peptides in red.
  - Modified sites in magenta (optional).
- Modification scope:
  - Selected peptides only.
  - All protein modifications.

### 5. Export styled session
- Exports rendered NGL HTML to `exports/<ACC>_styled_session.html`.
- Output preserves peptide/modification/intersection styling.
- Also downloadable in the browser: the styled HTML, the raw AlphaFold PDB, and a TSV
  of the mapped residues. The server-side `exports/` write is kept for parity, but on
  a shared server a fixed path is a cross-session collision and leaves user data on the
  host, so downloads are the preferred route.

### 6. Bring your own peptide table
Second source notebook: `notebooks/Peptide_mapper_fileupload_voila.ipynb`.

- Upload a TSV or CSV exported from a search engine. The delimiter is tab first, comma
  as a fallback; a single-column result is reported rather than left to fail later.
- Five columns are auto-detected — protein ID, peptide sequence, peptide start, peptide
  end, UniProt position — in three passes: case-insensitive exact match, then match on
  a punctuation-stripped key (so `Start position` matches `start`), then containment.
  Candidate names cover the notebook's own headers plus what MaxQuant, FragPipe and
  DIA-NN emit. Any field that cannot be identified is reported by name and left for the
  user to pick; it is never silently defaulted to the first column, which is what the
  notebook did.
- Protein identifiers are normalised: `sp|P07949|RET_HUMAN`, `P07949;Q12345` and padded
  values all reduce to the bare accession that AlphaFold's URL needs. Isoform suffixes
  such as `P07949-2` are kept. The notebook did none of this, so any FASTA-style
  identifier produced a download 404.
- Rows with non-numeric or reversed spans, zero-indexed starts, or missing fields are
  dropped, matching notebook parity.
- The protein picker lists each accession with its peptide and site counts.
- Everything downstream is shared with the Scop3P source: the search syntax, the
  selection list, mapping, the viewer and the exports.

## Functional Parity Mapping (Voila -> Shiny)
- `widgets.Text + Button` -> `input_text + input_action_button`.
- `ToggleButtons/Checkbox/SelectMultiple` -> `input_radio_buttons/input_checkbox/input_selectize`.
- `Output + display(view)` -> `output_ui` with embedded NGL HTML.
- Callback graph (`on_click`, `observe`) -> Shiny reactive effects/events.

## Peptide selector

"Map all (filtered)" can select every span at once (46 for O00571). Selectize renders
each selection as a chip on its own line, so the control grew to roughly 1400px and
pushed the structure viewer off screen. It is capped at 150px with internal scrolling,
long labels are clipped rather than widening the card, and a count below the control
reports how many peptides are selected now that the chips can scroll out of view.

## Shared UI conventions

This app uses the toolkit-wide vocabulary from `apps/common/ui_shell.py`: the accession
field is labelled **UniProtKB accession** (`ACCESSION_LABEL`), it sits in a
`scop3p_field_row` so its buttons share the input's baseline, and it carries a
**Load example** button wired to this app's `EXAMPLE_ACCESSION`. Result cards stretch to
the height of their controls card. See "Shared UI Vocabulary" in
[`apps/README.md`](../../apps/README.md).

## Accepted Deviations
- Viewer implementation is HTML-embedded NGL for Shiny compatibility instead of in-notebook `nglview` widget objects.
- UI layout follows Shiny conventions while preserving all notebook use cases and interaction semantics.
- The two sources share one app rather than being two apps. The tabs cover only how the
  peptide table is obtained; an uploaded table is normalised into exactly the column
  schema the Scop3P path produces, so `filter_peptides`, `build_options` and
  `map_selection` are reused unchanged.
- `build_options` omits the `maxScore=` label fragment when the table has no `score`
  column, which uploaded tables do not. Scop3P labels are byte-identical.
- Server code branches on the active source tab, never on which field happens to be
  filled: Shiny keeps an inactive tab's inputs readable, so a stale accession would
  otherwise override the uploaded protein.
- Session state is per-connection. It was previously a module-level singleton shared by
  every browser, which showed one user's peptide table to another — tolerable for
  Scop3P accessions, not for uploaded data.

## The peptide selector

Peptides are chosen from a multi-select dropdown rather than a list, because *Map all
(filtered)* can select 46 spans at once and selectize renders each as its own chip.

Two things about the chip needed fixing:

- **The remove control was unreachable.** Selectize renders the "x" as the last child inside
  the chip, and its own CSS makes the chip an inline-flex box. The label is then an anonymous
  flex item, which cannot shrink below its content width -- the flexbox `min-width:auto`
  trap -- so the "x" was pushed 215px past the clip edge on every one of 46 chips. They could
  be selected but never removed. The chip is now `inline-block` so `text-overflow` can
  ellipsize the label, with room reserved on the right for the "x" and the link positioned
  into it; `overflow:hidden` clips at the padding edge, so the reserved strip stays visible.
- **The label was longer than it needed to be.** `maxScore` carried full float precision
  (`278.273526557391`). It is now four significant figures (`278.3`) -- significant figures
  rather than fixed decimals, because the column holds both intensities in the hundreds and
  probabilities near zero, and a fixed `.2f` would render 1e-5 as "0.00" and lose it.
