# Peptide Mapper: Use-Case Analysis and Functional Mapping

## Source Notebook
- `notebooks/Peptide_mapper_scop3p_voila.ipynb`

## Core User Goal
Map phosphopeptide coverage from Scop3P onto AlphaFold structures for a UniProt accession and inspect overlap and modification context.

## Use Cases

### 1. Load peptide/modification data by accession
- Input: UniProt accession (example: `O00571`).
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

## Functional Parity Mapping (Voila -> Shiny)
- `widgets.Text + Button` -> `input_text + input_action_button`.
- `ToggleButtons/Checkbox/SelectMultiple` -> `input_radio_buttons/input_checkbox/input_selectize`.
- `Output + display(view)` -> `output_ui` with embedded NGL HTML.
- Callback graph (`on_click`, `observe`) -> Shiny reactive effects/events.

## Accepted Deviations
- Viewer implementation is HTML-embedded NGL for Shiny compatibility instead of in-notebook `nglview` widget objects.
- UI layout follows Shiny conventions while preserving all notebook use cases and interaction semantics.
