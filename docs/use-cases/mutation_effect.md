# Mutation Effect

## Source notebook
- `notebooks/Scop3P_b2b_mutation_effect_voila_app.ipynb`

## Functional scope
- Load a UniProt accession.
- Fetch the WT sequence from UniProt.
- Fetch Scop3P PTM annotations for the same accession.
- Run Bio2Byte single-sequence predictions for:
  - backbone dynamics
  - disorder (`disoMine`)
  - early folding
- Display WT predictions in an interactive line plot with PTM markers.
- Apply one or more user-defined mutations from comma-separated position and amino-acid inputs.
- Re-run Bio2Byte predictions on the mutant sequence.
- Overlay WT and mutant profiles in a second interactive plot.
- Show residue-level WT and WT-vs-mutant tables under the plots.
- Highlight mutated positions in the merged WT-vs-mutant table.
- Compute mutation-centric inference summaries for each predicted feature:
  - class at the mutation site
  - class over a +/-5 residue window
  - numeric deltas for site and window

## UI mapping to Shiny
- Notebook tab `WT prediction` -> Shiny tab `1) WT prediction`
- Notebook tab `Mutant prediction` -> Shiny tab `2) Mutant prediction`
- Notebook tab `Inference` -> Shiny tab `3) Inference`
- `ipywidgets.Text` for accession -> `ui.input_text("accession", ...)`
- `ipywidgets.Text` for mutation positions and target amino acids -> `ui.input_text(...)`
- `ipywidgets.Button` actions -> `ui.input_action_button(...)`
- Voilà output areas -> `render.ui` sections backed by reactive state

## Behavior preserved
- 1-indexed mutation application.
- Multi-mutation input via comma-separated values.
- WT and mutant predictions are calculated independently from the sequence given to `b2bTools.SingleSeq`.
- PTM markers are kept in both WT and mutant plots.
- Inference summaries keep the notebook thresholds:
  - backbone: `>1.0`, `>0.8`, `>0.69`
  - disorder: `>0.50`
  - early folding: `>0.169`

## Differences from Voilà
- Frontend is Shiny-native rather than `ipywidgets`.
- Bokeh plots are embedded in isolated iframes for reliable script execution in Shiny.
- Scrollable tables are rendered as HTML blocks rather than notebook display outputs.

## Validation targets
- Mutation parsing and application errors are surfaced clearly.
- Merged WT-vs-mutant table keeps mutated positions highlighted.
- Inference output is generated only after both WT and mutant predictions exist.
