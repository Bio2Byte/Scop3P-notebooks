# Structure Visualisation: Use-Case Analysis and Functional Mapping

## Source Notebook
- `notebooks/Scop3P_PTM_structure_viz_voila_app.ipynb`

## Core Goal
Explore PTMs, disease variants, structures, biophysical predictions, residue interaction networks, and structural alignment for one UniProt protein.

## Use Cases
1. Set current UniProt accession for the session.
2. Fetch Scop3P PTMs and inspect a table.
3. Fetch disease-associated variants and inspect a table.
4. Render 3D structure from PDB or AlphaFold with PTM highlights.
5. Fetch UniProt sequence and run Bio2Byte predictions (backbone/disoMine/earlyFolding).
6. Show prediction table and render metric-driven 3D coloring on structure.
7. Build RIN from uploaded/downloaded PDB and visualize interactive pyvis network (with PTM/variant emphasis).
8. Upload two PDB files, optionally limit chain/range, run TM-align, inspect output and aligned structure view.

## Parity Notes
- `ipywidgets` tabs/buttons are mapped to Shiny tabs and action buttons.
- 3D viewers are rendered as HTML in Shiny outputs (NGL-based rendering).
- Data tables are rendered in scrollable HTML table blocks.
- TM-align remains a server-side dependency and is executed via subprocess.
