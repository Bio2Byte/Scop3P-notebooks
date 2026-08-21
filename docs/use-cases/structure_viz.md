# Structure Visualisation: Use-Case Analysis and Functional Mapping

## Source Notebook
- `notebooks/Scop3P_PTM_structure_viz_voila_app.ipynb`

## Core Goal
Explore PTMs, disease variants, structures, biophysical predictions, residue interaction networks, and structural alignment for one UniProt protein.

## Use Cases
1. Set the current UniProtKB accession for the session, which also loads the PDB entries
   cross-referenced from its UniProt record into every structure picker.
2. Fetch Scop3P PTMs and inspect a table.
3. Fetch disease-associated variants and inspect a table.
4. Render 3D structure from AlphaFold, or from a PDB entry chosen from the accession's own
   cross-references, with PTM highlights.
5. Fetch UniProt sequence and run Bio2Byte predictions (backbone/disoMine/earlyFolding).
6. Show prediction table and render metric-driven 3D coloring on structure.
7. Build RIN from uploaded/downloaded PDB and visualize interactive pyvis network, coloured
   either by PTM/variant status or by a predicted Bio2Byte property.
8. Compare two structures -- uploaded, or picked from the accession's PDB entries -- with
   TM-align, optionally limited by chain/range, and inspect the aligned view.

## Parity Notes
- `ipywidgets` tabs/buttons are mapped to Shiny tabs and action buttons.
- 3D viewers are rendered as HTML in Shiny outputs (NGL-based rendering).
- Data tables are rendered in scrollable HTML table blocks.
- TM-align remains a server-side dependency and is executed via subprocess.

## PTM sources

- Scop3P modifications, plus optional UniProt PTM features via the
  **Include UniProt PTMs** checkbox on tab 1 (ported from the notebook, which grew this
  source after the first Shiny conversion).
- Merged on accession + three-letter residue + position, so a site both sources describe
  appears once with the Scop3P naming and the references combined.
- `residue` is normalised to a three-letter code, with the descriptive Scop3P text kept
  in `modification`. The 3D viewer's colour map and the UniProt site key both require
  the code form.
- The Bio2Byte tab is deliberately unchanged: the Shiny implementation, including its
  table management, is ahead of the notebook's.

## Choosing a PDB entry

Every field that names a PDB entry -- on the structure tab, the RIN tab, and both sides of
the TM-align tab -- is a **dropdown of the structures cross-referenced from the accession's
own UniProt entry**, not a free-text box. Setting the protein fetches
`https://rest.uniprot.org/uniprotkb/{accession}.json` once and reads its
`uniProtKBCrossReferences`, so the user picks from what actually exists for that protein
instead of typing a four-character code and finding out later that it belongs to a
different one.

Each option carries the information needed to choose between 34 entries for the same
protein: `2IVS - X-ray - 2.00 A - chains A, B`. Picking an entry then cascades into its
**chain** picker, whose labels show the UniProt range that chain covers -- `A (705-1013)`.
That range is what lets a structure covering only the kinase domain be sliced correctly
instead of being read as though it started at residue 1.

Notes on the parsing, which is fussier than it looks:

- `A/B=705-1013` is one range shared by two chains, not a chain named `A/B`.
- `A=10-50; 70-120` is one chain observed in two segments; the picker keeps the outer span
  so the slice covers everything present.
- A cross-reference with no `Chains` property still lists its entry, with no chains.

The UniProt lookup uses a short connect timeout rather than the service default. It sits on
the critical path of **Set protein**, and a synchronous reactive effect blocks the ASGI loop
for *every* connected session while it waits -- an unreachable UniProt froze the app for a
full minute before this was bounded. A failed lookup logs a warning and leaves the pickers
empty; it never takes the app down.

## Colouring the network by a Bio2Byte property

The RIN tab's **Node colour** control has two modes:

- *Site status (PTM / variant)* -- the original behaviour.
- *Bio2Byte property, site status on the border* -- the node **fill** carries a predicted
  value on a green-to-red scale while the **border** keeps carrying PTM/variant status, so
  a residue's predicted flexibility and its modification state are legible at once.

The property list is driven by whatever the Bio2Byte tab has actually predicted, so it
cannot offer a column that does not exist. Before any prediction has run it reads
**Run Bio2Byte first** rather than rendering as an empty box, and the legend says where to
go. **Show RIN** recolours from the cached graph rather than recomputing contacts, so
switching property is a redraw, not a rebuild.

Values are matched on the frame's `Position` column, never by row order: a PDB chain does
not necessarily start at residue 1 or run without gaps, so a positional lookup would
silently shift every value down the sequence.

The legend's interpretation bands (`flexible` / `context dependent` / `rigid` /
`membrane spanning`, and the DisoMine and EFoldMine bands) share their boundaries with
`common.mutation_effect`, which labels the same predictions in its own inference step.
A value sitting exactly on a boundary belongs to the **lower** band, matching that module's
`value > threshold` tests. A unit test pins the two together by comparing the values at
which each side switches band, so the wording can differ per app but the science cannot
drift.

## Numbering: SIFTS maps UniProt positions onto the structure

PTM, variant and Bio2Byte positions are UniProt-numbered. A PDB entry numbers its residues
however the depositors chose, and the two frequently disagree, so every position is
translated through **SIFTS** (the EBI's residue-level correspondence between structure and
sequence) before anything is drawn.

This is not a theoretical concern. On `1A3N` (haemoglobin) the author numbering is offset by
one from UniProt's because the initiator methionine is cleaved in the mature protein, and
**1 of 19** PTM marks landed on the residue the modification actually names. After mapping,
**19 of 19** do. The failure it prevents is silent: an unmapped mark renders perfectly and is
simply wrong.

Four tiers, best first, with the one that answered recorded on the result:

| Tier | Source | `source` value |
|---|---|---|
| 1 | PDBe SIFTS API residue-level mapping segments | `sifts-api` |
| 2 | PDBe "updated" mmCIF, `_atom_site.pdbx_sifts_xref_db_num` | `sifts-mmcif` |
| 3 | Offset inferred from UniProt's chain range against the residues present | `chain-range-offset` |
| 4 | Numbering assumed to agree | `direct` |

Tier 4 is the *correct* answer for AlphaFold models, which are built on the UniProt sequence,
so it is not treated as a failure there.

Two details that are easy to get wrong, both pinned by test:

- PDBe reports `author_residue_number` as **null** at one end of most real segments (3 of the
  5 segments across `2IVT` and `1A3N`). The missing end is reconstructed from the other end
  plus the UniProt span, since a SIFTS segment is colinear by construction.
- The adjacent `residue_number` field must **never** stand in for it. That is label (entity)
  numbering: for `2IVT` it runs 4-314 against author numbering 703-1013, so substituting it
  would shift every position by about 700.

A UniProt position that SIFTS does not place in the structure is **not drawn**, and the count
is reported -- a mark in the wrong place is worse than a mark missing. This is why selecting
`2IVT` for `P07949` draws 9 of 36 sites: the rest lie outside the kinase domain the entry
contains.

### What the user is told

The method is named on screen, not left in the log, because a structure figure is evidence
and its numbering provenance is part of the claim. The 3D viewer and the RIN tab both carry a
note, and the status line repeats it:

> ✓ Positions mapped via SIFTS residue-level mapping (PDBe API) for 2IVT chain A: 311
> residues aligned to UniProt numbering. PTM, variant and property positions are shown at the
> residues SIFTS assigns them in this entry.

When SIFTS is unavailable and tier 3 answers, the note says so and warns that the placement
is inferred rather than authoritative. A test asserts that the phrase "mapped via SIFTS"
appears **only** for tiers 1 and 2, so a guess can never present itself as authoritative.

### Independent confirmation

On `2IVT`, SIFTS maps the UniProt phosphotyrosine at 905 to author residue 905, which the
crystallographers modelled as `PTR` -- phosphotyrosine. Two independent sources agreeing on a
modified residue is the strongest available check that the mapping is right. Tiers 1 and 2
were also cross-checked against each other on `2IVT` and `1A3N`: 286 and 141 shared residues
respectively, with zero disagreements.

## Shared UI conventions

This app uses the toolkit-wide vocabulary from `apps/common/ui_shell.py`: the accession
field is labelled **UniProtKB accession** (`ACCESSION_LABEL`), it sits in a
`scop3p_field_row` so its buttons share the input's baseline, and it carries a
**Load example** button wired to this app's `EXAMPLE_ACCESSION`. Result cards stretch to
the height of their controls card. See "Shared UI Vocabulary" in
[`apps/README.md`](../../apps/README.md).
