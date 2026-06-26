# AlphaFold / PDBe-Style Topology Voila App

This project turns an AlphaFold DB mmCIF file into an interactive PDBe-style secondary-structure topology diagram inside a Voila/Jupyter app.

The app now derives topology from AlphaFold mmCIF data, converts it into the same API-shaped topology model used by PDBe (`helices`, `strands`, `coils`, `terms`, and path coordinates), and renders it with the official PDBe topology viewer plugin. DSSP flatfiles are still accepted as a fallback parser.

For AlphaFold models, helices and strands are read from mmCIF `_struct_conf`; strand-to-strand contacts are inferred from CA geometry and used to arrange sheet-like groups. The AlphaFold DB API is used to fetch the v6 mmCIF and metadata.

The topology is linked to Mol*, so hovering or clicking residues in the 2D topology can highlight/select the same residue in the 3D AlphaFold model. Residue-level pLDDT is exposed as a topology annotation/detail layer.

## Run

```powershell
pip install -r requirements.txt
voila dssp_topology_voila.ipynb
```

In JupyterLab, open `dssp_topology_voila.ipynb` and render it with Voila.
JupyterLab itself is optional and can be installed separately if you want a notebook editor.

## Use

1. Enter an **AFDB AC** value such as `P07949`.
2. Click **Fetch AFDB topology** to download the AlphaFold DB mmCIF and build the topology.
3. Hover or click topology residues to see pLDDT/secondary-structure details and link to Mol*.
4. Optionally upload a local `.cif`, `.mmcif`, or `.dssp` file and click **Visualize topology**.

If the original `P07949.dssp` path exists on the same machine, the app also enables **Load P07949 example**.

## Notes

PDBe’s exact production topology layout for PDB entries is precomputed server-side and served from the PDBe topology API. For AlphaFold entries there is no equivalent public PDBe topology endpoint, so this project generates compatible topology JSON locally from the AlphaFold mmCIF, then hands it to the official PDBe browser renderer.
