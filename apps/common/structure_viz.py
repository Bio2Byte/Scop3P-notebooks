from __future__ import annotations

import html
import json
import os
import subprocess
import tempfile
import uuid
from pathlib import Path

import networkx as nx
import pandas as pd
import pyvis.network
import requests
from Bio.PDB import PDBIO, PDBParser, Select
from b2bTools import SingleSeq, constants
from scipy.spatial import KDTree

from .services import Scop3PClient


B2B_METRIC_COLUMNS = (
    "backbone",
    "sidechain",
    "ppII",
    "coil",
    "sheet",
    "helix",
    "earlyFolding",
    "disoMine",
)
B2B_NORMALIZED_SUFFIX = "_normalized"


class StructureVizService:
    def __init__(self, workdir: Path, timeout: int = 60) -> None:
        self.workdir = workdir
        self.timeout = timeout
        self.scop3p_client = Scop3PClient(timeout=timeout)
        self.workdir.mkdir(parents=True, exist_ok=True)

    def fetch_ptms(self, accession: str) -> pd.DataFrame:
        dataframe = self.scop3p_client.fetch_modifications(accession)
        if dataframe.empty:
            return dataframe
        return dataframe

    def fetch_variants(self, accession: str) -> pd.DataFrame:
        url = f"https://www.ebi.ac.uk/proteins/api/variation/{accession}"
        response = requests.get(url, headers={"Accept": "application/json"}, timeout=self.timeout)
        response.raise_for_status()
        payload = response.json()
        rows: list[dict[str, object]] = []
        for feature in payload.get("features", []):
            if feature.get("type") != "VARIANT":
                continue
            for association in feature.get("association", []):
                if association.get("disease") is not True:
                    continue
                try:
                    position = int(feature.get("begin"))
                except Exception:
                    position = None
                rows.append(
                    {
                        "ACC_ID": accession,
                        "position": position,
                        "WT": feature.get("wildType"),
                        "MT": feature.get("mutatedType"),
                        "consequence": feature.get("consequenceType"),
                        "disease_name": association.get("name"),
                    }
                )
        return pd.DataFrame(rows)

    def fetch_sequence(self, accession: str) -> str:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        response = requests.get(url, timeout=self.timeout)
        response.raise_for_status()
        lines = response.text.splitlines()
        return "".join(line.strip() for line in lines if line and not line.startswith(">"))

    def predict_b2b(self, accession: str, sequence: str) -> pd.DataFrame:
        with tempfile.NamedTemporaryFile(prefix="seq_", suffix=".fasta", mode="w") as fasta_file:
            fasta_file.write(f">{accession}\n{sequence}\n")
            fasta_file.flush()
            predictor = SingleSeq(fasta_file.name)
            tools = []
            for name in ["TOOL_DYNAMINE", "TOOL_DISOMINE", "TOOL_EFOLDMINE"]:
                if hasattr(constants, name):
                    tools.append(getattr(constants, name))
            prediction = predictor.predict(tools=tools).get_all_predictions() if tools else predictor.predict().get_all_predictions()
        
        protein = prediction.get("proteins", {}).get(accession, {})
        return self._normalize_b2b_prediction(protein)

    @staticmethod
    def _normalize_b2b_prediction(protein: dict[str, object]) -> pd.DataFrame:
        import pprint;
        print("_normalize_b2b_prediction")
        pprint.pprint(protein, indent=4, sort_dicts=True)
        
        sequence = "".join(protein.get("seq", ""))
        print("SEQUENCE=", sequence, "length=", len(sequence))
        
        size = (
            len(sequence)
            or len(protein.get("backbone", []))
            or len(protein.get("sidechain", []))
            or len(protein.get("ppII", []))
            or len(protein.get("coil", []))
            or len(protein.get("sheet", []))
            or len(protein.get("helix", []))
            or len(protein.get("earlyFolding", []))
            or len(protein.get("disoMine", []))
        )
        print("SIZE=", size)

        def _coerce_series(value: object) -> list[object]:
            if isinstance(value, (list, tuple)):
                values = list(value)
            elif hasattr(value, "tolist"):
                values = list(value.tolist())  # type: ignore[call-arg]
            else:
                values = []

            if len(values) < size:
                values.extend([None] * (size - len(values)))
            elif len(values) > size:
                values = values[:size]
            return values

        dataframe = pd.DataFrame(
            {
                "Position": list(range(1, size + 1)),
                "Amino acid": protein.get("seq"),
                "backbone": _coerce_series(protein.get("backbone")),
                "sidechain": _coerce_series(protein.get("sidechain")),
                "ppII": _coerce_series(protein.get("ppII")),
                "coil": _coerce_series(protein.get("coil")),
                "sheet": _coerce_series(protein.get("sheet")),
                "helix": _coerce_series(protein.get("helix")),
                "earlyFolding": _coerce_series(protein.get("earlyFolding")),
                "disoMine": _coerce_series(protein.get("disoMine")),
            }
        )
        for column in B2B_METRIC_COLUMNS:
            dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce")
            dataframe[StructureVizService.b2b_metric_column(column, normalized=True)] = (
                StructureVizService._min_max_normalize_series(dataframe[column])
            )
        return dataframe

    @staticmethod
    def b2b_metric_column(metric: str, *, normalized: bool = False) -> str:
        return f"{metric}{B2B_NORMALIZED_SUFFIX}" if normalized else metric

    @staticmethod
    def _min_max_normalize_series(series: pd.Series) -> pd.Series:
        numeric = pd.to_numeric(series, errors="coerce")
        non_null = numeric.dropna()
        if non_null.empty:
            return pd.Series([pd.NA] * len(numeric), index=numeric.index, dtype="Float64")

        minimum = float(non_null.min())
        maximum = float(non_null.max())
        if minimum == maximum:
            return pd.Series(
                [0.0 if pd.notna(value) else pd.NA for value in numeric],
                index=numeric.index,
                dtype="Float64",
            )

        normalized = (numeric - minimum) / (maximum - minimum)
        return normalized.astype("Float64")

    def download_alphafold_pdb(self, accession: str) -> Path:
        out_path = self.workdir / f"AF-{accession}-F1-model_v6.pdb"
        url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.pdb"
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        out_path.write_bytes(response.content)
        return out_path

    def download_pdb(self, pdb_id: str) -> Path:
        pdb_key = pdb_id.strip().upper()
        out_path = self.workdir / f"{pdb_key}.pdb"
        if out_path.exists() and out_path.stat().st_size > 0:
            return out_path
        url = f"https://files.rcsb.org/download/{pdb_key}.pdb"
        response = requests.get(url, timeout=self.timeout)
        response.raise_for_status()
        out_path.write_bytes(response.content)
        return out_path

    def resolve_uploaded_or_remote_pdb(
        self,
        upload: list[dict[str, object]] | None,
        pdb_id: str | None = None,
        *,
        target_name: str | None = None,
    ) -> Path | None:
        if upload:
            upload_row = upload[0]
            source = Path(str(upload_row["datapath"]))
            filename = target_name or Path(str(upload_row["name"])).name
            target = self.workdir / filename
            target.write_bytes(source.read_bytes())
            return target

        if pdb_id and pdb_id.strip():
            downloaded = self.download_pdb(pdb_id)
            if target_name is None or downloaded.name == target_name:
                return downloaded
            target = self.workdir / target_name
            target.write_bytes(downloaded.read_bytes())
            return target

        return None


class ChainRangeSelect(Select):
    def __init__(self, chain_id: str, start: int | None, end: int | None) -> None:
        self.chain_id = chain_id
        self.start = start
        self.end = end

    def accept_residue(self, residue):  # noqa: ANN001
        if residue.parent.id != self.chain_id:
            return False
        seq_number = residue.id[1]
        if self.start is not None and seq_number < self.start:
            return False
        if self.end is not None and seq_number > self.end:
            return False
        return True


class StructureOps:
    @staticmethod
    def validate_pdb_id(pdb_id: str) -> str:
        pdb_key = pdb_id.strip().upper()
        if len(pdb_key) != 4 or not pdb_key.isalnum():
            raise ValueError(
                f"Invalid PDB ID '{pdb_id}'. Expected a 4-character RCSB identifier such as 2IVT."
            )
        return pdb_key

    @staticmethod
    def bfactor_pdb(pdb_path: Path, dataframe: pd.DataFrame, value_col: str, out_path: Path, chain: str | None = None) -> Path:
        values = dataframe[value_col].tolist()
        position_to_value = {index + 1: float(value) for index, value in enumerate(values) if value is not None and value == value}

        def set_bfactor(line: str, bfactor: float) -> str:
            return line[:60] + f"{bfactor:6.2f}" + line[66:]

        with pdb_path.open("r", encoding="utf-8", errors="ignore") as source, out_path.open("w", encoding="utf-8") as target:
            for line in source:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    target.write(line)
                    continue
                residue_chain = line[21].strip()
                if chain and residue_chain and residue_chain != chain:
                    target.write(line)
                    continue
                try:
                    residue_pos = int(line[22:26].strip())
                except Exception:
                    target.write(line)
                    continue
                target.write(set_bfactor(line, position_to_value[residue_pos]) if residue_pos in position_to_value else line)
        return out_path

    @staticmethod
    def chain_range_from_pdb(pdb_path: Path, chain_id: str) -> tuple[int, int]:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("chain", str(pdb_path))
        model = next(structure.get_models())
        chain = model[chain_id]
        positions = [res.id[1] for res in chain.get_residues() if res.id[0] == " "]
        return min(positions), max(positions)

    @staticmethod
    def chain_ranges_from_pdb(pdb_path: Path) -> dict[str, tuple[int, int]]:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("chains", str(pdb_path))
        model = next(structure.get_models())
        chain_ranges: dict[str, tuple[int, int]] = {}
        for chain in model.get_chains():
            positions = [res.id[1] for res in chain.get_residues() if res.id[0] == " "]
            if positions:
                chain_ranges[chain.id] = (min(positions), max(positions))
        return chain_ranges

    @staticmethod
    def save_chain_segment(source_pdb: Path, target_pdb: Path, chain_id: str, start: int | None, end: int | None) -> Path:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("segment", str(source_pdb))
        model = next(structure.get_models())
        if chain_id not in model:
            available = ", ".join(chain.id for chain in model.get_chains()) or "(none)"
            raise ValueError(
                f"Chain '{chain_id}' not found in {source_pdb.name}. Available chains: {available}."
            )

        chain = model[chain_id]
        positions = [res.id[1] for res in chain.get_residues() if res.id[0] == " "]
        if not positions:
            raise ValueError(f"Chain '{chain_id}' in {source_pdb.name} does not contain standard residues.")

        chain_start = min(positions)
        chain_end = max(positions)
        if start is not None and end is not None and start > end:
            raise ValueError(f"Invalid range for chain '{chain_id}': start {start} is greater than end {end}.")
        if start is not None and start > chain_end:
            raise ValueError(f"Start {start} is outside chain '{chain_id}' range {chain_start}-{chain_end}.")
        if end is not None and end < chain_start:
            raise ValueError(f"End {end} is outside chain '{chain_id}' range {chain_start}-{chain_end}.")

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(target_pdb), select=ChainRangeSelect(chain_id, start, end))
        segment_positions = [
            int(line[22:26].strip())
            for line in target_pdb.read_text(encoding="utf-8", errors="ignore").splitlines()
            if line.startswith(("ATOM", "HETATM")) and line[21].strip() == chain_id and line[22:26].strip()
        ]
        if not segment_positions:
            raise ValueError(
                f"Chain/range selection produced an empty segment for chain '{chain_id}' "
                f"with range {start or chain_start}-{end or chain_end} in {source_pdb.name}."
            )
        return target_pdb

    @staticmethod
    def run_tmalign(pdb1: Path, pdb2: Path, out_dir: Path, out_name: str = "aligned") -> tuple[Path, str]:
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = ["TM-align", str(pdb1.resolve()), str(pdb2.resolve()), "-o", out_name]
        try:
            process = subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as error:
            stderr = (error.stderr or "").strip()
            stdout = (error.stdout or "").strip()
            details = stderr or stdout or f"exit code {error.returncode}"
            raise RuntimeError(
                f"TM-align failed for {pdb1.name} vs {pdb2.name}: {details}"
            ) from error
        candidates = [
            out_dir / out_name,
            out_dir / f"{out_name}.pdb",
            out_dir / "TM_sup.pdb",
            out_dir / f"{out_name}.sup",
            out_dir / f"{out_name}.sup_all",
            out_dir / f"{out_name}.sup_all_atm",
        ]
        aligned_path = next((candidate for candidate in candidates if candidate.exists()), None)
        if aligned_path is None:
            files = ", ".join(sorted(path.name for path in out_dir.iterdir()))
            raise RuntimeError(f"No TM-align output found in {out_dir}. Files: {files}")
        return aligned_path, process.stdout

    @staticmethod
    def build_rin_graph(pdb_path: Path, chain: str = "A", cutoff: float = 8.0, atom_name: str = "CA") -> nx.Graph:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        model = next(structure.get_models())
        chain_obj = model[chain] if chain in model else next(model.get_chains())

        residues: list[int] = []
        coordinates: list[list[float]] = []
        for residue in chain_obj.get_residues():
            if residue.id[0] != " ":
                continue
            residue_number = residue.id[1]
            atom = residue[atom_name] if atom_name in residue else residue["CA"] if "CA" in residue else None
            if atom is None:
                continue
            residues.append(residue_number)
            coordinates.append(list(atom.get_coord()))

        graph = nx.Graph()
        for residue in residues:
            graph.add_node(residue)

        if not coordinates:
            return graph

        kdtree = KDTree(coordinates)
        pairs = kdtree.query_pairs(cutoff)
        for left, right in pairs:
            residue_left = residues[left]
            residue_right = residues[right]
            distance = float(((kdtree.data[left] - kdtree.data[right]) ** 2).sum() ** 0.5)
            graph.add_edge(residue_left, residue_right, distance=distance)

        graph.graph["chain"] = chain_obj.id
        return graph

    @staticmethod
    def rin_to_pyvis_html(
        graph: nx.Graph,
        out_html: Path,
        ptm_positions: list[int] | None = None,
        mutation_positions: list[int] | None = None,
    ) -> Path:
        ptm_set = set(ptm_positions or [])
        mutation_set = set(mutation_positions or [])
        network = pyvis.network.Network(height="650px", width="100%", bgcolor="#ffffff", font_color="#222")
        # Use CDN resources so the generated HTML does not depend on local /lib/bindings assets.
        network.cdn_resources = "remote"
        network.set_options(
            """
            {
              "nodes": {
                "borderWidth": 1
              },
              "interaction": {
                "hover": true,
                "selectConnectedEdges": true,
                "multiselect": false
              },
              "physics": {
                "stabilization": true
              }
            }
            """
        )
        network.from_nx(graph)
        for node in network.nodes:
            residue = int(node["id"])
            if residue in ptm_set and residue in mutation_set:
                color = "#9467bd"
                node["size"] = 40
            elif residue in ptm_set:
                color = "#1f77b4"
                node["size"] = 36
            elif residue in mutation_set:
                color = "#d62728"
                node["size"] = 36
            else:
                color = "#b0b0b0"
                node["size"] = 30
            node["color"] = {
                "background": color,
                "border": color,
                "highlight": {"background": color, "border": color},
                "hover": {"background": color, "border": color},
            }
        out_html.parent.mkdir(parents=True, exist_ok=True)
        network.write_html(str(out_html), open_browser=False, notebook=False)

        html_text = out_html.read_text(encoding="utf-8", errors="ignore")
        click_widget = """
<div id="rin-click-info" style="margin:8px 0 10px 0;padding:8px 10px;background:#f7f7f7;border:1px solid #ddd;border-radius:6px;font-family:monospace;font-size:12px;">
Click a node to inspect residue details.
</div>
<script type="text/javascript">
(() => {
  const info = document.getElementById("rin-click-info");
  const defaultMsg = "Click a node to inspect residue details.";
  const bindClick = () => {
    if (typeof network === "undefined" || !network) {
      window.setTimeout(bindClick, 80);
      return;
    }
    network.on("click", (params) => {
      if (!params.nodes || params.nodes.length === 0) {
        info.textContent = defaultMsg;
        return;
      }
      const node = params.nodes[0];
      let degreeText = "";
      try {
        if (typeof edges !== "undefined" && edges) {
          const degree = edges.get({
            filter: (edge) => edge.from === node || edge.to === node
          }).length;
          degreeText = ` | degree=${degree}`;
        }
      } catch (_err) {}
      info.textContent = `Selected residue: ${node}${degreeText}`;
    });
  };
  bindClick();
})();
</script>
"""
        if "</body>" in html_text:
            html_text = html_text.replace("</body>", click_widget + "\n</body>", 1)
            out_html.write_text(html_text, encoding="utf-8")

        return out_html


class StructureViewerBuilder:
    @staticmethod
    def ptm_html(pdb_text: str, accession: str, ptm_rows: list[dict[str, object]], chain: str | None = None) -> str:
        colors = {"SER": "#d62728", "THR": "#2ca02c", "TYR": "#ff7f0e"}
        payload = {"ptms": ptm_rows, "chain": chain or ""}
        safe_accession = html.escape(accession)
        uid = uuid.uuid4().hex
        panel_id = f"panel_{uid}"
        viewport_id = f"viewport_{uid}"
        return f"""<div style=\"position:relative;width:100%;height:700px;\"> 
  <div id=\"{panel_id}\" style=\"position:absolute;top:10px;left:10px;z-index:10;background:rgba(255,255,255,.9);padding:8px;border-radius:8px;\"> 
    <b>{safe_accession}</b><br/>PTM viewer (PDB/AlphaFold)
  </div>
  <div id=\"{viewport_id}\" style=\"width:100%;height:100%;\"></div>
</div>
<script>
(() => {{
  const pdbText = {json.dumps(pdb_text)};
  const payload = {json.dumps(payload)};
  const colorMap = {json.dumps(colors)};
  const stageEl = document.getElementById('{viewport_id}');
  const panelEl = document.getElementById('{panel_id}');
  if (!stageEl || !panelEl) return;

  const paint = () => {{
    const stage = new window.NGL.Stage(stageEl, {{ backgroundColor: 'white' }});
    const blob = new Blob([pdbText], {{type:'text/plain'}});
    stage.loadFile(blob, {{ext:'pdb'}}).then(comp => {{
      comp.addRepresentation('cartoon', {{color:'lightgrey'}});
      for (const row of payload.ptms) {{
        const pos = row.position;
        const residue = String(row.residue || '').trim();
        if (!pos) continue;
        const color = colorMap[residue] || '#7b241c';
        const chainPart = payload.chain ? ` AND :${{payload.chain}}` : '';
        const sele = `${{pos}}${{chainPart}}`;
        comp.addRepresentation('ball+stick', {{sele: sele, color: color}});
      }}
      comp.autoView();
    }});
  }};

  if (window.NGL) {{ paint(); return; }}
  const s = document.createElement('script');
  s.src = 'https://unpkg.com/ngl@latest/dist/ngl.js';
  s.onload = paint;
  s.onerror = () => {{ panelEl.innerHTML += '<div style="color:#b00020">Could not load NGL assets.</div>'; }};
  document.head.appendChild(s);
}})();
</script>
"""

    @staticmethod
    def b2b_html(pdb_text: str, accession: str, metric: str) -> str:
        safe_accession = html.escape(accession)
        uid = uuid.uuid4().hex
        panel_id = f"panel_{uid}"
        viewport_id = f"viewport_{uid}"
        return f"""<div style=\"position:relative;width:100%;height:700px;\"> 
  <div id=\"{panel_id}\" style=\"position:absolute;top:10px;left:10px;z-index:10;background:rgba(255,255,255,.9);padding:8px;border-radius:8px;\"> 
    <b>{safe_accession}</b><br/>Bio2Byte metric: {html.escape(metric)}
  </div>
  <div id=\"{viewport_id}\" style=\"width:100%;height:100%;\"></div>
</div>
<script>
(() => {{
  const pdbText = {json.dumps(pdb_text)};
  const stageEl = document.getElementById('{viewport_id}');
  const panelEl = document.getElementById('{panel_id}');
  if (!stageEl || !panelEl) return;
  const paint = () => {{
    const stage = new window.NGL.Stage(stageEl, {{ backgroundColor: 'white' }});
    const blob = new Blob([pdbText], {{type:'text/plain'}});
    stage.loadFile(blob, {{ext:'pdb'}}).then(comp => {{
      comp.addRepresentation('cartoon', {{colorScheme: 'bfactor'}});
      comp.autoView();
    }});
  }};
  if (window.NGL) {{ paint(); return; }}
  const s = document.createElement('script');
  s.src = 'https://unpkg.com/ngl@latest/dist/ngl.js';
  s.onload = paint;
  s.onerror = () => {{ panelEl.innerHTML += '<div style="color:#b00020">Could not load NGL assets.</div>'; }};
  document.head.appendChild(s);
}})();
</script>
"""
