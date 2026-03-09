from __future__ import annotations

import html
import json
from pathlib import Path


class NGLViewerBuilder:
    """Builds NGL HTML payload for inline render and export."""

    @staticmethod
    def build_html(
        *,
        accession: str,
        pdb_path: Path,
        union_ranges: list[tuple[int, int]],
        intersection_positions: list[int],
        modification_positions: list[int],
    ) -> str:
        pdb_text = pdb_path.read_text(encoding="utf-8", errors="ignore")
        payload = {
            "acc": accession,
            "union_ranges": union_ranges,
            "intersection": intersection_positions,
            "mods": modification_positions,
        }

        safe_accession = html.escape(accession)
        return f"""<!doctype html>
<html>
<head>
  <meta charset=\"utf-8\" />
  <title>{safe_accession} styled NGL session</title>
  <style>
    body {{ margin: 0; font-family: sans-serif; }}
    #ngl-root {{ position: relative; width: 100%; height: 700px; }}
    #viewport {{ width: 100%; height: 100%; }}
    #panel {{
      position: absolute; top: 10px; left: 10px; z-index: 10;
      background: rgba(255,255,255,0.9); padding: 10px; border-radius: 8px;
      max-width: 520px;
    }}
    code {{ font-size: 12px; }}
  </style>
</head>
<body>
  <div id=\"ngl-root\">
    <div id=\"panel\">
      <b>{safe_accession}</b><br/>
      <div>Grey: protein | Blue: peptides | Red: intersection | Magenta: mods</div>
      <div style=\"margin-top:6px;\"><code>union ranges: {len(union_ranges)} | mods: {len(modification_positions)} | intersection: {len(intersection_positions)}</code></div>
    </div>
    <div id=\"viewport\"></div>
  </div>

  <script>
    (() => {{
      const pdbText = {json.dumps(pdb_text)};
      const sessionPayload = {json.dumps(payload)};

      function rangesToSelection(ranges) {{
        if (!ranges || ranges.length === 0) return "";
        return ranges.map(r => `resi ${{r[0]}}-${{r[1]}}`).join(" OR ");
      }}

      function positionsToSelection(pos) {{
        if (!pos || pos.length === 0) return "";
        return pos.map(p => `resi ${{p}}`).join(" OR ");
      }}

      function showError(message) {{
        const panel = document.getElementById("panel");
        if (!panel) return;
        const err = document.createElement("div");
        err.style.marginTop = "8px";
        err.style.color = "#b00020";
        err.textContent = message;
        panel.appendChild(err);
      }}

      function renderStage() {{
        if (!window.NGL) {{
          showError("NGL failed to load in this browser session.");
          return;
        }}

        const stage = new window.NGL.Stage("viewport", {{ backgroundColor: "white" }});
        window.addEventListener("resize", () => stage.handleResize(), false);

        const blob = new Blob([pdbText], {{ type: "text/plain" }});
        stage.loadFile(blob, {{ ext: "pdb" }}).then(comp => {{
          comp.addRepresentation("cartoon", {{ color: "grey" }});

          const pepSel = rangesToSelection(sessionPayload.union_ranges);
          if (pepSel) {{
            comp.addRepresentation("cartoon", {{ sele: pepSel, color: "blue" }});
          }}

          const interSel = positionsToSelection(sessionPayload.intersection);
          if (interSel) {{
            comp.addRepresentation("ball+stick", {{ sele: interSel, color: "red" }});
          }}

          const modSel = positionsToSelection(sessionPayload.mods);
          if (modSel) {{
            comp.addRepresentation("ball+stick", {{ sele: modSel, color: "magenta" }});
          }}

          comp.autoView();
        }}).catch(() => {{
          showError("Unable to render the structure from the downloaded PDB.");
        }});
      }}

      if (window.NGL) {{
        renderStage();
        return;
      }}

      const nglScript = document.createElement("script");
      nglScript.src = "https://unpkg.com/ngl@latest/dist/ngl.js";
      nglScript.async = true;
      nglScript.onload = renderStage;
      nglScript.onerror = () => showError("Could not load NGL viewer assets.");
      document.head.appendChild(nglScript);
    }})();
  </script>
</body>
</html>
"""

    @staticmethod
    def export_html(path: Path, html_payload: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(html_payload, encoding="utf-8")
