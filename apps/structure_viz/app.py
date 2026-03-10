from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pandas as pd
from shiny import App, reactive, render, ui

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from apps.common.structure_viz import StructureOps, StructureViewerBuilder, StructureVizService  # noqa: E402
from apps.common.ui_shell import scop3p_card, scop3p_shell  # noqa: E402


class StructureVizController:
    def __init__(self) -> None:
        self.workdir = Path(tempfile.mkdtemp(prefix="scop3p_structure_viz_"))
        self.service = StructureVizService(self.workdir)

        self.accession = reactive.value("")
        self.status = reactive.value("Set a UniProt accession to start.")

        self.ptm_df = reactive.value(pd.DataFrame())
        self.var_df = reactive.value(pd.DataFrame())
        self.sequence = reactive.value("")
        self.b2b_df = reactive.value(pd.DataFrame())

        self.af_path = reactive.value(None)
        self.rin_path = reactive.value(None)

        self.viewer_html = reactive.value("")
        self.rin_html = reactive.value("")
        self.tm_report = reactive.value("")
        self.b2b_html = reactive.value("")
        self.tm_html = reactive.value("")


controller = StructureVizController()


def _scroll_df(dataframe: pd.DataFrame) -> ui.Tag:
    if dataframe is None or dataframe.empty:
        return ui.p("No rows.")
    css = """
    <style>
      .scroll-df-wrap { max-height: 420px; overflow:auto; border:1px solid #ddd; border-radius:6px; }
      .scroll-df-wrap table { min-width:100%; border-collapse: collapse; font-size:13px; }
      .scroll-df-wrap th,.scroll-df-wrap td { padding:6px 8px; border-bottom:1px solid #eee; white-space:nowrap; }
      .scroll-df-wrap thead th { position: sticky; top:0; background:#fafafa; z-index:2; }
    </style>
    """
    return ui.HTML(css + f"<div class='scroll-df-wrap'>{dataframe.to_html(index=False, escape=False)}</div>")


app_ui = scop3p_shell(
    "Structure Visualisation",
    "Inspect PTMs, disease variants, 3D structures, Bio2Byte overlays, residue interaction networks, and TM-align comparisons within one structure-centric workspace.",
    ui.div(
        scop3p_card(
            "Protein Setup",
            ui.layout_columns(
                ui.input_text("accession", "UniProt", value="", placeholder="e.g. P07949"),
                ui.input_action_button("set_accession", "Set protein", class_="btn-info"),
                col_widths=[8, 4],
            ),
        ),
        scop3p_card(
            "Session Status",
            ui.output_text_verbatim("status"),
            extra_class="scop3p-status",
        ),
        class_="scop3p-header-grid",
    ),
    ui.navset_tab(
        ui.nav_panel(
            "1) PTMs",
            scop3p_card(
                "PTM Table",
                ui.input_action_button("fetch_ptm", "Fetch PTMs", class_="btn-warning"),
                ui.output_ui("ptm_table"),
            ),
        ),
        ui.nav_panel(
            "2) Variants",
            scop3p_card(
                "Variant Table",
                ui.input_action_button("fetch_variants", "Fetch disease-associated variants", class_="btn-warning"),
                ui.output_ui("variant_table"),
            ),
        ),
        ui.nav_panel(
            "3) 3D Viewer",
            scop3p_card(
                "Structure Viewer",
                ui.layout_columns(
                    ui.input_radio_buttons("structure_source", "Source", {"pdb": "PDB", "af": "AlphaFold"}, selected="pdb", inline=True),
                    ui.input_text("pdb_id", "PDB ID", value="", placeholder="e.g. 2IVT"),
                    ui.input_text("chain", "Chain", value="", placeholder="A (optional)"),
                    col_widths=[4, 4, 4],
                ),
                ui.layout_columns(
                    ui.input_action_button("fetch_af", "Fetch AlphaFold", class_="btn-warning"),
                    ui.input_action_button("render_structure", "Show 3D", class_="btn-success"),
                    col_widths=[6, 6],
                ),
                ui.output_ui("structure_view"),
            ),
        ),
        ui.nav_panel(
            "4) Bio2Byte",
            scop3p_card(
                "Bio2Byte",
                ui.layout_columns(
                    ui.input_action_button("fetch_seq", "Fetch sequence", class_="btn-warning"),
                    ui.input_action_button("run_b2b", "Run predictions", class_="btn-danger"),
                    ui.input_action_button("render_b2b_3d", "Show 3D", class_="btn-success"),
                    col_widths=[4, 4, 4],
                ),
                ui.input_select("b2b_metric", "Color by", choices=[]),
                ui.output_ui("b2b_table"),
                ui.output_ui("b2b_view"),
            ),
        ),
        ui.nav_panel(
            "5) RIN",
            scop3p_card(
                "Residue Interaction Network",
                ui.layout_columns(
                    ui.input_action_button("rin_dl_af", "Download AlphaFold PDB", class_="btn-warning"),
                    ui.input_file("rin_upload", "Upload local PDB", accept=[".pdb"], multiple=False),
                    ui.input_text("rin_pdb_id", "Or fetch RCSB PDB ID", value="", placeholder="e.g. 2IVT"),
                    ui.input_text("rin_chain", "Chain", value="A"),
                    ui.input_slider("rin_cutoff", "Cutoff Å", min=4.0, max=12.0, value=8.0, step=0.5),
                    col_widths=[3, 3, 2, 2, 2],
                ),
                ui.layout_columns(
                    ui.input_action_button("build_rin", "Build RIN", class_="btn-danger"),
                    ui.input_action_button("show_rin", "Show RIN", class_="btn-success"),
                    col_widths=[6, 6],
                ),
                ui.output_ui("rin_view"),
            ),
        ),
        ui.nav_panel(
            "6) TM-align",
            scop3p_card(
                "TM-align",
                ui.layout_columns(
                    ui.input_file("tm_pdb1", "PDB 1: upload local", accept=[".pdb"], multiple=False),
                    ui.input_text("tm_pdb1_id", "PDB 1: or RCSB ID", value="", placeholder="e.g. 2IVT"),
                    col_widths=[6, 6],
                ),
                ui.layout_columns(
                    ui.input_file("tm_pdb2", "PDB 2: upload local", accept=[".pdb"], multiple=False),
                    ui.input_text("tm_pdb2_id", "PDB 2: or RCSB ID", value="", placeholder="e.g. 1CRN"),
                    col_widths=[6, 6],
                ),
                ui.layout_columns(
                    ui.input_text("tm_chain1", "Chain 1", value="A"),
                    ui.input_numeric("tm_start1", "Start 1", value=None),
                    ui.input_numeric("tm_end1", "End 1", value=None),
                    ui.input_text("tm_chain2", "Chain 2", value="A"),
                    ui.input_numeric("tm_start2", "Start 2", value=None),
                    ui.input_numeric("tm_end2", "End 2", value=None),
                    col_widths=[2, 2, 2, 2, 2, 2],
                ),
                ui.input_action_button("run_tmalign", "Align + Visualize", class_="btn-primary"),
                ui.output_text_verbatim("tm_output"),
                ui.output_ui("tm_view"),
            ),
        ),
    ),
)


def server(input, output, session):
    def require_accession() -> str | None:
        accession = input.accession().strip()
        if not accession:
            controller.status.set("Please enter a UniProt accession.")
            return None
        return accession

    @reactive.effect
    @reactive.event(input.set_accession)
    def _set_accession() -> None:
        accession = require_accession()
        if not accession:
            return
        controller.accession.set(accession)
        controller.status.set(f"Protein set: {accession} | session: {controller.workdir}")

    @reactive.effect
    @reactive.event(input.fetch_ptm)
    def _fetch_ptm() -> None:
        accession = controller.accession.get()
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return
        dataframe = controller.service.fetch_ptms(accession)
        controller.ptm_df.set(dataframe)
        controller.status.set(f"PTMs fetched: {len(dataframe)} rows.")

    @reactive.effect
    @reactive.event(input.fetch_variants)
    def _fetch_variants() -> None:
        accession = controller.accession.get()
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return
        dataframe = controller.service.fetch_variants(accession)
        controller.var_df.set(dataframe)
        controller.status.set(f"Variants fetched: {len(dataframe)} rows.")

    @reactive.effect
    @reactive.event(input.fetch_af)
    def _fetch_af() -> None:
        accession = controller.accession.get()
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return
        af_path = controller.service.download_alphafold_pdb(accession)
        controller.af_path.set(af_path)
        controller.status.set(f"AlphaFold downloaded: {af_path}")

    @reactive.effect
    @reactive.event(input.render_structure)
    def _render_structure() -> None:
        accession = controller.accession.get()
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return

        source = input.structure_source()
        chain = (input.chain() or "").strip().upper() or None

        if source == "af":
            pdb_path = controller.af_path.get()
            if pdb_path is None:
                controller.status.set("Fetch AlphaFold first.")
                return
        else:
            pdb_id = input.pdb_id().strip()
            if not pdb_id:
                controller.status.set("Provide a PDB ID for PDB source.")
                return
            pdb_path = controller.service.download_pdb(pdb_id)

        ptm_df = controller.ptm_df.get()
        ptm_rows = [] if ptm_df is None or ptm_df.empty else ptm_df[[c for c in ["position", "residue"] if c in ptm_df.columns]].fillna("").to_dict("records")
        html_payload = StructureViewerBuilder.ptm_html(
            pdb_text=Path(pdb_path).read_text(encoding="utf-8", errors="ignore"),
            accession=accession,
            ptm_rows=ptm_rows,
            chain=chain,
        )
        controller.viewer_html.set(html_payload)
        controller.status.set(f"Rendered 3D structure from: {Path(pdb_path).name}")

    @reactive.effect
    @reactive.event(input.fetch_seq)
    def _fetch_seq() -> None:
        accession = controller.accession.get()
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return
        sequence = controller.service.fetch_sequence(accession)
        controller.sequence.set(sequence)
        controller.status.set(f"Sequence fetched: {len(sequence)} aa")

    @reactive.effect
    @reactive.event(input.run_b2b)
    def _run_b2b() -> None:
        accession = controller.accession.get()
        sequence = controller.sequence.get()
        if not accession or not sequence:
            controller.status.set("Fetch sequence first.")
            return
        dataframe = controller.service.predict_b2b(accession, sequence)
        controller.b2b_df.set(dataframe)
        numeric = [column for column in dataframe.columns if pd.api.types.is_numeric_dtype(dataframe[column])]
        ui.update_select("b2b_metric", choices={column: column for column in numeric}, selected=numeric[0] if numeric else None)
        controller.status.set(f"Bio2Byte prediction completed ({len(dataframe)} rows).")

    @reactive.effect
    @reactive.event(input.render_b2b_3d)
    def _render_b2b() -> None:
        dataframe = controller.b2b_df.get()
        accession = controller.accession.get()
        metric = input.b2b_metric()
        af_path = controller.af_path.get()
        if dataframe is None or dataframe.empty or not metric:
            controller.status.set("Run predictions and choose a metric first.")
            return
        if af_path is None:
            controller.status.set("Fetch AlphaFold first (tab 3).")
            return
        out_pdb = controller.workdir / f"b2b_{metric}.pdb"
        bfactor_pdb = StructureOps.bfactor_pdb(Path(af_path), dataframe, metric, out_pdb)
        html_payload = StructureViewerBuilder.b2b_html(
            pdb_text=bfactor_pdb.read_text(encoding="utf-8", errors="ignore"),
            accession=accession,
            metric=metric,
        )
        controller.b2b_html.set(html_payload)
        controller.status.set(f"Rendered Bio2Byte 3D metric: {metric}")

    @reactive.effect
    @reactive.event(input.rin_dl_af)
    def _rin_dl_af() -> None:
        accession = controller.accession.get()
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return
        path = controller.service.download_alphafold_pdb(accession)
        controller.rin_path.set(path)
        controller.status.set(f"RIN input set to AlphaFold PDB: {path}")

    @reactive.effect
    @reactive.event(input.build_rin)
    def _build_rin() -> None:
        pdb_path = controller.service.resolve_uploaded_or_remote_pdb(
            input.rin_upload(),
            input.rin_pdb_id(),
        )
        if pdb_path is None and controller.rin_path.get() is not None:
            pdb_path = Path(controller.rin_path.get())

        if pdb_path is None:
            controller.status.set("Provide a local PDB upload, an RCSB PDB ID, or download AlphaFold first.")
            return

        graph = StructureOps.build_rin_graph(
            pdb_path,
            chain=(input.rin_chain() or "A").strip() or "A",
            cutoff=float(input.rin_cutoff()),
        )
        ptm_pos = []
        variant_pos = []
        if controller.ptm_df.get() is not None and not controller.ptm_df.get().empty and "position" in controller.ptm_df.get().columns:
            ptm_pos = controller.ptm_df.get()["position"].dropna().astype(int).tolist()
        if controller.var_df.get() is not None and not controller.var_df.get().empty and "position" in controller.var_df.get().columns:
            variant_pos = controller.var_df.get()["position"].dropna().astype(int).tolist()

        html_path = controller.workdir / f"rin_{controller.accession.get() or 'session'}.html"
        StructureOps.rin_to_pyvis_html(graph, html_path, ptm_pos, variant_pos)
        controller.rin_html.set(html_path.read_text(encoding="utf-8", errors="ignore"))
        controller.status.set(f"RIN built with {graph.number_of_nodes()} nodes / {graph.number_of_edges()} edges.")

    @reactive.effect
    @reactive.event(input.run_tmalign)
    def _run_tmalign() -> None:
        try:
            f1 = controller.service.resolve_uploaded_or_remote_pdb(
                input.tm_pdb1(),
                input.tm_pdb1_id(),
                target_name="tm_input_1.pdb",
            )
            f2 = controller.service.resolve_uploaded_or_remote_pdb(
                input.tm_pdb2(),
                input.tm_pdb2_id(),
                target_name="tm_input_2.pdb",
            )
            if f1 is None or f2 is None:
                controller.tm_report.set("Provide both structures via local upload or RCSB PDB ID.")
                controller.tm_html.set("")
                return

            chain1 = (input.tm_chain1() or "A").strip() or "A"
            chain2 = (input.tm_chain2() or "A").strip() or "A"
            start1 = int(input.tm_start1()) if input.tm_start1() is not None else None
            end1 = int(input.tm_end1()) if input.tm_end1() is not None else None
            start2 = int(input.tm_start2()) if input.tm_start2() is not None else None
            end2 = int(input.tm_end2()) if input.tm_end2() is not None else None

            seg1 = StructureOps.save_chain_segment(f1, controller.workdir / "seg1.pdb", chain1, start1, end1)
            seg2 = StructureOps.save_chain_segment(f2, controller.workdir / "seg2.pdb", chain2, start2, end2)
            aligned_path, report = StructureOps.run_tmalign(seg1, seg2, controller.workdir, out_name="aligned")

            first_line = report.splitlines()[0] if report.splitlines() else "TM-align completed."
            controller.tm_report.set(f"{first_line}\nAligned file: {aligned_path}")
            controller.tm_html.set(
                StructureViewerBuilder.ptm_html(
                    pdb_text=aligned_path.read_text(encoding="utf-8", errors="ignore"),
                    accession="TM-align result",
                    ptm_rows=[],
                    chain=None,
                )
            )
            controller.status.set("TM-align completed.")
        except Exception as error:
            controller.tm_html.set("")
            controller.tm_report.set(f"TM-align error: {error}")
            controller.status.set("TM-align failed.")

    @reactive.effect
    @reactive.event(input.show_rin)
    def _show_rin() -> None:
        if not controller.rin_html.get():
            controller.status.set("No RIN HTML yet. Build RIN first.")

    @render.text
    def status() -> str:
        return controller.status.get()

    @render.ui
    def ptm_table():
        return _scroll_df(controller.ptm_df.get())

    @render.ui
    def variant_table():
        return _scroll_df(controller.var_df.get())

    @render.ui
    def structure_view():
        payload = controller.viewer_html.get()
        if not payload:
            return ui.p("No structure rendered yet.")
        return ui.HTML(payload)

    @render.ui
    def b2b_table():
        return _scroll_df(controller.b2b_df.get())

    @render.ui
    def b2b_view():
        payload = controller.b2b_html.get()
        if not payload:
            return ui.p("No Bio2Byte 3D rendering yet.")
        return ui.HTML(payload)

    @render.ui
    def rin_view():
        if not controller.rin_html.get():
            return ui.p("No RIN built yet.")
        return ui.tags.iframe(
            srcdoc=controller.rin_html.get(),
            style="width:100%;height:700px;border:1px solid #ddd;border-radius:6px;",
        )

    @render.text
    def tm_output() -> str:
        return controller.tm_report.get()

    @render.ui
    def tm_view():
        payload = controller.tm_html.get()
        if not payload:
            return ui.p("No TM-align structure view yet.")
        return ui.HTML(payload)


app = App(app_ui, server)
