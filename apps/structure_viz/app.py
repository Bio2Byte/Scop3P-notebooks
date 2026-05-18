from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pandas as pd
from shiny import App, reactive, render, ui

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from common.structure_viz import (  # noqa: E402
    B2B_METRIC_COLUMNS,
    StructureOps,
    StructureViewerBuilder,
    StructureVizService,
)
from common.logging_utils import get_logger  # noqa: E402
from common.ui_shell import scop3p_card, scop3p_shell, scop3p_footer  # noqa: E402


LOGGER = get_logger("scop3p.structure_viz")


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
        self.tm_input_1 = reactive.value(None)
        self.tm_input_2 = reactive.value(None)
        self.tm_chain_ranges_1 = reactive.value({})
        self.tm_chain_ranges_2 = reactive.value({})
        self.tm_structures_loaded = reactive.value(False)
        self.tm_loaded_signature_1 = reactive.value(None)
        self.tm_loaded_signature_2 = reactive.value(None)


controller = StructureVizController()


def _tm_source_signature(upload, pdb_id: str) -> tuple[str, str] | None:  # noqa: ANN001
    if upload:
        row = upload[0]
        datapath = str(row.get("datapath", ""))
        name = str(row.get("name", ""))
        return ("upload", f"{datapath}|{name}")

    pdb_key = pdb_id.strip().upper()
    if pdb_key:
        return ("pdb", pdb_key)
    return None


def _scroll_df(dataframe: pd.DataFrame) -> ui.Tag:
    if dataframe is None or dataframe.empty:
        return ui.p("No rows.")
    css = """
    <style>
      .scroll-df-wrap { max-height: 420px; overflow:auto; border:1px solid #ddd; border-radius:6px; }
      .scroll-df-wrap table { min-width:100%; border-collapse: collapse; font-size:13px; }
      .scroll-df-wrap th,.scroll-df-wrap td { padding:6px 8px; border-bottom:1px solid #eee; white-space:nowrap; text-align:center; }
      .scroll-df-wrap thead th { position: sticky; top:0; background:#fafafa; z-index:2; }
    </style>
    """
    return ui.HTML(css + f"<div class='scroll-df-wrap'>{dataframe.to_html(index=False, escape=False)}</div>")


def _reset_b2b_state() -> None:
    controller.sequence.set("")
    controller.b2b_df.set(pd.DataFrame())
    controller.b2b_html.set("")
    ui.update_select("b2b_metric", choices={}, selected=None)


def _b2b_metric_names(dataframe: pd.DataFrame) -> list[str]:
    if dataframe is None or dataframe.empty:
        return []
    return [
        metric
        for metric in B2B_METRIC_COLUMNS
        if StructureVizService.b2b_metric_column(metric) in dataframe.columns
    ]


def _selected_b2b_metric_column(metric: str | None, *, normalized: bool) -> str | None:
    if not metric:
        return None
    return StructureVizService.b2b_metric_column(metric, normalized=normalized)


def _b2b_table_dataframe(dataframe: pd.DataFrame, *, normalized: bool) -> pd.DataFrame:
    if dataframe is None or dataframe.empty:
        return pd.DataFrame()
    column_names = ["Position", "Amino acid"] + [
        StructureVizService.b2b_metric_column(metric, normalized=normalized)
        for metric in B2B_METRIC_COLUMNS
    ]
    table = dataframe.loc[:, column_names].copy()
    table.columns = ["Position", "Amino acid", *B2B_METRIC_COLUMNS]
    return table


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
                    ui.input_action_button("reset_b2b", "Reset results", class_="btn-secondary"),
                    col_widths=[3, 3, 3, 3],
                ),
                ui.input_checkbox("b2b_normalized", "Show normalized values", value=False),
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
                    ui.input_select("tm_chain1", "Chain 1", choices={}),
                    ui.input_numeric("tm_start1", "Start 1", value=None),
                    ui.input_numeric("tm_end1", "End 1", value=None),
                    ui.input_select("tm_chain2", "Chain 2", choices={}),
                    ui.input_numeric("tm_start2", "Start 2", value=None),
                    ui.input_numeric("tm_end2", "End 2", value=None),
                    col_widths=[2, 2, 2, 2, 2, 2],
                ),
                ui.layout_columns(
                    ui.input_action_button("load_tmalign_structures", "Load structures", class_="btn-warning"),
                    ui.output_ui("tm_actions"),
                    col_widths=[6, 6],
                ),
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
        LOGGER.info("set_accession requested accession=%s", accession or "-", extra={"event": "set_accession"})
        if not accession:
            return
        controller.accession.set(accession)
        _reset_b2b_state()
        controller.status.set(f"Protein set: {accession} | session: {controller.workdir}")

    @reactive.effect
    @reactive.event(input.fetch_ptm)
    def _fetch_ptm() -> None:
        accession = controller.accession.get()
        LOGGER.info("fetch_ptm requested accession=%s", accession or "-", extra={"event": "fetch_ptm"})
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return
        dataframe = controller.service.fetch_ptms(accession)
        controller.ptm_df.set(dataframe)
        controller.status.set(f"PTMs fetched: {len(dataframe)} rows.")
        LOGGER.info("fetch_ptm completed rows=%s", len(dataframe), extra={"event": "fetch_ptm"})

    @reactive.effect
    @reactive.event(input.fetch_variants)
    def _fetch_variants() -> None:
        accession = controller.accession.get()
        LOGGER.info("fetch_variants requested accession=%s", accession or "-", extra={"event": "fetch_variants"})
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return
        dataframe = controller.service.fetch_variants(accession)
        controller.var_df.set(dataframe)
        controller.status.set(f"Variants fetched: {len(dataframe)} rows.")
        LOGGER.info("fetch_variants completed rows=%s", len(dataframe), extra={"event": "fetch_variants"})

    @reactive.effect
    @reactive.event(input.fetch_af)
    def _fetch_af() -> None:
        accession = controller.accession.get()
        LOGGER.info("fetch_af requested accession=%s", accession or "-", extra={"event": "fetch_af"})
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return
        af_path = controller.service.download_alphafold_pdb(accession)
        controller.af_path.set(af_path)
        controller.status.set(f"AlphaFold downloaded: {af_path}")
        LOGGER.info("fetch_af completed path=%s", af_path, extra={"event": "fetch_af"})

    @reactive.effect
    @reactive.event(input.render_structure)
    def _render_structure() -> None:
        accession = controller.accession.get()
        LOGGER.info("render_structure requested accession=%s source=%s", accession or "-", input.structure_source(), extra={"event": "render_structure"})
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
        LOGGER.info("render_structure completed pdb=%s chain=%s", pdb_path, chain or "-", extra={"event": "render_structure"})

    @reactive.effect
    @reactive.event(input.fetch_seq)
    def _fetch_seq() -> None:
        accession = controller.accession.get()
        LOGGER.info("fetch_seq requested accession=%s", accession or "-", extra={"event": "fetch_seq"})
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return
        controller.b2b_html.set("")
        sequence = controller.service.fetch_sequence(accession)
        controller.sequence.set(sequence)
        controller.status.set(f"Sequence fetched: {len(sequence)} aa")
        LOGGER.info("fetch_seq completed length=%s", len(sequence), extra={"event": "fetch_seq"})

    @reactive.effect
    @reactive.event(input.run_b2b)
    def _run_b2b() -> None:
        accession = controller.accession.get()
        sequence = controller.sequence.get()
        LOGGER.info("run_b2b requested accession=%s sequence_length=%s", accession or "-", len(sequence), extra={"event": "run_b2b"})
        if not accession or not sequence:
            controller.status.set("Fetch sequence first.")
            return
        
        controller.status.set("Predicting biophysical features, please wait...")
        dataframe = controller.service.predict_b2b(accession, sequence)
        controller.b2b_df.set(dataframe)
        controller.b2b_html.set("")
        metrics = _b2b_metric_names(dataframe)
        ui.update_select("b2b_metric", choices={metric: metric for metric in metrics}, selected=metrics[0] if metrics else None)
        controller.status.set(f"Bio2Byte prediction completed ({len(dataframe)} rows).")
        LOGGER.info("run_b2b completed rows=%s metrics=%s", len(dataframe), len(metrics), extra={"event": "run_b2b"})

    @reactive.effect
    @reactive.event(input.render_b2b_3d)
    def _render_b2b() -> None:
        dataframe = controller.b2b_df.get()
        accession = controller.accession.get()
        metric = input.b2b_metric()
        normalized = bool(input.b2b_normalized())
        metric_column = _selected_b2b_metric_column(metric, normalized=normalized)
        af_path = controller.af_path.get()
        LOGGER.info(
            "render_b2b requested accession=%s metric=%s normalized=%s",
            accession or "-",
            metric or "-",
            normalized,
            extra={"event": "render_b2b"},
        )
        if dataframe is None or dataframe.empty or not metric or metric_column is None:
            controller.status.set("Run predictions and choose a metric first.")
            return
        if af_path is None:
            controller.b2b_html.set("")
            controller.status.set("Fetch AlphaFold first (tab 3).")
            return
        out_pdb = controller.workdir / f"b2b_{metric_column}.pdb"
        bfactor_pdb = StructureOps.bfactor_pdb(Path(af_path), dataframe, metric_column, out_pdb)
        html_payload = StructureViewerBuilder.b2b_html(
            pdb_text=bfactor_pdb.read_text(encoding="utf-8", errors="ignore"),
            accession=accession,
            metric=metric_column,
        )
        controller.b2b_html.set(html_payload)
        controller.status.set(
            f"Rendered Bio2Byte 3D metric: {metric}"
            f"{' (normalized)' if normalized else ''}"
        )
        LOGGER.info("render_b2b completed metric=%s normalized=%s", metric, normalized, extra={"event": "render_b2b"})

    @reactive.effect
    @reactive.event(input.reset_b2b)
    def _reset_b2b() -> None:
        LOGGER.info("reset_b2b requested", extra={"event": "reset_b2b"})
        _reset_b2b_state()
        controller.status.set("Bio2Byte results cleared.")

    @reactive.effect
    @reactive.event(input.rin_dl_af)
    def _rin_dl_af() -> None:
        accession = controller.accession.get()
        LOGGER.info("rin_dl_af requested accession=%s", accession or "-", extra={"event": "rin_dl_af"})
        if not accession:
            controller.status.set("Set a UniProt accession first.")
            return
        path = controller.service.download_alphafold_pdb(accession)
        controller.rin_path.set(path)
        controller.status.set(f"RIN input set to AlphaFold PDB: {path}")
        LOGGER.info("rin_dl_af completed path=%s", path, extra={"event": "rin_dl_af"})

    @reactive.effect
    @reactive.event(input.build_rin)
    def _build_rin() -> None:
        LOGGER.info("build_rin requested", extra={"event": "build_rin"})
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
        LOGGER.info(
            "build_rin completed pdb=%s nodes=%s edges=%s",
            pdb_path,
            graph.number_of_nodes(),
            graph.number_of_edges(),
            extra={"event": "build_rin"},
        )

    @reactive.effect
    @reactive.event(input.run_tmalign)
    def _run_tmalign() -> None:
        LOGGER.info("run_tmalign requested", extra={"event": "run_tmalign"})
        try:
            current_signature_1 = _tm_source_signature(input.tm_pdb1(), input.tm_pdb1_id().strip())
            current_signature_2 = _tm_source_signature(input.tm_pdb2(), input.tm_pdb2_id().strip())
            f1 = controller.tm_input_1.get()
            f2 = controller.tm_input_2.get()
            if f1 is None or f2 is None:
                controller.tm_report.set("Load both structures first.")
                controller.tm_html.set("")
                return
            if (
                current_signature_1 != controller.tm_loaded_signature_1.get()
                or current_signature_2 != controller.tm_loaded_signature_2.get()
            ):
                controller.tm_structures_loaded.set(False)
                controller.tm_report.set("TM-align inputs changed. Reload both structures first.")
                controller.tm_html.set("")
                return

            chain1 = input.tm_chain1() or "A"
            chain2 = input.tm_chain2() or "A"
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
            LOGGER.info("run_tmalign completed aligned=%s", aligned_path, extra={"event": "run_tmalign"})
        except Exception as error:
            LOGGER.exception("run_tmalign failed", extra={"event": "run_tmalign"})
            controller.tm_html.set("")
            controller.tm_report.set(f"TM-align error: {error}")
            controller.status.set("TM-align failed.")

    @reactive.effect
    @reactive.event(input.load_tmalign_structures)
    def _load_tmalign_structures() -> None:
        LOGGER.info("load_tmalign_structures requested", extra={"event": "load_tmalign_structures"})
        controller.tm_html.set("")
        controller.tm_structures_loaded.set(False)
        try:
            tm_pdb1_id = input.tm_pdb1_id().strip()
            tm_pdb2_id = input.tm_pdb2_id().strip()
            if tm_pdb1_id:
                StructureOps.validate_pdb_id(tm_pdb1_id)
            if tm_pdb2_id:
                StructureOps.validate_pdb_id(tm_pdb2_id)

            f1 = controller.service.resolve_uploaded_or_remote_pdb(
                input.tm_pdb1(),
                tm_pdb1_id,
                target_name="tm_input_1.pdb",
            )
            f2 = controller.service.resolve_uploaded_or_remote_pdb(
                input.tm_pdb2(),
                tm_pdb2_id,
                target_name="tm_input_2.pdb",
            )
            if f1 is None or f2 is None:
                controller.tm_input_1.set(None)
                controller.tm_input_2.set(None)
                controller.tm_chain_ranges_1.set({})
                controller.tm_chain_ranges_2.set({})
                controller.tm_report.set("Provide both structures via local upload or RCSB PDB ID.")
                return

            ranges_1 = StructureOps.chain_ranges_from_pdb(f1)
            ranges_2 = StructureOps.chain_ranges_from_pdb(f2)
            if not ranges_1 or not ranges_2:
                raise ValueError("Could not find any standard-residue chains in one of the loaded structures.")

            controller.tm_input_1.set(f1)
            controller.tm_input_2.set(f2)
            controller.tm_chain_ranges_1.set(ranges_1)
            controller.tm_chain_ranges_2.set(ranges_2)
            controller.tm_loaded_signature_1.set(_tm_source_signature(input.tm_pdb1(), tm_pdb1_id))
            controller.tm_loaded_signature_2.set(_tm_source_signature(input.tm_pdb2(), tm_pdb2_id))
            controller.tm_structures_loaded.set(True)

            first_chain_1 = next(iter(ranges_1))
            first_chain_2 = next(iter(ranges_2))
            ui.update_select("tm_chain1", choices={chain: chain for chain in ranges_1}, selected=first_chain_1)
            ui.update_select("tm_chain2", choices={chain: chain for chain in ranges_2}, selected=first_chain_2)
            start_1, end_1 = ranges_1[first_chain_1]
            start_2, end_2 = ranges_2[first_chain_2]
            ui.update_numeric("tm_start1", value=start_1, min=start_1, max=end_1)
            ui.update_numeric("tm_end1", value=end_1, min=start_1, max=end_1)
            ui.update_numeric("tm_start2", value=start_2, min=start_2, max=end_2)
            ui.update_numeric("tm_end2", value=end_2, min=start_2, max=end_2)
            controller.tm_report.set(
                "Loaded TM-align structures:\n"
                f"1) {f1.name}: chains {', '.join(ranges_1)}\n"
                f"2) {f2.name}: chains {', '.join(ranges_2)}"
            )
            controller.status.set("TM-align structures loaded.")
            LOGGER.info(
                "load_tmalign_structures completed structure1=%s chains1=%s structure2=%s chains2=%s",
                f1,
                list(ranges_1),
                f2,
                list(ranges_2),
                extra={"event": "load_tmalign_structures"},
            )
        except Exception as error:
            LOGGER.exception("load_tmalign_structures failed", extra={"event": "load_tmalign_structures"})
            controller.tm_input_1.set(None)
            controller.tm_input_2.set(None)
            controller.tm_chain_ranges_1.set({})
            controller.tm_chain_ranges_2.set({})
            controller.tm_loaded_signature_1.set(None)
            controller.tm_loaded_signature_2.set(None)
            controller.tm_report.set(f"TM-align load error: {error}")
            controller.status.set("TM-align structure load failed.")

    @reactive.effect
    def _invalidate_loaded_tmalign_inputs() -> None:
        current_signature_1 = _tm_source_signature(input.tm_pdb1(), input.tm_pdb1_id().strip())
        current_signature_2 = _tm_source_signature(input.tm_pdb2(), input.tm_pdb2_id().strip())
        loaded_signature_1 = controller.tm_loaded_signature_1.get()
        loaded_signature_2 = controller.tm_loaded_signature_2.get()

        if loaded_signature_1 is None and loaded_signature_2 is None:
            return
        if current_signature_1 == loaded_signature_1 and current_signature_2 == loaded_signature_2:
            return

        controller.tm_structures_loaded.set(False)
        controller.tm_input_1.set(None)
        controller.tm_input_2.set(None)
        controller.tm_chain_ranges_1.set({})
        controller.tm_chain_ranges_2.set({})
        controller.tm_loaded_signature_1.set(None)
        controller.tm_loaded_signature_2.set(None)
        controller.tm_html.set("")
        controller.tm_report.set("TM-align inputs changed. Reload both structures first.")
        LOGGER.info("tmalign inputs invalidated after source change", extra={"event": "load_tmalign_structures"})

    @reactive.effect
    def _sync_tm_chain1_range() -> None:
        chain_ranges = controller.tm_chain_ranges_1.get()
        selected_chain = input.tm_chain1()
        if not chain_ranges or selected_chain not in chain_ranges:
            return
        start, end = chain_ranges[selected_chain]
        current_start = input.tm_start1()
        current_end = input.tm_end1()
        next_start = start if current_start is None or current_start < start or current_start > end else current_start
        next_end = end if current_end is None or current_end < start or current_end > end else current_end
        if next_start > next_end:
            next_start, next_end = start, end
        ui.update_numeric("tm_start1", value=next_start, min=start, max=end)
        ui.update_numeric("tm_end1", value=next_end, min=start, max=end)

    @reactive.effect
    def _sync_tm_chain2_range() -> None:
        chain_ranges = controller.tm_chain_ranges_2.get()
        selected_chain = input.tm_chain2()
        if not chain_ranges or selected_chain not in chain_ranges:
            return
        start, end = chain_ranges[selected_chain]
        current_start = input.tm_start2()
        current_end = input.tm_end2()
        next_start = start if current_start is None or current_start < start or current_start > end else current_start
        next_end = end if current_end is None or current_end < start or current_end > end else current_end
        if next_start > next_end:
            next_start, next_end = start, end
        ui.update_numeric("tm_start2", value=next_start, min=start, max=end)
        ui.update_numeric("tm_end2", value=next_end, min=start, max=end)

    @reactive.effect
    @reactive.event(input.show_rin)
    def _show_rin() -> None:
        LOGGER.info("show_rin requested has_html=%s", bool(controller.rin_html.get()), extra={"event": "show_rin"})
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
        return _scroll_df(
            _b2b_table_dataframe(
                controller.b2b_df.get(),
                normalized=bool(input.b2b_normalized()),
            )
        )

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
    def tm_actions():
        return ui.input_action_button(
            "run_tmalign",
            "Align + Visualize",
            class_="btn-primary",
            disabled=not controller.tm_structures_loaded.get(),
        )

    @render.ui
    def tm_view():
        payload = controller.tm_html.get()
        if not payload:
            return ui.p("No TM-align structure view yet.")
        return ui.HTML(payload)


content_ui = ui.div(
    app_ui, scop3p_footer()
)

app = App(content_ui, server)
