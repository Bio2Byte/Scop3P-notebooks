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
    build_ptm_table,
    merge_ptm_tables,
    b2b_legend_html,
    b2b_value_range,
    chain_choices_for_pdb,
    identity_mapping,
    numeric_b2b_columns,
    parse_tmalign_report,
    pdb_entry_choices,
    remap_positions,
    remap_site_rows,
    uniprot_range_for_chain,
)
from common.busy import (  # noqa: E402
    INERT_CSS,
    background,
    background_task_button,
    busy_indicators,
    finish_task,
    gate,
    task_button,
    task_outcome,
)
from common.vendor import enable_compression, static_assets  # noqa: E402
from common.logging_utils import get_logger, new_trail  # noqa: E402
from common.structure_labels import (  # noqa: E402
    ALL_CHAINS_PLACEHOLDER,
    NO_PROTEIN_PLACEHOLDER,
    NO_STRUCTURES_LOADED_PLACEHOLDER,
)
from common.ui_shell import (  # noqa: E402
    ACCESSION_LABEL,
    scop3p_card,
    scop3p_example_button,
    scop3p_field_row,
    scop3p_footer,
    scop3p_shell,
    scop3p_structure_picker,
)


LOGGER = get_logger("scop3p.structure_viz")

#: Shown in the Bio2Byte-property picker before any prediction has run.
B2B_METRIC_PLACEHOLDER = "Run Bio2Byte first"

#: Worked example: RET, with PTMs, disease variants and many PDB entries.
EXAMPLE_ACCESSION = "P07949"
EXAMPLE_PDB_ID = "2IVT"


class StructureVizController:
    def __init__(self) -> None:
        self.workdir = Path(tempfile.mkdtemp(prefix="scop3p_structure_viz_"))
        self.service = StructureVizService(self.workdir)

        self.accession = reactive.value("")
        self.status = reactive.value("Set a UniProtKB accession to start.")

        # PDB entries cross-referenced from the accession, and the per-chain UniProt
        # ranges that come with them. Populated when the protein is set, so every PDB
        # picker in the app offers only structures that actually contain this protein.
        self.pdb_xrefs = reactive.value([])

        self.ptm_df = reactive.value(pd.DataFrame())
        self.var_df = reactive.value(pd.DataFrame())
        self.sequence = reactive.value("")
        self.b2b_df = reactive.value(pd.DataFrame())

        self.af_path = reactive.value(None)
        self.rin_path = reactive.value(None)

        # How the current structure's residue numbers line up with UniProt's, and which
        # source said so. Shown to the user, because "SIFTS placed these marks" and
        # "we assumed the numbering agrees" are different claims about the same figure.
        self.position_mapping = reactive.value(None)

        self.viewer_html = reactive.value("")
        self.rin_html = reactive.value("")
        # The built network is cached so "Show RIN" can recolour without rebuilding it,
        # which is what makes flipping between Bio2Byte properties cheap.
        self.rin_graph = reactive.value(None)
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


def _tm_describe(pdb_id: str, chain: str, start, end) -> str:
    """Name one side of the comparison the way the user chose it."""
    bits = [pdb_id.strip().upper() or "uploaded file"]
    if chain:
        bits.append(f"chain {chain}")
    if start is not None and end is not None:
        bits.append(f"{int(start)}-{int(end)}")
    return " ".join(bits)


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
    # Render missing values as blank. pandas writes None as the literal string "None"
    # and a float NA as "NaN", which reads like data rather than absence -- and the
    # Scop3P v1 payload legitimately leaves evidence and reference null.
    display = dataframe.astype(object).where(dataframe.notna(), "")
    display = display.replace({None: "", "None": "", "nan": "", "<NA>": ""})
    return ui.HTML(css + f"<div class='scroll-df-wrap'>{display.to_html(index=False, escape=False)}</div>")


def _reset_b2b_state(controller: StructureVizController) -> None:
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
    ui.tags.style(
        """
.sv-button-row {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  margin-bottom: 12px;
}
"""
    ),
    ui.div(
        scop3p_card(
            "Protein Setup",
            scop3p_field_row(
                ui.input_text(
                    "accession",
                    ACCESSION_LABEL,
                    value="",
                    placeholder=f"e.g. {EXAMPLE_ACCESSION}",
                ),
                task_button(
                        "set_accession", "Set protein", class_="btn btn-info"),
                scop3p_example_button("load_example"),
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
                ui.div(
                    task_button(
                        "fetch_ptm", "Fetch PTMs", class_="btn btn-warning"),
                    class_="sv-button-row",
                ),
                ui.input_checkbox(
                    "include_uniprot_ptms",
                    "Include UniProt PTMs",
                    value=True,
                ),
                ui.p(
                    "Scop3P contributes experimentally observed modifications; UniProt "
                    "adds annotated PTMs of every kind (acetylation, methylation, "
                    "glycosylation and so on). A site described by both is listed once, "
                    "keeping the Scop3P naming and merging the references.",
                    class_="scop3p-note",
                ),
                ui.output_ui("ptm_table"),
            ),
        ),
        ui.nav_panel(
            "2) Variants",
            scop3p_card(
                "Variant Table",
                task_button(
                        "fetch_variants", "Fetch disease-associated variants", class_="btn-warning"),
                ui.output_ui("variant_table"),
            ),
        ),
        ui.nav_panel(
            "3) 3D Viewer",
            scop3p_card(
                "Structure Viewer",
                ui.layout_columns(
                    ui.input_radio_buttons("structure_source", "Source", {"pdb": "PDB", "af": "AlphaFold"}, selected="pdb", inline=True),
                    scop3p_structure_picker("pdb_id", "PDB entry", {"": NO_PROTEIN_PLACEHOLDER}),
                    ui.input_select("chain", "Chain", choices={"": ALL_CHAINS_PLACEHOLDER}),
                    col_widths=[4, 4, 4],
                ),
                ui.layout_columns(
                    task_button(
                        "fetch_af", "Fetch AlphaFold", class_="btn-warning"),
                    task_button(
                        "render_structure", "Show 3D", class_="btn-success"),
                    col_widths=[6, 6],
                ),
                ui.output_ui("numbering_note"),
                ui.output_ui("structure_view"),
            ),
        ),
        ui.nav_panel(
            "4) Bio2Byte",
            scop3p_card(
                "Bio2Byte",
                ui.layout_columns(
                    task_button(
                        "fetch_seq", "Fetch sequence", class_="btn-warning"),
                    background_task_button("run_b2b", "Run predictions", class_="btn-danger"),
                    task_button(
                        "render_b2b_3d", "Show 3D", class_="btn-success"),
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
                    task_button(
                        "rin_dl_af", "Download AlphaFold PDB", class_="btn-warning"),
                    ui.input_file("rin_upload", "Upload local PDB", accept=[".pdb"], multiple=False),
                    scop3p_structure_picker("rin_pdb_id", "Or a PDB entry", {"": NO_PROTEIN_PLACEHOLDER}),
                    ui.input_select("rin_chain", "Chain", choices={"A": "A"}),
                    ui.input_slider("rin_cutoff", "Cutoff Å", min=4.0, max=12.0, value=8.0, step=0.5),
                    col_widths=[3, 3, 2, 2, 2],
                ),
                # Node colouring. "Show RIN" re-renders from the cached graph, so
                # changing the property recolours without rebuilding the network.
                ui.layout_columns(
                    ui.input_select(
                        "rin_color_mode",
                        "Node colour",
                        choices={
                            "site": "Site status (PTM / variant)",
                            "b2b": "Bio2Byte property, site status on the border",
                        },
                        selected="site",
                    ),
                    ui.input_select(
                        "rin_b2b_metric",
                        "Bio2Byte property",
                        choices={"": B2B_METRIC_PLACEHOLDER},
                    ),
                    col_widths=[6, 6],
                ),
                ui.output_ui("rin_numbering_note"),
                ui.output_ui("rin_legend"),
                ui.layout_columns(
                    task_button(
                        "build_rin", "Build RIN", class_="btn-danger"),
                    task_button(
                        "show_rin", "Show RIN", class_="btn-success"),
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
                    scop3p_structure_picker("tm_pdb1_id", "PDB 1: or a PDB entry", {"": NO_PROTEIN_PLACEHOLDER}),
                    col_widths=[6, 6],
                ),
                ui.layout_columns(
                    ui.input_file("tm_pdb2", "PDB 2: upload local", accept=[".pdb"], multiple=False),
                    scop3p_structure_picker("tm_pdb2_id", "PDB 2: or a PDB entry", {"": NO_PROTEIN_PLACEHOLDER}),
                    col_widths=[6, 6],
                ),
                ui.layout_columns(
                    ui.input_select(
                        "tm_chain1",
                        "Chain 1",
                        choices={"": NO_STRUCTURES_LOADED_PLACEHOLDER},
                    ),
                    ui.input_numeric("tm_start1", "Start 1", value=None),
                    ui.input_numeric("tm_end1", "End 1", value=None),
                    ui.input_select(
                        "tm_chain2",
                        "Chain 2",
                        choices={"": NO_STRUCTURES_LOADED_PLACEHOLDER},
                    ),
                    ui.input_numeric("tm_start2", "Start 2", value=None),
                    ui.input_numeric("tm_end2", "End 2", value=None),
                    col_widths=[2, 2, 2, 2, 2, 2],
                ),
                ui.layout_columns(
                    task_button(
                        "load_tmalign_structures", "Load structures", class_="btn-warning"),
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
    # One trail per browser session: step numbers must not interleave across
    # sessions, and a module-level trail would be shared by every user.
    # One controller per session. It was created at module scope, which made every
    # browser session share the same accession, PTM table, sequence, Bio2Byte frame and
    # working directory -- one user's protein rendered in another user's browser. The
    # per-session extended tasks below would have crossed sessions the same way.
    controller = StructureVizController()

    trail = new_trail()
    trail.opened("Structure Visualisation")

    # Slow work runs in a worker thread. Nothing in here may touch reactive values or
    # `input`: it executes outside any reactive context. See common.busy.
    _b2b_task = background(
        lambda accession, sequence: controller.service.predict_b2b(accession, sequence)
    )

    def require_accession() -> str | None:
        accession = input.accession().strip()
        if not accession:
            trail.blocked("missing accession")
            controller.status.set("Please enter a UniProtKB accession.")
            return None
        return accession

    @reactive.effect
    @reactive.event(input.load_example)
    def _load_example() -> None:
        trail.clicked("Load example")
        ui.update_text("accession", value=EXAMPLE_ACCESSION)
        controller.status.set(
            f"Example accession {EXAMPLE_ACCESSION} loaded. Click Set protein."
        )

    @reactive.effect
    @reactive.event(input.set_accession)
    def _set_accession() -> None:
        accession = require_accession()
        # The value is recorded before the click so the record reads in the order the
        # user acted: they typed an accession, then pressed the button.
        trail.entered(ACCESSION_LABEL, accession)
        trail.clicked("Set protein")
        LOGGER.info("set_accession requested accession=%s", accession or "-", extra={"event": "set_accession"})
        if not accession:
            return
        controller.accession.set(accession)
        _reset_b2b_state(controller)

        # Offer the accession's own PDB entries everywhere an entry can be chosen.
        lookup_failed = False
        try:
            refs = controller.service.fetch_pdb_xrefs(accession)
        except Exception as error:  # noqa: BLE001 - the rest of the app still works
            LOGGER.warning(
                "pdb xref lookup failed accession=%s error=%s",
                accession, error, extra={"event": "set_accession"},
            )
            refs = []
            lookup_failed = True
        controller.pdb_xrefs.set(refs)

        # "0 entries" and "the lookup failed" must not read the same. The first is an
        # answer about the protein; the second is a broken run, and recording it at INFO
        # as a result would quietly turn an outage into a finding.
        if lookup_failed:
            trail.blocked(f"PDB entry lookup failed for {accession}; structure pickers are empty")
        else:
            trail.produced(f"{len(refs)} PDB entries cross-referenced from {accession}")
        entry_choices = pdb_entry_choices(refs, lookup_failed=lookup_failed)
        for input_id in ("pdb_id", "rin_pdb_id", "tm_pdb1_id", "tm_pdb2_id"):
            # update_selectize, not update_select: these are selectize widgets now, and
            # update_select would be accepted and then do nothing.
            ui.update_selectize(input_id, choices=entry_choices, selected="")

        suffix = (
            f" | {len(refs)} PDB entr{'y' if len(refs) == 1 else 'ies'}"
            if refs
            else " | no PDB entries cross-referenced"
        )
        controller.status.set(
            f"Protein set: {accession} | session: {controller.workdir}{suffix}"
        )

    @reactive.effect
    # ignore_init: reactive.event fires once for a select's initial value, which would
    # open every record with selections the user never made.
    @reactive.event(input.pdb_id, ignore_init=True)
    def _pdb_entry_changed() -> None:
        """Offer only the chains of this entry that carry the protein."""
        if input.pdb_id():
            trail.selected("PDB entry", input.pdb_id())
        refs = controller.pdb_xrefs.get()
        choices = chain_choices_for_pdb(refs, input.pdb_id())
        if not choices:
            ui.update_select("chain", choices={"": ALL_CHAINS_PLACEHOLDER}, selected="")
            return
        ui.update_select(
            "chain", choices={"": ALL_CHAINS_PLACEHOLDER, **choices}, selected=next(iter(choices))
        )

    @reactive.effect
    @reactive.event(input.rin_pdb_id, ignore_init=True)
    def _rin_pdb_entry_changed() -> None:
        if input.rin_pdb_id():
            trail.selected("PDB entry (RIN)", input.rin_pdb_id())
        refs = controller.pdb_xrefs.get()
        choices = chain_choices_for_pdb(refs, input.rin_pdb_id())
        if not choices:
            ui.update_select("rin_chain", choices={"A": "A"}, selected="A")
            return
        ui.update_select("rin_chain", choices=choices, selected=next(iter(choices)))

    @reactive.effect
    @reactive.event(input.chain, ignore_init=True)
    def _chain_changed() -> None:
        if input.chain():
            trail.selected("Chain", input.chain())

    @reactive.effect
    @reactive.event(input.rin_chain, ignore_init=True)
    def _rin_chain_changed() -> None:
        if input.rin_chain():
            trail.selected("Chain (RIN)", input.rin_chain())

    @reactive.effect
    @reactive.event(input.rin_color_mode, ignore_init=True)
    def _rin_colour_mode_changed() -> None:
        trail.selected("RIN node colour", input.rin_color_mode())

    @reactive.effect
    @reactive.event(input.rin_b2b_metric, ignore_init=True)
    def _rin_metric_changed() -> None:
        if input.rin_b2b_metric():
            trail.selected("RIN Bio2Byte property", input.rin_b2b_metric())

    @reactive.effect
    def _sync_rin_b2b_metrics() -> None:
        """Keep the property list in step with whatever Bio2Byte has predicted."""
        frame = controller.b2b_df.get()
        columns = numeric_b2b_columns(frame)
        if not columns:
            # An empty select renders as an empty box, which reads as a broken control
            # rather than as "nothing to choose yet".
            ui.update_select(
                "rin_b2b_metric",
                choices={"": B2B_METRIC_PLACEHOLDER},
                selected="",
            )
            return
        current = input.rin_b2b_metric()
        ui.update_select(
            "rin_b2b_metric",
            choices={column: column for column in columns},
            selected=current if current in columns else columns[0],
        )

    def _numbering_note_ui():
        """Tell the user, in the figure's own panel, how positions were placed.

        A structure figure is evidence, and its numbering provenance is part of the
        claim: marks placed by SIFTS and marks placed by assuming two numbering schemes
        agree deserve different confidence. So the method is named on screen rather than
        left in the log.
        """
        mapping = controller.position_mapping.get()
        if mapping is None:
            return ui.p(
                "Positions are UniProt-numbered. On a PDB entry they are translated to "
                "the structure's own numbering via SIFTS when the structure is drawn.",
                class_="scop3p-note",
            )
        if mapping.is_sifts:
            return ui.p(
                f"\u2713 {mapping.describe()} PTM, variant and property positions are "
                "shown at the residues SIFTS assigns them in this entry.",
                class_="scop3p-note",
            )
        if mapping.is_identity:
            return ui.p(
                f"{mapping.describe()} Correct for AlphaFold models, which are built on "
                "the UniProt sequence.",
                class_="scop3p-note",
            )
        return ui.p(
            f"\u26a0 {mapping.describe()} SIFTS was unavailable for this entry, so this "
            "is inferred rather than authoritative; check a known site before relying on "
            "the placement.",
            class_="scop3p-note",
        )

    @output(suspend_when_hidden=False)
    @render.ui
    def numbering_note():
        return _numbering_note_ui()

    @output(suspend_when_hidden=False)
    @render.ui
    def rin_numbering_note():
        return _numbering_note_ui()

    # See the note on variant_table: this output shares that nav_panel's suspension
    # problem, and without the guard the colour bar never appears.
    @output(suspend_when_hidden=False)
    @render.ui
    def rin_legend():
        """Colour bar for the selected property, with its interpretation bands."""
        if input.rin_color_mode() != "b2b":
            return None
        metric = input.rin_b2b_metric()
        frame = controller.b2b_df.get()
        if not metric or frame is None or frame.empty:
            return ui.p(
                "Run the Bio2Byte prediction on tab 4 first; its properties then become "
                "available here.",
                class_="scop3p-note",
            )
        span = b2b_value_range(frame, metric)
        if span is None:
            return ui.p(f"No numeric values for {metric}.", class_="scop3p-note")
        return ui.HTML(b2b_legend_html(metric, span[0], span[1]))

    @reactive.effect
    @reactive.event(input.fetch_ptm)
    def _fetch_ptm() -> None:
        accession = controller.accession.get()
        trail.clicked("Fetch PTMs")
        LOGGER.info("fetch_ptm requested accession=%s", accession or "-", extra={"event": "fetch_ptm"})
        if not accession:
            trail.blocked("accession not set")
            controller.status.set("Set a UniProtKB accession first.")
            return
        try:
            scop3p_table = build_ptm_table(
                controller.service.fetch_ptms(accession), accession, source="Scop3P"
            )
        except Exception as error:  # noqa: BLE001
            LOGGER.exception("fetch_ptm scop3p failed accession=%s", accession, extra={"event": "fetch_ptm"})
            trail.failed("Scop3P PTM fetch failed", error=type(error).__name__)
            controller.status.set(f"Scop3P PTM error: {error}")
            return

        uniprot_table = None
        uniprot_note = ""
        if input.include_uniprot_ptms():
            try:
                uniprot_table = controller.service.fetch_uniprot_ptms(accession)
            except Exception as error:  # noqa: BLE001
                # UniProt is the secondary source here, so losing it degrades the
                # result rather than failing it.
                LOGGER.warning(
                    "uniprot ptm fetch failed accession=%s error=%s",
                    accession,
                    error,
                    extra={"event": "fetch_ptm"},
                )
                uniprot_note = f" UniProt PTMs unavailable ({error})."

        dataframe = merge_ptm_tables(scop3p_table, uniprot_table)
        controller.ptm_df.set(dataframe)

        counts = f"Scop3P {len(scop3p_table)}"
        if uniprot_table is not None:
            counts += f" + UniProt {len(uniprot_table)}"
        controller.status.set(f"PTMs fetched: {len(dataframe)} sites ({counts}).{uniprot_note}")
        trail.produced(
            f"{len(dataframe)} PTM sites",
            scop3p=len(scop3p_table),
            uniprot="off" if uniprot_table is None else len(uniprot_table),
        )
        LOGGER.info(
            "fetch_ptm completed sites=%s scop3p=%s uniprot=%s",
            len(dataframe),
            len(scop3p_table),
            "off" if uniprot_table is None else len(uniprot_table),
            extra={"event": "fetch_ptm"},
        )

    @reactive.effect
    @reactive.event(input.fetch_variants)
    def _fetch_variants() -> None:
        accession = controller.accession.get()
        trail.clicked("Fetch disease-associated variants")
        LOGGER.info("fetch_variants requested accession=%s", accession or "-", extra={"event": "fetch_variants"})
        if not accession:
            trail.blocked("accession not set")
            controller.status.set("Set a UniProtKB accession first.")
            return
        dataframe = controller.service.fetch_variants(accession)
        controller.var_df.set(dataframe)
        controller.status.set(f"Variants fetched: {len(dataframe)} rows.")
        trail.produced(f"{len(dataframe)} disease-associated variants")
        LOGGER.info("fetch_variants completed rows=%s", len(dataframe), extra={"event": "fetch_variants"})

    @reactive.effect
    @reactive.event(input.fetch_af)
    def _fetch_af() -> None:
        accession = controller.accession.get()
        trail.clicked("Fetch AlphaFold")
        LOGGER.info("fetch_af requested accession=%s", accession or "-", extra={"event": "fetch_af"})
        if not accession:
            trail.blocked("accession not set")
            controller.status.set("Set a UniProtKB accession first.")
            return
        af_path = controller.service.download_alphafold_pdb(accession)
        controller.af_path.set(af_path)
        controller.status.set(f"AlphaFold downloaded: {af_path}")
        LOGGER.info("fetch_af completed path=%s", af_path, extra={"event": "fetch_af"})

    @reactive.effect
    @reactive.event(input.render_structure)
    def _render_structure() -> None:
        accession = controller.accession.get()
        trail.clicked("Show 3D (structure)")
        LOGGER.info("render_structure requested accession=%s source=%s", accession or "-", input.structure_source(), extra={"event": "render_structure"})
        if not accession:
            trail.blocked("accession not set")
            controller.status.set("Set a UniProtKB accession first.")
            return

        source = input.structure_source()
        chain = (input.chain() or "").strip().upper() or None

        if source == "af":
            pdb_path = controller.af_path.get()
            if pdb_path is None:
                trail.blocked("alphafold missing")
                controller.status.set("Fetch AlphaFold first.")
                return
        else:
            pdb_id = input.pdb_id().strip()
            if not pdb_id:
                trail.blocked("pdb id missing")
                controller.status.set("Provide a PDB ID for PDB source.")
                return
            pdb_path = controller.service.download_pdb(pdb_id)

        ptm_df = controller.ptm_df.get()
        ptm_rows = [] if ptm_df is None or ptm_df.empty else ptm_df[[c for c in ["position", "residue"] if c in ptm_df.columns]].fillna("").to_dict("records")

        # PTM positions are UniProt-numbered. A PDB entry numbers its residues however
        # the depositors chose, so they are translated through SIFTS before being handed
        # to the viewer -- otherwise a mark lands on whatever residue happens to carry
        # that number. AlphaFold models need no translation.
        mapping = _resolve_mapping(
            "" if source == "af" else pdb_id, chain or "", Path(pdb_path)
        )
        requested = len(ptm_rows)
        ptm_rows = remap_site_rows(ptm_rows, mapping)
        dropped = requested - len({row.get("uniprot_position", row.get("position")) for row in ptm_rows}) if requested else 0

        html_payload = StructureViewerBuilder.ptm_html(
            pdb_text=Path(pdb_path).read_text(encoding="utf-8", errors="ignore"),
            accession=accession,
            ptm_rows=ptm_rows,
            chain=chain,
        )
        controller.viewer_html.set(html_payload)
        note = f" {mapping.describe()}"
        if dropped > 0:
            note += (
                f" {dropped} site(s) are not present in this structure and are not drawn."
            )
        controller.status.set(f"Rendered 3D structure from: {Path(pdb_path).name}.{note}")
        trail.produced(
            f"3D structure rendered from {Path(pdb_path).name}",
            numbering=mapping.source, sites_drawn=len(ptm_rows), sites_not_in_structure=dropped,
        )
        LOGGER.info(
            "render_structure completed pdb=%s chain=%s mapping=%s sites=%s dropped=%s",
            pdb_path, chain or "-", mapping.source, len(ptm_rows), dropped,
            extra={"event": "render_structure"},
        )

    @reactive.effect
    @reactive.event(input.fetch_seq)
    def _fetch_seq() -> None:
        accession = controller.accession.get()
        trail.clicked("Fetch sequence")
        LOGGER.info("fetch_seq requested accession=%s", accession or "-", extra={"event": "fetch_seq"})
        if not accession:
            trail.blocked("accession not set")
            controller.status.set("Set a UniProtKB accession first.")
            return
        controller.b2b_html.set("")
        sequence = controller.service.fetch_sequence(accession)
        controller.sequence.set(sequence)
        controller.status.set(f"Sequence fetched: {len(sequence)} aa")
        trail.produced(f"sequence fetched, {len(sequence)} residues")
        LOGGER.info("fetch_seq completed length=%s", len(sequence), extra={"event": "fetch_seq"})

    @reactive.effect
    @reactive.event(input.run_b2b)
    def _run_b2b() -> None:
        accession = controller.accession.get()
        sequence = controller.sequence.get()
        trail.clicked("Run Bio2Byte predictions")
        LOGGER.info("run_b2b requested accession=%s sequence_length=%s", accession or "-", len(sequence), extra={"event": "run_b2b"})
        if not accession or not sequence:
            trail.blocked("sequence missing")
            controller.status.set("Fetch sequence first.")
            # The click already disabled the button in the browser; without this it stays
            # dead for the rest of the session.
            finish_task("run_b2b")
            return

        controller.status.set(
            f"Predicting biophysical features for {accession} "
            f"({len(sequence)} residues). This takes a few seconds."
        )
        _b2b_task(accession, sequence)

    @reactive.effect
    def _run_b2b_done() -> None:
        """Apply the prediction once the worker thread has finished with it."""

        def succeeded(dataframe) -> None:  # noqa: ANN001
            controller.b2b_df.set(dataframe)
            controller.b2b_html.set("")
            metrics = _b2b_metric_names(dataframe)
            ui.update_select(
                "b2b_metric",
                choices={metric: metric for metric in metrics},
                selected=metrics[0] if metrics else None,
            )
            controller.status.set(f"Bio2Byte prediction completed ({len(dataframe)} rows).")
            trail.produced(
                f"Bio2Byte predictions over {len(dataframe)} residues",
                properties=len(metrics),
            )
            LOGGER.info(
                "run_b2b completed rows=%s metrics=%s",
                len(dataframe), len(metrics), extra={"event": "run_b2b"},
            )

        def failed(error: Exception) -> None:
            LOGGER.exception("run_b2b failed", exc_info=error, extra={"event": "run_b2b"})
            trail.failed("run_b2b failed", error=type(error).__name__)
            controller.status.set(f"Bio2Byte prediction error: {error}")

        task_outcome(
            _b2b_task,
            on_success=succeeded,
            on_error=failed,
            on_finished=lambda: finish_task("run_b2b"),
        )

    @reactive.effect
    @reactive.event(input.render_b2b_3d)
    def _render_b2b() -> None:
        dataframe = controller.b2b_df.get()
        accession = controller.accession.get()
        metric = input.b2b_metric()
        normalized = bool(input.b2b_normalized())
        metric_column = _selected_b2b_metric_column(metric, normalized=normalized)
        af_path = controller.af_path.get()
        trail.clicked("Show 3D (Bio2Byte)")
        LOGGER.info(
            "render_b2b requested accession=%s metric=%s normalized=%s",
            accession or "-",
            metric or "-",
            normalized,
            extra={"event": "render_b2b"},
        )
        if dataframe is None or dataframe.empty or not metric or metric_column is None:
            trail.blocked("prediction or metric missing")
            controller.status.set("Run predictions and choose a metric first.")
            return
        if af_path is None:
            trail.blocked("alphafold missing")
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
        trail.clicked("Reset Bio2Byte results")
        LOGGER.info("reset_b2b requested", extra={"event": "reset_b2b"})
        _reset_b2b_state(controller)
        controller.status.set("Bio2Byte results cleared.")

    @reactive.effect
    @reactive.event(input.rin_dl_af)
    def _rin_dl_af() -> None:
        accession = controller.accession.get()
        trail.clicked("Download AlphaFold PDB for the RIN")
        LOGGER.info("rin_dl_af requested accession=%s", accession or "-", extra={"event": "rin_dl_af"})
        if not accession:
            trail.blocked("accession not set")
            controller.status.set("Set a UniProtKB accession first.")
            return
        path = controller.service.download_alphafold_pdb(accession)
        controller.rin_path.set(path)
        controller.status.set(f"RIN input set to AlphaFold PDB: {path}")
        LOGGER.info("rin_dl_af completed path=%s", path, extra={"event": "rin_dl_af"})

    def _resolve_mapping(pdb_id: str, chain: str, pdb_path: Path | None):
        """The UniProt-to-author numbering map for a structure, recorded for the UI.

        AlphaFold models are built on the UniProt sequence, so their numbering already
        agrees and an identity mapping is the correct answer rather than a fallback.
        """
        accession = controller.accession.get()
        if not pdb_id:
            mapping = identity_mapping()
        else:
            mapping = controller.service.build_position_mapping(
                pdb_id,
                accession,
                chain,
                uniprot_range=uniprot_range_for_chain(
                    controller.pdb_xrefs.get(), pdb_id, chain
                ),
                pdb_path=pdb_path,
            )
        controller.position_mapping.set(mapping)
        return mapping

    def _remap_frame_positions(frame, mapping):  # noqa: ANN001
        """A copy of a UniProt-numbered frame renumbered into author numbering.

        The Bio2Byte frame is indexed by UniProt position while the network's nodes carry
        author numbers, so without this the property lookup misses on exactly the entries
        where the two numberings differ -- which is every entry this mapping exists for.
        """
        if frame is None or frame.empty or mapping.is_identity:
            return frame
        if "Position" not in frame.columns:
            return frame
        renumbered = frame.copy()
        renumbered["Position"] = [
            (mapping.to_pdb(position) or [None])[0]
            for position in pd.to_numeric(renumbered["Position"], errors="coerce").fillna(-1).astype(int)
        ]
        return renumbered.dropna(subset=["Position"])

    def _site_positions(mapping=None) -> tuple[list[int], list[int]]:
        ptm_frame = controller.ptm_df.get()
        variant_frame = controller.var_df.get()
        ptm_positions: list[int] = []
        variant_positions: list[int] = []
        if ptm_frame is not None and not ptm_frame.empty and "position" in ptm_frame.columns:
            ptm_positions = ptm_frame["position"].dropna().astype(int).tolist()
        if variant_frame is not None and not variant_frame.empty and "position" in variant_frame.columns:
            variant_positions = variant_frame["position"].dropna().astype(int).tolist()
        if mapping is not None and not mapping.is_identity:
            # Sites are UniProt-numbered; the network's nodes are author-numbered.
            ptm_positions = remap_positions(ptm_positions, mapping)
            variant_positions = remap_positions(variant_positions, mapping)
        return ptm_positions, variant_positions

    def _render_rin(graph) -> None:  # noqa: ANN001
        """Write the network HTML using the colour settings as they stand."""
        mapping = controller.position_mapping.get() or identity_mapping()
        ptm_positions, variant_positions = _site_positions(mapping)
        by_property = input.rin_color_mode() == "b2b"
        metric = input.rin_b2b_metric() if by_property else None
        frame = _remap_frame_positions(controller.b2b_df.get(), mapping) if by_property else None

        stem = f"rin_{controller.accession.get() or 'session'}"
        if by_property and metric:
            stem = f"{stem}_{str(metric).replace('/', '_').replace(' ', '_')}"
        html_path = controller.workdir / f"{stem}.html"

        StructureOps.rin_to_pyvis_html(
            graph,
            html_path,
            ptm_positions,
            variant_positions,
            b2b_frame=frame,
            b2b_metric=metric,
        )
        controller.rin_html.set(html_path.read_text(encoding="utf-8", errors="ignore"))

    @reactive.effect
    @reactive.event(input.build_rin)
    def _build_rin() -> None:
        trail.clicked("Build RIN")
        LOGGER.info("build_rin requested", extra={"event": "build_rin"})
        pdb_path = controller.service.resolve_uploaded_or_remote_pdb(
            input.rin_upload(),
            input.rin_pdb_id(),
        )
        if pdb_path is None and controller.rin_path.get() is not None:
            pdb_path = Path(controller.rin_path.get())

        if pdb_path is None:
            trail.blocked("pdb source missing")
            controller.status.set("Provide a local PDB upload, an RCSB PDB ID, or download AlphaFold first.")
            return

        chain = (input.rin_chain() or "A").strip() or "A"

        # When the structure came from a PDB entry, build the network over the residue
        # range UniProt says that chain covers. Otherwise a network for a 2000-residue
        # asymmetric unit is built to describe a 300-residue protein, and the site
        # overlays land on residues belonging to something else.
        # Only claim a SIFTS mapping when the structure really is that PDB entry. An
        # upload or an AlphaFold model can sit alongside a stale entry in the picker, and
        # mapping one structure's numbering with another entry's SIFTS data would be
        # worse than not mapping at all.
        upload = input.rin_upload()
        from_pdb_entry = not upload and pdb_path == controller.service.pdb_path_for(
            input.rin_pdb_id()
        )
        rin_pdb_id = input.rin_pdb_id().strip() if from_pdb_entry else ""
        span = uniprot_range_for_chain(controller.pdb_xrefs.get(), rin_pdb_id, chain)
        mapping = _resolve_mapping(rin_pdb_id, chain, pdb_path)

        rin_input = pdb_path
        range_note = ""
        if span is not None:
            # save_chain_segment cuts on the numbers written in the file, which are author
            # numbers. Translate UniProt's span through SIFTS first; slicing on the raw
            # UniProt bounds only works where the two numberings happen to coincide.
            author_bounds = remap_positions([span[0], span[1]], mapping)
            if author_bounds:
                low, high = min(author_bounds), max(author_bounds)
            else:
                low, high = span
            trimmed = controller.workdir / f"rin_{pdb_path.stem}_{chain}_{low}_{high}.pdb"
            try:
                rin_input = StructureOps.save_chain_segment(pdb_path, trimmed, chain, low, high)
                range_note = f" (chain {chain}, UniProt {span[0]}-{span[1]})"
                if (low, high) != (span[0], span[1]):
                    range_note += f" = residues {low}-{high} in {rin_pdb_id}"
            except Exception as error:  # noqa: BLE001 - fall back to the whole file
                LOGGER.warning(
                    "chain range extraction failed pdb=%s chain=%s error=%s",
                    pdb_path, chain, error, extra={"event": "build_rin"},
                )

        graph = StructureOps.build_rin_graph(
            rin_input, chain=chain, cutoff=float(input.rin_cutoff())
        )
        controller.rin_graph.set(graph)
        _render_rin(graph)
        controller.status.set(
            f"RIN built with {graph.number_of_nodes()} nodes / "
            f"{graph.number_of_edges()} edges{range_note}. {mapping.describe()}"
        )
        trail.produced(
            f"RIN built for chain {chain}",
            nodes=graph.number_of_nodes(), edges=graph.number_of_edges(),
            cutoff=float(input.rin_cutoff()), numbering=mapping.source,
        )
        LOGGER.info(
            "build_rin completed pdb=%s chain=%s range=%s nodes=%s edges=%s",
            pdb_path,
            chain,
            span,
            graph.number_of_nodes(),
            graph.number_of_edges(),
            extra={"event": "build_rin"},
        )

    @reactive.effect
    @reactive.event(input.run_tmalign)
    def _run_tmalign() -> None:
        trail.clicked("Run TM-align")
        LOGGER.info("run_tmalign requested", extra={"event": "run_tmalign"})
        try:
            current_signature_1 = _tm_source_signature(input.tm_pdb1(), input.tm_pdb1_id().strip())
            current_signature_2 = _tm_source_signature(input.tm_pdb2(), input.tm_pdb2_id().strip())
            f1 = controller.tm_input_1.get()
            f2 = controller.tm_input_2.get()
            if f1 is None or f2 is None:
                trail.blocked("structures not loaded")
                controller.tm_report.set("Load both structures first.")
                controller.tm_html.set("")
                return
            if (
                current_signature_1 != controller.tm_loaded_signature_1.get()
                or current_signature_2 != controller.tm_loaded_signature_2.get()
            ):
                controller.tm_structures_loaded.set(False)
                trail.blocked("stale inputs")
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
            alignment = StructureOps.run_tmalign(
                seg1, seg2, controller.workdir, out_name="aligned"
            )
            aligned_path, report = alignment.superposed, alignment.report

            # The scores are the result. Previously only report.splitlines()[0] was shown,
            # which is TM-align's blank first line, so the user saw a temp-file path and
            # nothing about how similar the two structures are.
            result = parse_tmalign_report(report)
            controller.tm_report.set(
                f"{result.summary()}\n\n"
                f"Aligned structure: {aligned_path.name}\n\n"
                f"--- full TM-align output ---\n{report.strip()}"
            )
            # Both structures, distinctly coloured. Rendering only the superposed one
            # showed a single shape, which is not a superposition and cannot answer the
            # question the protocol is for.
            controller.tm_html.set(
                StructureViewerBuilder.superposition_html(
                    superposed_text=alignment.superposed.read_text(
                        encoding="utf-8", errors="ignore"
                    ),
                    reference_text=alignment.reference.read_text(
                        encoding="utf-8", errors="ignore"
                    ),
                    label_superposed=_tm_describe(
                        input.tm_pdb1_id(), chain1, start1, end1
                    )
                    + " - superposed",
                    label_reference=_tm_describe(
                        input.tm_pdb2_id(), chain2, start2, end2
                    )
                    + " - reference",
                )
            )
            controller.status.set("TM-align completed.")
            LOGGER.info("run_tmalign completed aligned=%s", aligned_path, extra={"event": "run_tmalign"})
        except Exception as error:
            LOGGER.exception("run_tmalign failed", extra={"event": "run_tmalign"})
            trail.failed("run_tmalign failed", error=type(error).__name__)
            controller.tm_html.set("")
            controller.tm_report.set(f"TM-align error: {error}")
            controller.status.set("TM-align failed.")

    @reactive.effect
    @reactive.event(input.load_tmalign_structures)
    def _load_tmalign_structures() -> None:
        trail.clicked("Load TM-align structures")
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
                trail.blocked("missing structure input")
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
            ui.update_select(
                "tm_chain1",
                choices={chain: chain for chain in ranges_1}
                or {"": NO_STRUCTURES_LOADED_PLACEHOLDER},
                selected=first_chain_1,
            )
            ui.update_select(
                "tm_chain2",
                choices={chain: chain for chain in ranges_2}
                or {"": NO_STRUCTURES_LOADED_PLACEHOLDER},
                selected=first_chain_2,
            )
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
            trail.failed("load_tmalign_structures failed", error=type(error).__name__)
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
        """Re-render the cached network with the current colour settings.

        Recolouring does not need the graph rebuilt, so switching Bio2Byte property is
        a redraw rather than a re-run of the contact calculation.
        """
        trail.clicked("Show RIN")
        graph = controller.rin_graph.get()
        if graph is None:
            trail.blocked("no cached graph")
            controller.status.set("No RIN yet. Click Build RIN first.")
            return

        by_property = input.rin_color_mode() == "b2b"
        metric = input.rin_b2b_metric()
        frame = controller.b2b_df.get()
        if by_property and (not metric or frame is None or frame.empty):
            controller.status.set(
                "Run the Bio2Byte prediction on tab 4 before colouring the network by a "
                "property."
            )
            return

        try:
            _render_rin(graph)
        except Exception as error:  # noqa: BLE001
            LOGGER.exception("show_rin failed", extra={"event": "show_rin"})
            trail.failed("show_rin failed", error=type(error).__name__)
            controller.status.set(f"RIN render error: {error}")
            return

        controller.status.set(
            f"RIN recoloured by {metric}." if by_property else "RIN coloured by site status."
        )
        LOGGER.info(
            "show_rin completed mode=%s metric=%s",
            input.rin_color_mode(),
            metric or "-",
            extra={"event": "show_rin"},
        )

    @render.text
    def status() -> str:
        return controller.status.get()

    @render.ui
    def ptm_table():
        return _scroll_df(controller.ptm_df.get())

    # suspend_when_hidden=False because this output lives in a nav_panel that is not
    # the initially-active tab. Shiny decides suspension from the client-reported
    # ".clientdata_output_<id>_hidden" value, and Session._is_hidden() treats "never
    # reported" as hidden, so such an output is suspended at page load and is never
    # woken when the user opens its tab: it sits at "recalculating" forever with no
    # error logged anywhere. Verified against shiny 1.7.0, which requirements-shiny.txt
    # permits (shiny>=1.1,<2).
    @output(suspend_when_hidden=False)
    @render.ui
    def variant_table():
        return _scroll_df(controller.var_df.get())

    @output(suspend_when_hidden=False)
    @render.ui
    def structure_view():
        payload = controller.viewer_html.get()
        if not payload:
            return ui.p("No structure rendered yet.")
        return ui.HTML(payload)

    @output(suspend_when_hidden=False)
    @render.ui
    def b2b_table():
        return _scroll_df(
            _b2b_table_dataframe(
                controller.b2b_df.get(),
                normalized=bool(input.b2b_normalized()),
            )
        )

    @output(suspend_when_hidden=False)
    @render.ui
    def b2b_view():
        payload = controller.b2b_html.get()
        if not payload:
            return ui.p("No Bio2Byte 3D rendering yet.")
        return ui.HTML(payload)

    @output(suspend_when_hidden=False)
    @render.ui
    def rin_view():
        if not controller.rin_html.get():
            return ui.p("No RIN built yet.")
        return ui.tags.iframe(
            srcdoc=controller.rin_html.get(),
            style="width:100%;height:700px;border:1px solid #ddd;border-radius:6px;",
        )

    @output(suspend_when_hidden=False)
    @render.text
    def tm_output() -> str:
        return controller.tm_report.get()

    @output(suspend_when_hidden=False)
    @render.ui
    def tm_actions():
        return gate(
            task_button("run_tmalign", "Align + Visualize", class_="btn-primary"),
            ready=controller.tm_structures_loaded.get(),
            hint="Load both structures first.",
        )

    @output(suspend_when_hidden=False)
    @render.ui
    def tm_view():
        payload = controller.tm_html.get()
        if not payload:
            return ui.p("No TM-align structure view yet.")
        return ui.HTML(payload)


content_ui = ui.div(
    # Spinners on recalculating outputs. Unlike a status message this reaches the
    # browser even while a synchronous handler runs.
    busy_indicators(),
    ui.tags.style(INERT_CSS),
    app_ui, scop3p_footer()
)

# static_assets serves the vendored browser libraries; every app mounts the same prefix,
# so /vendor/... resolves whichever app the portal is serving. enable_compression is not
# optional cosmetics: Shiny sends static files raw, and molstar.js is 5 MB uncompressed
# against 1.45 MB gzipped, so without it vendoring would put more bytes on the wire.
app = App(content_ui, server, static_assets=static_assets())
enable_compression(app)
