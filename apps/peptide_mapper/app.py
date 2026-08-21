from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

import pandas as pd
from shiny import App, reactive, render, ui

from common.busy import busy_indicators, task_button
from common.logging_utils import get_logger, new_trail
from common.models import PeptideSelectionMode
from common.peptide_mapper import (
    PeptideColumnMapping,
    PeptideMapperService,
    build_upload_mapping,
    detect_peptide_columns,
    map_selection,
    mapped_residue_rows,
    peptides_for_protein,
    protein_choices,
    read_peptide_table,
)
from common.services import AlphaFoldService, Scop3PClient
from common.ui_shell import (
    ACCESSION_LABEL,
    scop3p_card,
    scop3p_example_button,
    scop3p_field_row,
    scop3p_footer,
    scop3p_shell,
)
from common.viewer import NGLViewerBuilder


LOGGER = get_logger("scop3p.peptide_mapper")

#: Title of the upload source tab. Compared against input.source_tabs() to decide
#: which source is active; a navset with an id reports the ACTIVE PANEL TITLE.
UPLOAD_TAB = "Upload your own"

#: Worked example: DDX3X, a well-covered phosphoprotein in Scop3P.
EXAMPLE_ACCESSION = "O00571"


# Stateless and safe to share: the Scop3P client holds no request state, and the
# AlphaFold service is a disk cache we *want* shared across sessions.
_CLIENT = Scop3PClient()
_AF_SERVICE = AlphaFoldService(cache_dir=Path("af_cache"))


class PeptideMapperController:
    """Stateful coordinator for Peptide Mapper app behavior.

    Instantiated once per Shiny session, inside ``server()``. Do not hoist this
    back to module level: every attribute below is a ``reactive.value``, so a
    module-level instance is shared by every connected browser and one user's
    peptide table renders in another user's viewer.
    """

    def __init__(self) -> None:
        self.client = _CLIENT
        self.af_service = _AF_SERVICE

        self.dataframe = reactive.value(pd.DataFrame())
        self.filtered_dataframe = reactive.value(pd.DataFrame())
        self.viewer_html = reactive.value("")
        self.summary_text = reactive.value("Load an accession to start.")
        self.status_text = reactive.value("")

        # Upload mode. `dataframe` above always holds the peptide table currently in
        # play, whatever its source, so everything downstream is source-agnostic.
        self.upload_raw = reactive.value(pd.DataFrame())
        self.upload_mapped = reactive.value(pd.DataFrame())
        self.upload_name = reactive.value("")

        self.last_pdb_path = reactive.value(None)
        self.last_union_ranges = reactive.value([])
        self.last_intersection_positions = reactive.value([])
        self.last_modification_positions = reactive.value([])

    def clear_render_state(self) -> None:
        self.viewer_html.set("")
        self.summary_text.set("Load data and select peptides to render.")
        self.last_pdb_path.set(None)
        self.last_union_ranges.set([])
        self.last_intersection_positions.set([])
        self.last_modification_positions.set([])


def _as_selectize_choices(options: list[tuple[str, str]]) -> dict[str, str]:
    # Shiny selectize choices must be value -> label.
    return {value: label for label, value in options}


app_ui = scop3p_shell(
    "Peptide Mapper",
    "Enter an accession, load peptides from Scop3P, filter and select mapped spans, then visualize coverage and modification sites on the AlphaFold structure.",
    ui.tags.style(
        """
        /* Scoped to the action grid. When this matched every .btn in the card it also
           hit the buttons inside .scop3p-field-row, forcing them full width and
           breaking the input/button baseline. */
        .pm-button-grid .btn {
          width: 100%;
          min-height: 54px;
          white-space: normal;
        }
        .pm-main-grid {
          display: grid;
          grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
          gap: 18px;
          align-items: stretch;
        }
        .pm-main-grid > .scop3p-card {
          height: 100%;
        }
        .pm-top-row {
          display: grid;
          grid-template-columns: minmax(180px, 1fr) 140px;
          gap: 14px;
          align-items: end;
          margin-bottom: 14px;
        }
        .pm-actions-row {
          display: grid;
          grid-template-columns: minmax(260px, 320px) minmax(280px, 1fr);
          gap: 18px;
          align-items: start;
        }
        .pm-button-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
          gap: 12px;
        }
        /* The peptide selector is a dropdown with multi-select, not a list. Cap the
           control so a large selection scrolls inside it rather than growing the card,
           and clip long labels so the card never scrolls sideways. */
        .pm-peptides .selectize-input {
          max-height: 150px;
          overflow-y: auto;
        }
        /* The chip must always show its remove control, however long the label.
           Selectize renders the "x" as the last child *inside* .item, and its own CSS
           makes .item an inline-flex box. The label is then an anonymous flex item, which
           cannot shrink below its content width (the flexbox min-width:auto trap), so the
           "x" was pushed 215px past the clip edge and disappeared -- the chip could be
           selected but never removed.

           Fixed by dropping back to inline-block, so text-overflow can ellipsize the
           label, and reserving room on the right for the "x" via padding with the link
           absolutely positioned into it. overflow:hidden clips at the padding edge, so
           anything inside that reserved strip stays visible. */
        .pm-peptides .selectize-input > .item {
          position: relative;
          display: inline-block;
          max-width: 100%;
          padding-right: 22px;
          overflow: hidden;
          text-overflow: ellipsis;
          white-space: nowrap;
          vertical-align: top;
        }
        .pm-peptides .selectize-input > .item > .remove {
          position: absolute;
          top: 0;
          right: 0;
          width: 20px;
          text-align: center;
          border-left: 1px solid rgba(0, 0, 0, 0.12);
        }
        .pm-peptides .selectize-dropdown-content {
          max-height: 260px;
        }
        .pm-selection-count {
          margin: 6px 0 0;
          font-size: 0.85rem;
          color: var(--scop3p-muted);
        }
        .pm-upload-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
          gap: 10px;
          margin: 10px 0;
        }
        .pm-mods-row .form-group,
        .pm-list-block .form-group {
          margin-bottom: 0;
        }
        @media (max-width: 1200px) {
          .pm-main-grid,
          .pm-top-row,
          .pm-actions-row,
          .pm-button-grid {
            grid-template-columns: 1fr;
          }
        }
        """
    ),
    ui.div(
        scop3p_card(
            "Session Status",
            ui.output_text_verbatim("status"),
            extra_class="scop3p-status",
        ),
        scop3p_card(
            "Selection Summary",
            ui.output_text_verbatim("summary"),
            extra_class="scop3p-status",
        ),
        class_="scop3p-header-grid",
    ),
    ui.div(
        scop3p_card(
            "Controls",
            # These tabs cover the SOURCE of the peptide table only. Search,
            # selection, mapping, the viewer and the exports live outside the navset
            # and are shared, so nothing below is duplicated per source.
            ui.navset_tab(
                ui.nav_panel(
                    "Scop3P peptides",
                    scop3p_field_row(
                        ui.input_text(
                            "accession",
                            ACCESSION_LABEL,
                            value="",
                            placeholder=f"e.g. {EXAMPLE_ACCESSION}",
                        ),
                        task_button(
                        "load_btn", "Load", class_="btn btn-primary"),
                        scop3p_example_button("load_example"),
                    ),
                ),
                ui.nav_panel(
                    "Upload your own",
                    ui.input_file(
                        "peptide_upload",
                        "Peptide table (TSV/CSV)",
                        accept=[".tsv", ".txt", ".csv"],
                        multiple=False,
                    ),
                    task_button(
                        "load_upload", "Load file", class_="btn-primary"),
                    ui.div(
                        ui.input_select("col_protein", "Protein ID", choices={}),
                        ui.input_select("col_sequence", "Peptide sequence", choices={}),
                        ui.input_select("col_start", "Peptide start", choices={}),
                        ui.input_select("col_end", "Peptide end", choices={}),
                        ui.input_select("col_position", "UniProt position", choices={}),
                        class_="pm-upload-grid",
                    ),
                    task_button(
                        "build_mapping", "Build mapping", class_="btn-success"),
                    ui.input_select("protein", "Protein", choices={}),
                    ui.output_ui("upload_preview"),
                    ui.p(
                        "Columns are detected automatically; correct any wrong guess, "
                        "then build the mapping.",
                        class_="scop3p-note",
                    ),
                ),
                id="source_tabs",
            ),
            ui.div(
                ui.input_radio_buttons(
                    "list_mode",
                    "List",
                    choices=[PeptideSelectionMode.UNIQUE_SPANS.value, PeptideSelectionMode.ALL_ROWS.value],
                    selected=PeptideSelectionMode.UNIQUE_SPANS.value,
                    inline=False,
                ),
                class_="pm-list-block",
            ),
            ui.input_text(
                "search",
                "Search",
                placeholder="Filter: substring (SSFG), range (70-90), >=150, <=300, or single pos (154)",
            ),
            # Wrapped so the CSS can cap the control's height. "Map all (filtered)"
            # can select 46 spans at once, and selectize renders every selection as a
            # chip on its own line, so the control grew past the height of the card
            # and pushed the viewer off screen.
            ui.div(
                ui.input_selectize(
                    "peptides",
                    "Peptides",
                    choices={},
                    multiple=True,
                    options={"placeholder": "Select peptide entries"},
                ),
                ui.output_ui("peptide_selection_count"),
                class_="pm-peptides",
            ),
            ui.div(
                ui.div(
                    task_button(
                        "map_all", "Map all (filtered)", class_="btn-warning"),
                    task_button(
                        "export_html", "Export styled HTML", class_="btn-info"),
                    ui.download_button("download_html", "Download HTML", class_="btn-info"),
                    ui.download_button("download_pdb", "Download PDB", class_="btn-secondary"),
                    ui.download_button("download_residues", "Residues (TSV)", class_="btn-secondary"),
                    class_="pm-button-grid",
                ),
                ui.div(
                    ui.input_checkbox("show_mods", "Show modified sites (magenta)", value=True),
                    ui.div(
                        ui.input_radio_buttons(
                            "mods_scope",
                            "Mods",
                            choices=["Selected peptides only", "All protein mods"],
                            selected="Selected peptides only",
                            inline=False,
                        ),
                        class_="pm-mods-row",
                    ),
                ),
                class_="pm-actions-row",
            ),
            extra_class="pm-controls-card",
        ),
        scop3p_card(
            "Structure Viewer",
            ui.output_ui("viewer"),
        ),
        class_="pm-main-grid",
    ),
)


def server(input, output, session):
    # One trail per browser session: step numbers must not interleave across
    # sessions, and a module-level trail would be shared by every user.
    trail = new_trail()
    trail.opened("Peptide Mapper")

    controller = PeptideMapperController()

    def active_accession() -> str:
        """The accession to fetch a structure for, per the selected source tab.

        Branch on the tab, not on "whichever field is non-empty": Shiny keeps the
        inputs of an inactive tab readable, so a stale accession left in the Scop3P
        field would otherwise silently override the uploaded protein.
        """
        if input.source_tabs() == UPLOAD_TAB:
            return (input.protein() or "").strip()
        return input.accession().strip()

    def _column_mapping() -> PeptideColumnMapping:
        return PeptideColumnMapping(
            protein=input.col_protein() or None,
            sequence=input.col_sequence() or None,
            start=input.col_start() or None,
            end=input.col_end() or None,
            position=input.col_position() or None,
        )

    @reactive.effect
    @reactive.event(input.load_example)
    def _load_example() -> None:
        trail.clicked("Load example")
        ui.update_text("accession", value=EXAMPLE_ACCESSION)
        controller.status_text.set(
            f"Example accession {EXAMPLE_ACCESSION} loaded. Click Load."
        )

    @reactive.effect
    @reactive.event(input.load_upload)
    def _load_upload() -> None:
        trail.clicked("Load file")
        uploaded = input.peptide_upload()
        if not uploaded:
            controller.status_text.set("Choose a peptide table (TSV or CSV) first.")
            return
        item = uploaded[0]
        name = item.get("name") or "uploaded"
        try:
            raw, delimiter = read_peptide_table(item["datapath"])
        except Exception as error:
            LOGGER.exception("upload read failed name=%s", name, extra={"event": "load_upload"})
            trail.failed("upload read failed", error=type(error).__name__)
            controller.status_text.set(f"Could not read {name}: {error}")
            return

        controller.upload_raw.set(raw)
        controller.upload_name.set(name)

        columns = [str(column) for column in raw.columns]
        choices = {column: column for column in columns}
        detected = detect_peptide_columns(columns)
        for input_id, field in (
            ("col_protein", "protein"),
            ("col_sequence", "sequence"),
            ("col_start", "start"),
            ("col_end", "end"),
            ("col_position", "position"),
        ):
            ui.update_select(input_id, choices=choices, selected=getattr(detected, field))

        message = (
            f"Loaded {name} ({delimiter}-separated): {len(raw)} rows x {len(columns)} columns."
        )
        if detected.missing:
            message += (
                " Could not identify: " + ", ".join(detected.missing)
                + ". Pick those columns by hand, then click Build mapping."
            )
        else:
            message += " Columns detected. Click Build mapping."
        controller.status_text.set(message)
        trail.produced(
            f"uploaded table read: {name}",
            rows=len(raw), columns=len(columns), unmapped=",".join(detected.missing) or "-",
        )
        LOGGER.info(
            "load_upload completed name=%s rows=%s columns=%s missing=%s",
            name,
            len(raw),
            len(columns),
            ",".join(detected.missing) or "-",
            extra={"event": "load_upload"},
        )

    @reactive.effect
    @reactive.event(input.build_mapping)
    def _build_mapping() -> None:
        trail.clicked("Build mapping")
        raw = controller.upload_raw.get()
        if raw is None or raw.empty:
            controller.status_text.set("Load a peptide table first.")
            return
        try:
            mapped = build_upload_mapping(raw, _column_mapping())
        except Exception as error:
            LOGGER.exception("build_mapping failed", extra={"event": "build_mapping"})
            trail.failed("build_mapping failed", error=type(error).__name__)
            controller.status_text.set(f"Mapping error: {error}")
            return

        if mapped.empty:
            controller.status_text.set(
                "No usable rows after mapping. Check that start, end and position are "
                "numeric and that start is 1-indexed."
            )
            return

        choices = protein_choices(mapped)
        first = next(iter(choices))
        # `dataframe` is the single source everything downstream reads, so the upload
        # path ends here: one protein's rows go in and the Scop3P pipeline continues.
        controller.upload_mapped.set(mapped)
        controller.dataframe.set(peptides_for_protein(mapped, first))
        controller.filtered_dataframe.set(controller.dataframe.get())
        controller.clear_render_state()
        ui.update_select("protein", choices=choices, selected=first)

        controller.status_text.set(
            f"Mapped {len(mapped)} rows across {len(choices)} protein(s). "
            f"Showing {first}."
        )
        LOGGER.info(
            "build_mapping completed rows=%s proteins=%s",
            len(mapped),
            len(choices),
            extra={"event": "build_mapping"},
        )

    @reactive.effect
    @reactive.event(input.protein)
    def _protein_changed() -> None:
        mapped = controller.upload_mapped.get()
        accession = (input.protein() or "").strip()
        if mapped is None or mapped.empty or not accession:
            return
        subset = peptides_for_protein(mapped, accession)
        controller.dataframe.set(subset)
        controller.filtered_dataframe.set(subset)
        controller.clear_render_state()
        ui.update_selectize("peptides", choices={}, selected=[])
        controller.status_text.set(f"Showing {len(subset)} peptide rows for {accession}.")

    @reactive.effect
    @reactive.event(input.load_btn)
    def _load_data() -> None:
        accession = input.accession().strip()
        trail.entered(ACCESSION_LABEL, accession or "-")
        trail.clicked("Load")
        LOGGER.info("load requested accession=%s", accession or "-", extra={"event": "load_btn"})
        if not accession:
            trail.blocked("missing accession")
            controller.status_text.set("Enter an accession (e.g., O00571), then click Load.")
            return

        try:
            dataframe = controller.client.fetch_peptides_modifications(accession)
        except Exception as error:
            LOGGER.exception("load failed accession=%s", accession, extra={"event": "load_btn"})
            trail.failed("load failed", error=type(error).__name__)
            controller.status_text.set(f"Scop3P API error: {error}")
            controller.dataframe.set(pd.DataFrame())
            controller.filtered_dataframe.set(pd.DataFrame())
            controller.clear_render_state()
            ui.update_selectize("peptides", choices={}, selected=[])
            return

        controller.dataframe.set(dataframe)
        controller.filtered_dataframe.set(dataframe)
        controller.clear_render_state()

        if dataframe.empty:
            LOGGER.warning("load completed without rows accession=%s", accession, extra={"event": "load_btn"})
            controller.status_text.set(f"No peptides returned for {accession}.")
            ui.update_selectize("peptides", choices={}, selected=[])
            return

        mode = PeptideSelectionMode(input.list_mode())
        options = PeptideMapperService.build_options(dataframe, mode)
        ui.update_selectize("peptides", choices=_as_selectize_choices(options), selected=[])
        controller.status_text.set(f"Loaded {len(dataframe)} peptide-mod rows for {accession}.")
        trail.produced(f"{len(dataframe)} peptides loaded", mode=mode.value)
        LOGGER.info(
            "load completed accession=%s rows=%s mode=%s",
            accession,
            len(dataframe),
            mode.value,
            extra={"event": "load_btn"},
        )

    @reactive.effect
    def _update_filter_and_choices() -> None:
        dataframe = controller.dataframe.get()
        if dataframe is None or dataframe.empty:
            return

        filtered = PeptideMapperService.filter_peptides(dataframe, input.search())
        controller.filtered_dataframe.set(filtered)

        mode = PeptideSelectionMode(input.list_mode())
        options = PeptideMapperService.build_options(filtered, mode)

        selected = list(input.peptides())
        valid_values = {value for _, value in options}
        restored = [value for value in selected if value in valid_values]

        ui.update_selectize("peptides", choices=_as_selectize_choices(options), selected=restored)
        LOGGER.info(
            "filter updated query=%r filtered_rows=%s options=%s restored=%s",
            input.search(),
            len(filtered),
            len(options),
            len(restored),
            extra={"event": "filter_update"},
        )

    @reactive.effect
    @reactive.event(input.map_all)
    def _map_all() -> None:
        filtered = controller.filtered_dataframe.get()
        trail.clicked("Map all (filtered)")
        LOGGER.info("map_all requested", extra={"event": "map_all"})
        if filtered is None or filtered.empty:
            trail.blocked("no filtered rows")
            controller.status_text.set("No filtered rows available. Load data first.")
            return

        mode = PeptideSelectionMode(input.list_mode())
        options = PeptideMapperService.build_options(filtered, mode)
        values = [value for _, value in options]
        ui.update_selectize("peptides", selected=values)
        LOGGER.info("map_all selected_count=%s", len(values), extra={"event": "map_all"})

    @reactive.effect
    @reactive.event(input.peptides, input.show_mods, input.mods_scope)
    def _render_selection() -> None:
        selected_values = list(input.peptides())
        if not selected_values:
            return

        accession = active_accession()
        LOGGER.info(
            "render requested accession=%s selected=%s show_mods=%s scope=%s",
            accession or "-",
            len(selected_values),
            input.show_mods(),
            input.mods_scope(),
            extra={"event": "render_selection"},
        )
        if not accession:
            trail.blocked("missing accession")
            controller.status_text.set(
                "Choose a protein in the Upload tab first."
                if input.source_tabs() == UPLOAD_TAB
                else "Enter an accession and click Load first."
            )
            return

        dataframe_all = controller.dataframe.get()
        dataframe_filtered = controller.filtered_dataframe.get()
        if dataframe_all is None or dataframe_all.empty:
            trail.blocked("no data loaded")
            controller.status_text.set("No data loaded. Click Load first.")
            return

        try:
            pdb_path = controller.af_service.download_pdb(accession)
        except Exception as error:
            LOGGER.exception("alphafold download failed accession=%s", accession, extra={"event": "render_selection"})
            trail.failed("AlphaFold download failed", error=type(error).__name__)
            controller.status_text.set(f"AlphaFold download error: {error}")
            return

        try:
            union_ranges, intersection_positions, modification_positions = map_selection(
                dataframe_all=dataframe_all,
                dataframe_filtered=dataframe_filtered,
                selected_keys=selected_values,
                mode=PeptideSelectionMode(input.list_mode()),
                mods_scope=input.mods_scope(),
            )
        except Exception as error:
            LOGGER.exception("selection mapping failed accession=%s", accession, extra={"event": "render_selection"})
            trail.failed("selection mapping failed", error=type(error).__name__)
            controller.status_text.set(f"Selection mapping error: {error}")
            return

        if not input.show_mods():
            modification_positions = []

        html_payload = NGLViewerBuilder.build_html(
            accession=accession,
            pdb_path=pdb_path,
            union_ranges=union_ranges,
            intersection_positions=intersection_positions,
            modification_positions=modification_positions,
        )

        controller.viewer_html.set(html_payload)
        controller.last_pdb_path.set(pdb_path)
        controller.last_union_ranges.set(union_ranges)
        controller.last_intersection_positions.set(intersection_positions)
        controller.last_modification_positions.set(modification_positions)

        coverage_start = min(start for start, _ in union_ranges) if union_ranges else "-"
        coverage_end = max(end for _, end in union_ranges) if union_ranges else "-"
        controller.summary_text.set(
            "\n".join(
                [
                    f"Accession: {accession}",
                    f"AlphaFold model: {pdb_path}",
                    f"Selected entries: {len(selected_values)}",
                    f"Coverage: {coverage_start} -> {coverage_end}",
                    f"Intersection (red): {len(intersection_positions)} residues",
                    f"Modified sites (magenta): {len(set(modification_positions))} unique positions",
                ]
            )
        )
        LOGGER.info(
            "render completed accession=%s pdb=%s ranges=%s intersection=%s mods=%s",
            accession,
            pdb_path,
            len(union_ranges),
            len(intersection_positions),
            len(set(modification_positions)),
            extra={"event": "render_selection"},
        )

    @reactive.effect
    @reactive.event(input.export_html)
    def _export_html() -> None:
        accession = active_accession()
        trail.clicked("Export styled HTML")
        LOGGER.info("export requested accession=%s", accession or "-", extra={"event": "export_html"})
        if not accession:
            trail.blocked("missing accession")
            controller.status_text.set("Enter an accession first.")
            return

        if not controller.viewer_html.get() or controller.last_pdb_path.get() is None:
            trail.blocked("no rendered selection")
            controller.status_text.set("Render a selection first before exporting.")
            return

        export_path = Path("exports") / f"{accession}_styled_session.html"
        NGLViewerBuilder.export_html(export_path, controller.viewer_html.get())
        controller.status_text.set(f"Exported styled HTML to: {export_path.resolve()}")
        trail.exported(f"styled HTML session: {export_path.name}")
        LOGGER.info("export completed path=%s", export_path.resolve(), extra={"event": "export_html"})

    @render.text
    def status() -> str:
        return controller.status_text.get()

    @render.text
    def summary() -> str:
        return controller.summary_text.get()

    @render.ui
    def peptide_selection_count():
        """How many peptides are selected, since the chips can now be scrolled away."""
        selected = list(input.peptides())
        if not selected:
            return None
        # `df or []` raises: a DataFrame has no truth value. Check for None explicitly.
        filtered = controller.filtered_dataframe.get()
        total = 0 if filtered is None else len(filtered)
        noun = "peptide" if len(selected) == 1 else "peptides"
        return ui.p(
            f"{len(selected)} {noun} selected"
            + (f" of {total} filtered rows" if total else ""),
            class_="pm-selection-count",
        )

    @render.ui
    def viewer():
        payload = controller.viewer_html.get()
        if not payload:
            return ui.p("No structure rendered yet.")
        return ui.HTML(payload)

    # suspend_when_hidden=False because this output sits in the "Upload your own"
    # tab, which is not the initially-active one; see the note in
    # apps/rinalign/app.py for why such outputs never wake up otherwise.
    @output(suspend_when_hidden=False)
    @render.ui
    def upload_preview():
        raw = controller.upload_raw.get()
        if raw is None or raw.empty:
            return ui.p(
                "Upload a TSV or CSV exported from your search engine.",
                class_="scop3p-note",
            )
        preview = raw.head(5)
        return ui.TagList(
            ui.p(
                f"{controller.upload_name.get()} - first {len(preview)} of {len(raw)} rows",
                class_="scop3p-note",
            ),
            ui.HTML(
                "<div style='overflow-x:auto;max-width:100%;'>"
                + preview.to_html(index=False, border=0, classes="table table-sm")
                + "</div>"
            ),
        )

    # Downloads rather than server-side files. The notebook wrote into exports/
    # because Voila had no download primitive; on a shared server that path is a
    # cross-session collision and leaves user data on the host.
    @render.download_button(filename=lambda: f"{active_accession() or 'peptides'}_session.html")
    def download_html():
        payload = controller.viewer_html.get()
        if not payload:
            yield "<!doctype html><p>Render a selection before downloading.</p>"
            return
        yield payload

    @render.download_button(filename=lambda: f"{active_accession() or 'structure'}.pdb")
    def download_pdb():
        path = controller.last_pdb_path.get()
        if not path:
            yield "REMARK  Render a selection first; no structure has been fetched.\n"
            return
        yield Path(path).read_text(encoding="utf-8", errors="replace")

    @render.download_button(filename=lambda: f"{active_accession() or 'peptides'}_mapped_residues.tsv")
    def download_residues():
        table = mapped_residue_rows(
            controller.last_union_ranges.get(),
            controller.last_intersection_positions.get(),
            controller.last_modification_positions.get(),
        )
        yield table.to_csv(sep="\t", index=False)


content_ui = ui.div(
    busy_indicators(),
    app_ui, scop3p_footer()
)

app = App(content_ui, server)
