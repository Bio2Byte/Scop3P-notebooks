from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

import pandas as pd
from shiny import App, reactive, render, ui

from common.logging_utils import get_logger
from common.models import PeptideSelectionMode
from common.peptide_mapper import PeptideMapperService, map_selection
from common.services import AlphaFoldService, Scop3PClient
from common.ui_shell import scop3p_card, scop3p_shell, scop3p_footer
from common.viewer import NGLViewerBuilder


LOGGER = get_logger("scop3p.peptide_mapper")


class PeptideMapperController:
    """Stateful coordinator for Peptide Mapper app behavior."""

    def __init__(self) -> None:
        self.client = Scop3PClient()
        self.af_service = AlphaFoldService(cache_dir=Path("af_cache"))

        self.dataframe = reactive.value(pd.DataFrame())
        self.filtered_dataframe = reactive.value(pd.DataFrame())
        self.viewer_html = reactive.value("")
        self.summary_text = reactive.value("Load an accession to start.")
        self.status_text = reactive.value("")

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


controller = PeptideMapperController()


def _as_selectize_choices(options: list[tuple[str, str]]) -> dict[str, str]:
    # Shiny selectize choices must be value -> label.
    return {value: label for label, value in options}


app_ui = scop3p_shell(
    "Peptide Mapper",
    "Enter an accession, load peptides from Scop3P, filter and select mapped spans, then visualize coverage and modification sites on the AlphaFold structure.",
    ui.tags.style(
        """
        .pm-controls-card .btn {
          width: 100%;
          min-height: 54px;
          white-space: normal;
        }
        .pm-main-grid {
          display: grid;
          grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
          gap: 18px;
          align-items: start;
        }
        .pm-top-row {
          display: grid;
          grid-template-columns: minmax(180px, 1fr) 140px minmax(280px, 1.2fr);
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
          grid-template-columns: 1fr 1fr;
          gap: 12px;
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
            ui.div(
                ui.input_text("accession", "ACC_ID (UniProt accession number)", value="", placeholder="e.g. O00571"),
                ui.input_action_button("load_btn", "Load", class_="btn-primary"),
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
                class_="pm-top-row",
            ),
            ui.input_text(
                "search",
                "Search",
                placeholder="Filter: substring (SSFG), range (70-90), >=150, <=300, or single pos (154)",
            ),
            ui.input_selectize(
                "peptides",
                "Peptides",
                choices={},
                multiple=True,
                options={"placeholder": "Select peptide entries"},
            ),
            ui.div(
                ui.div(
                    ui.input_action_button("map_all", "Map all (filtered)", class_="btn-warning"),
                    ui.input_action_button("export_html", "Export styled HTML", class_="btn-info"),
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
    @reactive.effect
    @reactive.event(input.load_btn)
    def _load_data() -> None:
        accession = input.accession().strip()
        LOGGER.info("load requested accession=%s", accession or "-", extra={"event": "load_btn"})
        if not accession:
            controller.status_text.set("Enter an accession (e.g., O00571), then click Load.")
            return

        try:
            dataframe = controller.client.fetch_peptides_modifications(accession)
        except Exception as error:
            LOGGER.exception("load failed accession=%s", accession, extra={"event": "load_btn"})
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
            controller.status_text.set(f"No peptides returned for {accession}.")
            ui.update_selectize("peptides", choices={}, selected=[])
            return

        mode = PeptideSelectionMode(input.list_mode())
        options = PeptideMapperService.build_options(dataframe, mode)
        ui.update_selectize("peptides", choices=_as_selectize_choices(options), selected=[])
        controller.status_text.set(f"Loaded {len(dataframe)} peptide-mod rows for {accession}.")
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
        LOGGER.info("map_all requested", extra={"event": "map_all"})
        if filtered is None or filtered.empty:
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

        accession = input.accession().strip()
        LOGGER.info(
            "render requested accession=%s selected=%s show_mods=%s scope=%s",
            accession or "-",
            len(selected_values),
            input.show_mods(),
            input.mods_scope(),
            extra={"event": "render_selection"},
        )
        if not accession:
            controller.status_text.set("Enter an accession and click Load first.")
            return

        dataframe_all = controller.dataframe.get()
        dataframe_filtered = controller.filtered_dataframe.get()
        if dataframe_all is None or dataframe_all.empty:
            controller.status_text.set("No data loaded. Click Load first.")
            return

        try:
            pdb_path = controller.af_service.download_pdb(accession)
        except Exception as error:
            LOGGER.exception("alphafold download failed accession=%s", accession, extra={"event": "render_selection"})
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
                    f"ACC_ID: {accession}",
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
        accession = input.accession().strip()
        LOGGER.info("export requested accession=%s", accession or "-", extra={"event": "export_html"})
        if not accession:
            controller.status_text.set("Enter an accession first.")
            return

        if not controller.viewer_html.get() or controller.last_pdb_path.get() is None:
            controller.status_text.set("Render a selection first before exporting.")
            return

        export_path = Path("exports") / f"{accession}_styled_session.html"
        NGLViewerBuilder.export_html(export_path, controller.viewer_html.get())
        controller.status_text.set(f"Exported styled HTML to: {export_path.resolve()}")
        LOGGER.info("export completed path=%s", export_path.resolve(), extra={"event": "export_html"})

    @render.text
    def status() -> str:
        return controller.status_text.get()

    @render.text
    def summary() -> str:
        return controller.summary_text.get()

    @render.ui
    def viewer():
        payload = controller.viewer_html.get()
        if not payload:
            return ui.p("No structure rendered yet.")
        return ui.HTML(payload)


content_ui = ui.div(
    app_ui, scop3p_footer()
)

app = App(content_ui, server)
