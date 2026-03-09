from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

import pandas as pd
from shiny import App, reactive, render, ui

from apps.common.models import PeptideSelectionMode
from apps.common.peptide_mapper import PeptideMapperService, map_selection
from apps.common.services import AlphaFoldService, Scop3PClient
from apps.common.viewer import NGLViewerBuilder


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


app_ui = ui.page_fluid(
    ui.h2("Scop3P Peptide Mapper (Shiny)") ,
    ui.p(
        "Enter an accession, load peptides from Scop3P, filter/select peptides, "
        "and visualize mapped regions on AlphaFold structure."
    ),
    ui.layout_columns(
        ui.input_text("accession", "ACC_ID", value="", placeholder="e.g. O00571"),
        ui.input_action_button("load_btn", "Load", class_="btn-primary"),
        ui.input_radio_buttons(
            "list_mode",
            "List",
            choices=[PeptideSelectionMode.UNIQUE_SPANS.value, PeptideSelectionMode.ALL_ROWS.value],
            selected=PeptideSelectionMode.UNIQUE_SPANS.value,
            inline=True,
        ),
        col_widths=[3, 2, 7],
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
    ui.layout_columns(
        ui.input_action_button("map_all", "Map all (filtered)", class_="btn-warning"),
        ui.input_checkbox("show_mods", "Show modified sites (magenta)", value=True),
        ui.input_radio_buttons(
            "mods_scope",
            "Mods",
            choices=["Selected peptides only", "All protein mods"],
            selected="Selected peptides only",
            inline=True,
        ),
        ui.input_action_button("export_html", "Export styled HTML", class_="btn-info"),
        col_widths=[2, 3, 5, 2],
    ),
    ui.hr(),
    ui.h4("Status"),
    ui.output_text_verbatim("status"),
    ui.h4("Selection Summary"),
    ui.output_text_verbatim("summary"),
    ui.h4("Structure Viewer"),
    ui.output_ui("viewer"),
)


def server(input, output, session):
    @reactive.effect
    @reactive.event(input.load_btn)
    def _load_data() -> None:
        accession = input.accession().strip()
        if not accession:
            controller.status_text.set("Enter an accession (e.g., O00571), then click Load.")
            return

        try:
            dataframe = controller.client.fetch_peptides_modifications(accession)
        except Exception as error:
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

    @reactive.effect
    @reactive.event(input.map_all)
    def _map_all() -> None:
        filtered = controller.filtered_dataframe.get()
        if filtered is None or filtered.empty:
            controller.status_text.set("No filtered rows available. Load data first.")
            return

        mode = PeptideSelectionMode(input.list_mode())
        options = PeptideMapperService.build_options(filtered, mode)
        values = [value for _, value in options]
        ui.update_selectize("peptides", selected=values)

    @reactive.effect
    @reactive.event(input.peptides, input.show_mods, input.mods_scope)
    def _render_selection() -> None:
        selected_values = list(input.peptides())
        if not selected_values:
            return

        accession = input.accession().strip()
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

    @reactive.effect
    @reactive.event(input.export_html)
    def _export_html() -> None:
        accession = input.accession().strip()
        if not accession:
            controller.status_text.set("Enter an accession first.")
            return

        if not controller.viewer_html.get() or controller.last_pdb_path.get() is None:
            controller.status_text.set("Render a selection first before exporting.")
            return

        export_path = Path("exports") / f"{accession}_styled_session.html"
        NGLViewerBuilder.export_html(export_path, controller.viewer_html.get())
        controller.status_text.set(f"Exported styled HTML to: {export_path.resolve()}")

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


app = App(app_ui, server)
