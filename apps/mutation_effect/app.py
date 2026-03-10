from __future__ import annotations

import sys
from pathlib import Path

from shiny import App, reactive, render, ui

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from apps.common.mutation_effect import (  # noqa: E402
    MutationEffectInference,
    MutationEffectService,
    MutationEffectViews,
)


service = MutationEffectService()


def _bokeh_iframe(html_doc: str) -> ui.Tag:
    return ui.tags.iframe(
        srcdoc=html_doc,
        style="width:100%;height:380px;border:0;border-radius:10px;background:#fff;",
    )


def _panel(title: str, *children) -> ui.Tag:
    return ui.div(
        ui.h4(title, style="margin-top:0;margin-bottom:0.75rem;"),
        *children,
        style="background:#fff;border:1px solid #d9dee8;border-radius:14px;padding:1rem;",
    )


app_ui = ui.page_fluid(
    ui.tags.style(
        """
        :root {
          --bg: #f6f4ef;
          --panel: #ffffff;
          --ink: #1f2937;
          --muted: #617184;
          --accent: #1f6fb2;
          --accent-2: #d2872c;
          --line: #d9dee8;
        }
        body {
          background:
            radial-gradient(circle at top right, rgba(210,135,44,0.08), transparent 26%),
            linear-gradient(180deg, #f8f6f1 0%, #eef2f7 100%);
          color: var(--ink);
        }
        .app-shell {
          max-width: 1600px;
          margin: 0 auto;
          padding: 18px 10px 32px;
        }
        .hero {
          display: grid;
          grid-template-columns: 1.2fr 0.8fr;
          gap: 18px;
          margin-bottom: 18px;
        }
        .hero-copy {
          background: linear-gradient(135deg, #11263a 0%, #234d73 100%);
          color: #fff;
          border-radius: 18px;
          padding: 22px 24px;
          box-shadow: 0 18px 40px rgba(17,38,58,0.18);
        }
        .hero-copy h2 {
          margin: 0 0 8px;
          font-size: 2rem;
          font-weight: 700;
          letter-spacing: -0.02em;
        }
        .hero-copy p {
          margin: 0;
          max-width: 70ch;
          color: rgba(255,255,255,0.84);
        }
        .status-card {
          background: rgba(255,255,255,0.9);
          border: 1px solid rgba(255,255,255,0.65);
          border-radius: 18px;
          padding: 18px;
          box-shadow: 0 18px 40px rgba(17,38,58,0.08);
        }
        .status-card pre {
          margin: 0;
          white-space: pre-wrap;
          font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
          font-size: 13px;
        }
        .section-grid {
          display: grid;
          grid-template-columns: 420px minmax(0, 1fr);
          gap: 18px;
          align-items: start;
        }
        .note {
          color: var(--muted);
          font-size: 0.95rem;
        }
        .shiny-input-container {
          width: 100%;
        }
        .nav-tabs {
          margin-top: 8px;
        }
        .nav-tabs .nav-link {
          color: #274560;
          font-weight: 600;
        }
        .nav-tabs .nav-link.active {
          background: #fff;
          border-color: var(--line) var(--line) #fff;
        }
        .action-row {
          display: flex;
          gap: 10px;
          align-items: end;
        }
        .summary-banner {
          margin-bottom: 12px;
          padding: 10px 12px;
          border-radius: 10px;
          background: #fbf6e8;
          border: 1px solid #ecd7ab;
          color: #684a17;
          font-weight: 600;
        }
        @media (max-width: 1000px) {
          .hero, .section-grid {
            grid-template-columns: 1fr;
          }
        }
        """
    ),
    ui.div(
        ui.div(
            ui.div(
                ui.h2("Mutation Effect", style="margin:0;"),
                ui.p(
                    "Compare WT and mutant Bio2Byte predictions, overlay Scop3P PTMs, and derive mutation-centric inference summaries.",
                    style="margin:8px 0 0;",
                ),
                class_="hero-copy",
            ),
            ui.div(ui.output_text_verbatim("status"), class_="status-card"),
            class_="hero",
        ),
        ui.navset_tab(
            ui.nav_panel(
                "1) WT prediction",
                ui.div(
                    _panel(
                        "Protein Setup",
                        ui.input_text("accession", "UniProt accession", value="P07949"),
                        ui.input_action_button("run_wt", "Fetch + Predict WT", class_="btn btn-primary"),
                        ui.p("This runs the wild-type UniProt fetch, Scop3P PTM fetch, and Bio2Byte prediction.", class_="note"),
                    ),
                    _panel(
                        "WT Results",
                        ui.output_ui("wt_results"),
                    ),
                    class_="section-grid",
                ),
            ),
            ui.nav_panel(
                "2) Mutant prediction",
                ui.div(
                    _panel(
                        "Mutation Setup",
                        ui.input_text("positions", "Positions", value="606", placeholder="e.g. 10,25,100"),
                        ui.input_text("mut_aas", "To amino acid", value="A", placeholder="e.g. A,V,G"),
                        ui.input_action_button("run_mut", "Apply + Predict", class_="btn btn-warning"),
                        ui.p("Use comma-separated 1-indexed positions and amino-acid targets.", class_="note"),
                    ),
                    _panel(
                        "Mutant Results",
                        ui.output_ui("mut_results"),
                    ),
                    class_="section-grid",
                ),
            ),
            ui.nav_panel(
                "3) Inference",
                ui.div(
                    _panel(
                        "Label Shift Analysis",
                        ui.input_action_button("run_inf", "Run inference", class_="btn btn-info"),
                        ui.p("Summarize class shifts at the mutation site and within a +/-5 residue window.", class_="note"),
                    ),
                    _panel(
                        "Inference Results",
                        ui.output_ui("inf_results"),
                    ),
                    class_="section-grid",
                ),
            ),
        ),
        class_="app-shell",
    ),
    title="Mutation Effect (Shiny)",
)


def server(input, output, session):
    status_text = reactive.Value("Ready.")
    wt_plot_html = reactive.Value("")
    wt_table_html = reactive.Value("")
    mut_summary_html = reactive.Value("")
    mut_plot_html = reactive.Value("")
    mut_table_html = reactive.Value("")
    inf_sections = reactive.Value([])

    accession = reactive.Value("")
    sequence = reactive.Value("")
    mods_df = reactive.Value(None)
    wt_df = reactive.Value(None)
    mut_df = reactive.Value(None)
    mutations = reactive.Value([])

    @output
    @render.text
    def status() -> str:
        return status_text.get()

    @reactive.effect
    @reactive.event(input.run_wt)
    def _run_wt() -> None:
        try:
            accession_value = input.accession().strip()
            status_text.set("Fetching UniProt sequence and Scop3P PTMs...")
            sequence_value = service.fetch_uniprot_sequence(accession_value)
            mods = service.fetch_scop3p_modifications(accession_value)

            status_text.set("Running WT biophysical prediction...")
            prediction = service.predict_biophysical(accession_value, sequence_value)
            dataframe = service.prediction_to_df(prediction, accession_value)
            dataframe["seq"] = list(sequence_value)

            accession.set(accession_value)
            sequence.set(sequence_value)
            mods_df.set(mods)
            wt_df.set(dataframe)
            mut_df.set(None)
            mutations.set([])
            mut_summary_html.set("")
            mut_plot_html.set("")
            mut_table_html.set("")
            inf_sections.set([])

            wt_plot_html.set(MutationEffectViews.make_wt_plot(dataframe, mods))
            wt_table_html.set(
                MutationEffectViews.scrollable_table_html(
                    MutationEffectViews.make_wt_table(dataframe, mods),
                    title="WT predicted features (per residue)",
                    height_px=420,
                    width="100%",
                    sticky_cols=0,
                    col_widths_px=[90, 70, 60],
                )
            )
            status_text.set(f"WT prediction ready. PTMs: {0 if mods is None else len(mods)}.")
        except Exception as error:
            status_text.set(f"WT error: {error}")

    @reactive.effect
    @reactive.event(input.run_mut)
    def _run_mut() -> None:
        try:
            if wt_df.get() is None or not sequence.get():
                raise ValueError("Run WT prediction first.")

            parsed_mutations = service.parse_mutations(input.positions(), input.mut_aas())
            status_text.set("Applying mutations and predicting mutant...")

            mutant_sequence = service.apply_mutations(sequence.get(), parsed_mutations)
            prediction = service.predict_biophysical(accession.get(), mutant_sequence)
            dataframe = service.prediction_to_df(prediction, accession.get())
            dataframe["seq"] = list(mutant_sequence)

            mut_df.set(dataframe)
            mutations.set(parsed_mutations)

            labels = service.build_mutation_labels(sequence.get(), parsed_mutations)
            mut_summary_html.set(
                f"<div class='summary-banner'>Predicted properties for mutants: {', '.join(labels) if labels else '(none)'}</div>"
            )
            mut_plot_html.set(MutationEffectViews.make_mut_plot(wt_df.get(), dataframe, mods_df.get()))
            mut_table_html.set(
                MutationEffectViews.scrollable_table_html(
                    MutationEffectViews.make_wt_mut_merged_table(wt_df.get(), dataframe, mods_df.get()),
                    title="WT vs Mutant predicted features (aligned by seqpos)",
                    height_px=420,
                    width="100%",
                    sticky_cols=0,
                    highlight_seqpos=[mutation.position for mutation in parsed_mutations],
                    seqpos_col="seqpos",
                )
            )
            inf_sections.set([])
            status_text.set("Mutant prediction ready.")
        except Exception as error:
            status_text.set(f"Mutant error: {error}")

    @reactive.effect
    @reactive.event(input.run_inf)
    def _run_inf() -> None:
        try:
            if wt_df.get() is None or mut_df.get() is None:
                raise ValueError("Run WT prediction and Mutant prediction first.")

            status_text.set("Running inference...")
            sections = []
            for feature, (_, title) in MutationEffectInference.LABEL_FUNCS.items():
                dataframe = MutationEffectInference.mutation_effect_table_with_label_shift(
                    wt_df=wt_df.get(),
                    mut_df=mut_df.get(),
                    feature=feature,
                    mutations=mutations.get(),
                    window=5,
                )
                sections.append(
                    (
                        title,
                        MutationEffectViews.scrollable_table_html(
                            dataframe,
                            title=title,
                            height_px=320,
                            width="100%",
                            sticky_cols=0,
                        ),
                    )
                )
            inf_sections.set(sections)
            status_text.set("Inference ready.")
        except Exception as error:
            status_text.set(f"Inference error: {error}")

    @output
    @render.ui
    def wt_results():
        if not wt_plot_html.get():
            return ui.p("No WT prediction yet.", class_="note")
        return ui.TagList(
            ui.HTML(MutationEffectViews.track_guide_html()),
            _bokeh_iframe(wt_plot_html.get()),
            ui.hr(),
            ui.HTML(wt_table_html.get()),
        )

    @output
    @render.ui
    def mut_results():
        if not mut_plot_html.get():
            return ui.p("No mutant prediction yet.", class_="note")
        return ui.TagList(
            ui.HTML(mut_summary_html.get()),
            ui.HTML(MutationEffectViews.track_guide_html()),
            _bokeh_iframe(mut_plot_html.get()),
            ui.hr(),
            ui.HTML(mut_table_html.get()),
        )

    @output
    @render.ui
    def inf_results():
        sections = inf_sections.get()
        if not sections:
            return ui.p("No inference results yet.", class_="note")
        return ui.TagList(
            *[
                ui.TagList(ui.h5(title), ui.HTML(table_html))
                for title, table_html in sections
            ]
        )


app = App(app_ui, server)
