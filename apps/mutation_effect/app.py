from __future__ import annotations

import sys
from pathlib import Path

from shiny import App, reactive, render, ui

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from common.mutation_effect import (  # noqa: E402
    MutationEffectInference,
    MutationEffectService,
    MutationEffectViews,
)
from common.busy import (
    background,
    background_task_button,
    busy_indicators,
    finish_task,
    task_button,
    task_outcome,
)
from common.vendor import enable_compression, static_assets  # noqa: E402
from common.logging_utils import get_logger, new_trail  # noqa: E402
from common.ui_shell import (  # noqa: E402
    ACCESSION_LABEL,
    scop3p_card,
    scop3p_example_button,
    scop3p_field_row,
    scop3p_footer,
    scop3p_shell,
)


service = MutationEffectService()
LOGGER = get_logger("scop3p.mutation_effect")

#: Worked example: RET, with well-characterised phosphosites and disease mutations.
EXAMPLE_ACCESSION = "P07949"
EXAMPLE_POSITIONS = "606"
EXAMPLE_MUTATIONS = "A"


def _bokeh_iframe(html_doc: str) -> ui.Tag:
    return ui.tags.iframe(
        srcdoc=html_doc,
        style="width:100%;height:380px;border:0;border-radius:10px;background:#fff;",
    )


app_ui = scop3p_shell(
    "Mutation Effect",
    "Compare wild-type and mutant Bio2Byte predictions, keep Scop3P PTM context visible, and derive mutation-centric inference summaries from the same workflow.",
    ui.div(
        scop3p_card(
            "Session Status",
            ui.output_text_verbatim("status"),
            extra_class="scop3p-status",
        ),
        class_="scop3p-header-grid",
    ),
    ui.navset_tab(
        ui.nav_panel(
            "1) WT prediction",
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
                        scop3p_example_button("load_example"),
                    ),
                    background_task_button("run_wt", "Fetch + Predict WT", class_="btn btn-primary"),
                    ui.p("This runs the wild-type UniProt fetch, Scop3P PTM fetch, and Bio2Byte prediction.", class_="scop3p-note"),
                ),
                scop3p_card(
                    "WT Results",
                    ui.output_ui("wt_results"),
                ),
                class_="scop3p-two-col",
            ),
        ),
        ui.nav_panel(
            "2) Mutant prediction",
            ui.div(
                scop3p_card(
                    "Mutation Setup",
                    ui.input_text("positions", "Positions", value="", placeholder="e.g. 10,25,100"),
                    scop3p_field_row(
                        ui.input_text("mut_aas", "To amino acid", value="", placeholder="e.g. A,V,G"),
                        scop3p_example_button("load_example_mut"),
                    ),
                    task_button(
                        "run_mut", "Apply + Predict", class_="btn btn-warning"),
                    ui.p("Use comma-separated 1-indexed positions and amino-acid targets.", class_="scop3p-note"),
                ),
                scop3p_card(
                    "Mutant Results",
                    ui.output_ui("mut_results"),
                ),
                class_="scop3p-two-col",
            ),
        ),
        ui.nav_panel(
            "3) Inference",
            ui.div(
                scop3p_card(
                    "Label Shift Analysis",
                    task_button(
                        "run_inf", "Run inference", class_="btn btn-info"),
                    ui.p("Summarize class shifts at the mutation site and within a +/-5 residue window.", class_="scop3p-note"),
                ),
                scop3p_card(
                    "Inference Results",
                    ui.output_ui("inf_results"),
                ),
                class_="scop3p-two-col",
            ),
        ),
    ),
)


def server(input, output, session):
    # One trail per browser session: step numbers must not interleave across
    # sessions, and a module-level trail would be shared by every user.
    trail = new_trail()
    trail.opened("Mutation Effect")

    def _predict_wt(accession_value: str) -> dict:
        """The blocking half of the WT run: two API calls and a Bio2Byte prediction.

        Runs in a worker thread, so it must not touch reactive values -- everything it
        needs arrives as an argument and everything it produces comes back in this dict.
        """
        sequence_value = service.fetch_uniprot_sequence(accession_value)
        mods = service.fetch_scop3p_modifications(accession_value)
        prediction = service.predict_biophysical(accession_value, sequence_value)
        dataframe = service.prediction_to_df(prediction, accession_value)
        dataframe["seq"] = list(sequence_value)
        return {
            "accession": accession_value,
            "sequence": sequence_value,
            "mods": mods,
            "frame": dataframe,
        }

    _wt_task = background(_predict_wt)

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
    @reactive.event(input.load_example)
    def _load_example() -> None:
        trail.clicked("Load example")
        ui.update_text("accession", value=EXAMPLE_ACCESSION)
        status_text.set(
            f"Example accession {EXAMPLE_ACCESSION} loaded. Click Fetch + Predict WT."
        )

    @reactive.effect
    @reactive.event(input.load_example_mut)
    def _load_example_mut() -> None:
        trail.clicked("Load example mutation")
        ui.update_text("positions", value=EXAMPLE_POSITIONS)
        ui.update_text("mut_aas", value=EXAMPLE_MUTATIONS)
        status_text.set(
            f"Example mutation {EXAMPLE_ACCESSION} "
            f"{EXAMPLE_POSITIONS}->{EXAMPLE_MUTATIONS} loaded. Click Apply + Predict."
        )

    @reactive.effect
    @reactive.event(input.run_wt)
    def _run_wt() -> None:
        accession_value = input.accession().strip()
        trail.entered(ACCESSION_LABEL, accession_value or "-")
        trail.clicked("Fetch + Predict WT")
        LOGGER.info("run_wt requested accession=%s", accession_value or "-", extra={"event": "run_wt"})
        if not accession_value:
            trail.blocked("accession missing")
            status_text.set("Enter a UniProtKB accession first.")
            finish_task("run_wt")
            return
        status_text.set(
            f"Fetching {accession_value} and running the WT biophysical prediction. "
            "This takes a few seconds."
        )
        _wt_task(accession_value)

    @reactive.effect
    def _run_wt_done() -> None:
        def succeeded(payload: dict) -> None:
            dataframe, mods = payload["frame"], payload["mods"]
            accession.set(payload["accession"])
            sequence.set(payload["sequence"])
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
            trail.produced(
                f"WT prediction over {len(dataframe)} residues",
                ptms=0 if mods is None else len(mods),
            )
            LOGGER.info(
                "run_wt completed accession=%s seq_len=%s ptms=%s rows=%s",
                payload["accession"], len(payload["sequence"]),
                0 if mods is None else len(mods), len(dataframe),
                extra={"event": "run_wt"},
            )

        def failed(error: Exception) -> None:
            LOGGER.exception("run_wt failed", exc_info=error, extra={"event": "run_wt"})
            trail.failed("run_wt failed", error=type(error).__name__)
            status_text.set(f"WT error: {error}")

        task_outcome(
            _wt_task,
            on_success=succeeded,
            on_error=failed,
            on_finished=lambda: finish_task("run_wt"),
        )

    @reactive.effect
    @reactive.event(input.run_mut)
    def _run_mut() -> None:
        try:
            trail.clicked("Apply + Predict mutant")
            LOGGER.info("run_mut requested positions=%s aas=%s", input.positions(), input.mut_aas(), extra={"event": "run_mut"})
            if wt_df.get() is None or not sequence.get():
                trail.blocked("wt prediction missing")
                raise ValueError("Run WT prediction first.")

            parsed_mutations = service.parse_mutations(input.positions(), input.mut_aas())
            status_text.set("Applying mutations and predicting mutant...")

            mutant_sequence = service.apply_mutations(sequence.get(), parsed_mutations)
            # Not cached: a mutant is seen once, and caching it would evict the
            # wild-type prediction that every comparison needs.
            prediction = service.predict_biophysical(
                accession.get(), mutant_sequence, wild_type=False
            )
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
            LOGGER.exception("run_mut failed", extra={"event": "run_mut"})
            trail.failed("run_mut failed", error=type(error).__name__)
            status_text.set(f"Mutant error: {error}")
            return
        LOGGER.info(
            "run_mut completed mutations=%s rows=%s",
            len(parsed_mutations),
            len(dataframe),
            extra={"event": "run_mut"},
        )

    @reactive.effect
    @reactive.event(input.run_inf)
    def _run_inf() -> None:
        try:
            trail.clicked("Run inference")
            LOGGER.info("run_inf requested mutations=%s", len(mutations.get()), extra={"event": "run_inf"})
            if wt_df.get() is None or mut_df.get() is None:
                trail.blocked("required predictions missing")
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
            LOGGER.exception("run_inf failed", extra={"event": "run_inf"})
            trail.failed("run_inf failed", error=type(error).__name__)
            status_text.set(f"Inference error: {error}")
            return
        trail.produced(f"inference produced {len(sections)} section(s)")
        LOGGER.info("run_inf completed sections=%s", len(sections), extra={"event": "run_inf"})

    @output
    @render.ui
    def wt_results():
        if not wt_plot_html.get():
            return ui.p("No WT prediction yet.", class_="scop3p-note")
        return ui.TagList(
            ui.HTML(MutationEffectViews.track_guide_html()),
            _bokeh_iframe(wt_plot_html.get()),
            ui.hr(),
            ui.HTML(wt_table_html.get()),
        )

    # suspend_when_hidden=False because this output lives in a nav_panel that is not
    # the initially-active tab. Shiny decides suspension from the client-reported
    # ".clientdata_output_<id>_hidden" value, and Session._is_hidden() treats "never
    # reported" as hidden, so such an output is suspended at page load and is never
    # woken when the user opens its tab: it sits at "recalculating" forever with no
    # error logged anywhere. Verified against shiny 1.7.0, which requirements-shiny.txt
    # permits (shiny>=1.1,<2).
    @output(suspend_when_hidden=False)
    @render.ui
    def mut_results():
        if not mut_plot_html.get():
            return ui.p("No mutant prediction yet.", class_="scop3p-note")
        return ui.TagList(
            ui.HTML(mut_summary_html.get()),
            ui.HTML(MutationEffectViews.track_guide_html()),
            _bokeh_iframe(mut_plot_html.get()),
            ui.hr(),
            ui.HTML(mut_table_html.get()),
        )

    @output(suspend_when_hidden=False)
    @render.ui
    def inf_results():
        sections = inf_sections.get()
        if not sections:
            return ui.p("No inference results yet.", class_="scop3p-note")
        return ui.TagList(
            *[
                ui.TagList(ui.h5(title), ui.HTML(table_html))
                for title, table_html in sections
            ]
        )

content_ui = ui.div(
    busy_indicators(),
    app_ui, scop3p_footer()
)

# static_assets serves the vendored browser libraries; every app mounts the same prefix,
# so /vendor/... resolves whichever app the portal is serving. enable_compression is not
# optional cosmetics: Shiny sends static files raw, and molstar.js is 5 MB uncompressed
# against 1.45 MB gzipped, so without it vendoring would put more bytes on the wire.
app = App(content_ui, server, static_assets=static_assets())
enable_compression(app)
