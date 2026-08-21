"""RINAlign: residue interaction network alignment and comparison, as a Shiny app.

Ported from ``notebooks/RINAlign_align_and compare_networks.ipynb``. The science is
in :mod:`common.rinalign` and the browser code in :mod:`common.rinalign_views`; this
module is wiring only.

The notebook's ``_voila_iframe`` becomes ``ui.tags.iframe(srcdoc=...)``, the same
mechanism ``apps/mutation_effect/app.py`` already uses for Bokeh. Each script-bearing
view gets its own iframe so their duplicated element ids and ``window[...]`` helpers
cannot collide -- with one exception that matters: the linked view is a SINGLE iframe,
because its D3 force graph and its NGL stage must share one ``window`` for the
``__RIN_HL`` / ``__RIN_ONSELECT`` handshake to work.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

from shiny import App, reactive, render, ui

_APPS_DIR = Path(__file__).resolve().parents[1]
if str(_APPS_DIR) not in sys.path:
    sys.path.append(str(_APPS_DIR))

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from common import rinalign_views as views  # noqa: E402
from common.structure_labels import LOOKUP_FAILED_RETRY_FETCH
from common.busy import INERT_CSS, busy_indicators, gate, task_button
from common.vendor import enable_compression, static_assets  # noqa: E402
from common.logging_utils import get_logger, new_trail  # noqa: E402
from common.rinalign import (  # noqa: E402
    RINAlignService,
    align_rins,
    build_rin,
    diff_rins,
    rin_html,
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


LOGGER = get_logger("scop3p.rinalign")

AFDB_KEY = "af"
#: Worked example: p53, richly represented in the PDB and in Scop3P.
EXAMPLE_ACCESSION = "P04637"
VIEW_HEIGHTS = {"contact_map": 800, "aligned": 620, "force": 760, "linked": 880}

#: Views produced by each mode. Offering an alignment's node mapping while a diff is
#: loaded (or vice versa) only ever leads to a "not available in this mode" panel.
DIFF_VIEW_CHOICES = {
    "contact_map": "Contact map",
    "aligned": "Aligned overlay",
    "force": "Force network",
    "linked": "Linked 3D",
}
ALIGN_VIEW_CHOICES = {"mapping": "Node mapping"}


def _view_iframe(html_document: str, height_px: int) -> ui.Tag:
    """One view, one document, one browser ``window``.

    Generalises ``apps/mutation_effect/app.py::_bokeh_iframe``. Shiny escapes the
    ``srcdoc`` attribute itself, which is why ``rinalign_views.html_document`` does
    not escape -- doing both would double-escape the whole view.
    """
    return ui.tags.iframe(
        srcdoc=html_document,
        style=(
            f"width:100%;height:{height_px}px;border:0;"
            "border-radius:10px;background:#fff;"
        ),
    )


def _empty(message: str) -> ui.Tag:
    return ui.p(message, class_="scop3p-note")


_INTRO = (
    "Build residue interaction networks from two structures and compare them: diff "
    "two models of the same protein into conserved, lost and gained contacts, or "
    "align two different proteins by graph topology. Scop3P PTMs and UniProt disease "
    "variants can be overlaid on the network views."
)


app_ui = scop3p_shell(
    "RIN Alignment",
    _INTRO,
    ui.tags.style(
        """
.ra-button-row {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  align-items: center;
  margin: 10px 0;
}
.ra-rin-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 14px;
}
.ra-main-grid {
  display: grid;
  grid-template-columns: 420px minmax(0, 1fr);
  gap: 18px;
  align-items: stretch;
}
.ra-main-grid > .scop3p-card {
  height: 100%;
}
@media (max-width: 1200px) {
  .ra-main-grid { grid-template-columns: 1fr; }
}
/* Numbered steps inside the controls card, so the two stages of the workflow read
   as a sequence now that they are no longer separate tabs. */
.ra-step {
  margin: 0 0 12px;
  font-size: 0.82rem;
  text-transform: uppercase;
  letter-spacing: 0.08em;
  color: var(--scop3p-muted);
}
.ra-view-disabled {
  opacity: 0.45;
  pointer-events: none;
  user-select: none;
}
.ra-disabled-hint {
  margin: 4px 0 0;
  font-style: italic;
}
.ra-step-rule {
  margin: 20px 0 16px;
  border: 0;
  border-top: 1px solid var(--scop3p-line);
}
.ra-controls-card .scop3p-status pre,
.ra-controls-card pre {
  margin: 0;
  white-space: pre-wrap;
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 12px;
  background: rgba(16, 38, 60, 0.04);
  border-radius: 8px;
  padding: 8px 10px;
}
"""
    ),
    ui.div(
        scop3p_card(
            "Session Status",
            ui.output_text_verbatim("status"),
            extra_class="scop3p-status",
        ),
        scop3p_card("Protein", ui.output_ui("protein_info")),
        class_="scop3p-header-grid",
    ),
    # Two cards, not three tabs: the whole workflow is "set up on the left, look at
    # the result on the right", and tabs hid step 2 from anyone who had not thought to
    # look for it. The view selector lives with the networks it switches, and stays
    # disabled until a comparison exists to switch between.
    ui.div(
        scop3p_card(
            "Controls",
            ui.h5("1. Structure selection", class_="ra-step"),
            ui.input_radio_buttons(
                "mode",
                "Mode",
                {
                    "same": "Same protein (diff)",
                    "align": "Different proteins (align)",
                },
                selected="same",
                inline=False,
            ),
            scop3p_field_row(
                ui.input_text(
                    "accession",
                    ACCESSION_LABEL,
                    value="",
                    placeholder=f"e.g. {EXAMPLE_ACCESSION}",
                ),
                task_button(
                        "fetch", "Fetch", class_="btn btn-primary"),
                scop3p_example_button("load_example"),
            ),
            scop3p_structure_picker("left_structure", "Left (Model A)", {}),
            scop3p_structure_picker("right_structure", "Right (Model B)", {}),
            ui.hr(class_="ra-step-rule"),
            ui.h5("2. Annotation sources", class_="ra-step"),
            ui.input_checkbox("include_uniprot_ptms", "Include UniProt PTMs", value=True),
            ui.div(
                ui.input_action_button(
                    "fetch_ptm", "Fetch PTMs (Scop3P)", class_="btn btn-warning"
                ),
                ui.input_action_button(
                    "fetch_variants", "Fetch disease variants", class_="btn btn-info"
                ),
                class_="ra-button-row",
            ),
            ui.p(
                "Optional. Fetch before building, and the marks are overlaid on the "
                "network views.",
                class_="scop3p-note",
            ),
            ui.output_text_verbatim("ptmvar_status"),
            # Last, because this is the step you take once the two structures and any
            # annotations are chosen: pick the contact distance, build, compare.
            ui.hr(class_="ra-step-rule"),
            ui.h5("3. Build and compare", class_="ra-step"),
            ui.input_slider(
                "cutoff", "Contact cutoff (A)", min=4.0, max=14.0, value=8.0, step=0.5
            ),
            ui.div(
                task_button(
                        "generate", "Generate RINs", class_="btn btn-info"),
                ui.output_ui("compare_actions"),
                class_="ra-button-row",
            ),
            ui.p(
                "Changing a structure or the cutoff invalidates the existing networks, "
                "so Compare always reflects what the controls say.",
                class_="scop3p-note",
            ),
            extra_class="ra-controls-card",
        ),
        scop3p_card(
            "Networks",
            ui.div(
                ui.output_ui("left_rin_stats"),
                ui.output_ui("right_rin_stats"),
                class_="ra-rin-grid",
            ),
            ui.output_ui("summary_view"),
            # Rendered rather than static so it can be disabled until a comparison
            # exists. A radio selector rather than a nested navset: Shiny suspends
            # hidden outputs, and a navset nested inside another navset's panel does
            # not propagate its show event to Shiny's output-visibility bookkeeping,
            # so the panel becomes visible while the output stays "recalculating"
            # forever. One always-visible output switched by an input avoids that,
            # and still keeps exactly one heavy iframe alive at a time.
            ui.output_ui("view_selector"),
            ui.output_ui("comparison_view"),
        ),
        class_="ra-main-grid",
    ),
)

def server(input, output, session):
    # One trail per browser session: step numbers must not interleave across
    # sessions, and a module-level trail would be shared by every user.
    trail = new_trail()
    trail.opened("RIN Alignment")

    # Per session, not per process: the notebook's single ``app_state`` dict would
    # let two users overwrite each other's networks and annotation sets.
    service = RINAlignService(
        workdir=Path(tempfile.mkdtemp(prefix="scop3p_rinalign_")), timeout=30
    )

    status_text = reactive.value("Enter a UniProtKB accession to start.")
    ptmvar_text = reactive.value("No PTM or variant data fetched.")

    accession = reactive.value("")
    protein_html = reactive.value("")
    structure_entries = reactive.value({})

    graph_left = reactive.value(None)
    graph_right = reactive.value(None)
    rin_left_html = reactive.value("")
    rin_right_html = reactive.value("")
    rins_ready = reactive.value(False)
    rin_signature = reactive.value(None)
    labels = reactive.value(("L", "R"))

    left_structure_text = reactive.value(None)
    left_structure_format = reactive.value("pdb")
    left_chain = reactive.value(None)

    ptm_positions = reactive.value(frozenset())
    variant_positions = reactive.value(frozenset())

    # One dict written atomically, so a half-updated mixture of views can never
    # render: all five come from the same comparison or none do.
    compare_payload = reactive.value(None)

    def say(message: str) -> None:
        status_text.set(message)

    def fail(event: str, error: Exception) -> None:
        LOGGER.exception("%s failed", event, extra={"event": event})
        trail.failed(f"{event} failed", error=type(error).__name__)
        status_text.set(f"{type(error).__name__}: {error}")
        if event == "fetch_protein":
            # Say it in the pickers too. The status card explains what happened, but the
            # empty dropdown is where the user is looking, and an empty dropdown alone
            # gives no hint that pressing Fetch again is worth trying.
            for input_id in ("left_structure", "right_structure"):
                ui.update_selectize(
                    input_id, choices={"": LOOKUP_FAILED_RETRY_FETCH}, selected=""
                )

    def clear_networks(message: str | None = None) -> None:
        graph_left.set(None)
        graph_right.set(None)
        rin_left_html.set("")
        rin_right_html.set("")
        rins_ready.set(False)
        rin_signature.set(None)
        left_structure_text.set(None)
        compare_payload.set(None)
        if message:
            say(message)

    def selected_entries():
        entries = structure_entries.get()
        left = entries.get(input.left_structure())
        right = entries.get(input.right_structure())
        return left, right

    # -- handlers -----------------------------------------------------------

    @reactive.effect
    @reactive.event(input.load_example)
    def _load_example() -> None:
        trail.clicked("Load example")
        ui.update_text("accession", value=EXAMPLE_ACCESSION)
        say(f"Example accession {EXAMPLE_ACCESSION} loaded. Click Fetch.")

    @reactive.effect
    @reactive.event(input.fetch)
    def _fetch() -> None:
        trail.clicked("Fetch")
        value = input.accession().strip().upper()
        trail.entered(ACCESSION_LABEL, value or "-")
        if not value:
            say(f"Enter a UniProtKB accession, for example {EXAMPLE_ACCESSION}.")
            return
        try:
            say(f"Fetching {value}...")
            info = service.fetch_uniprot_info(value)
            accession.set(value)
            protein_html.set(
                "<div class='summary-banner'>"
                f"<b>{info['protein_name'] or value}</b><br>"
                f"Gene: <b>{info['gene_name'] or '-'}</b> | "
                f"Length: <b>{info['length']}</b> aa | "
                f"<i>{info['organism'] or '-'}</i> | <code>{info['accession']}</code>"
                "</div>"
            )

            say("Fetching PDB structures...")
            structures = service.fetch_pdb_structures(value, uniprot_data=info.get("_raw"))
            say("Checking AlphaFold...")
            alphafold = service.check_alphafold(value)

            ordered = ([alphafold] if alphafold else []) + list(structures)
            if not ordered:
                structure_entries.set({})
                ui.update_selectize("left_structure", choices={})
                ui.update_selectize("right_structure", choices={})
                clear_networks()
                say(
                    f"No structures found for {value}. It has no PDB entries and no "
                    "AlphaFold model."
                )
                return

            # Shiny select values must be plain strings; the notebook stored whole
            # dicts as Dropdown values, which has no equivalent here.
            keyed = {str(index): entry for index, entry in enumerate(ordered)}
            structure_entries.set(keyed)
            choices = {key: entry.label for key, entry in keyed.items()}

            ui.update_selectize("left_structure", choices=choices, selected="0")
            ui.update_selectize(
                "right_structure",
                choices=choices,
                selected="1" if (alphafold and len(ordered) > 1) else "0",
            )

            clear_networks()
            say(
                f"Found {len(structures)} PDB chain(s)"
                + (" + AlphaFold" if alphafold else "")
                + ". Select two structures and click Generate RINs."
            )
        except Exception as error:  # noqa: BLE001
            fail("fetch_protein", error)
            return
        LOGGER.info(
            "fetch_protein completed accession=%s structures=%s alphafold=%s",
            value,
            len(structures),
            bool(alphafold),
            extra={"event": "fetch_protein"},
        )

    @reactive.effect
    @reactive.event(input.generate)
    def _generate() -> None:
        trail.clicked("Generate RINs")
        entries = structure_entries.get()
        if not entries:
            say("Fetch a UniProtKB accession first.")
            return

        left_entry, right_entry = selected_entries()
        if left_entry is None or right_entry is None:
            say("Select a structure for both Left and Right.")
            return

        cutoff = float(input.cutoff())
        try:
            built = {}
            for side, entry in (("left", left_entry), ("right", right_entry)):
                say(f"Downloading {entry.pdb_id}...")
                path, text = service.download_structure(entry)
                say(f"Building RIN for {entry.display} at {cutoff} A...")
                graph, _residues = build_rin(path, chain_id=entry.chain_id, cutoff=cutoff)
                built[side] = (graph, entry, path, text)

            left_graph, left_entry, left_path, left_text = built["left"]
            right_graph, right_entry, _rp, _rt = built["right"]

            graph_left.set(left_graph)
            graph_right.set(right_graph)
            rin_left_html.set(rin_html(left_graph, f"Left · {left_entry.display}"))
            rin_right_html.set(rin_html(right_graph, f"Right · {right_entry.display}"))
            labels.set((left_entry.display, right_entry.display))

            # The linked view renders the LEFT structure in 3D.
            left_structure_text.set(left_text)
            left_structure_format.set("cif" if left_path.suffix == ".cif" else "pdb")
            left_chain.set(left_entry.chain_id)

            # Only now: the notebook enabled Compare after the left side succeeded,
            # then crashed on GR.number_of_nodes() if the right side had failed.
            rins_ready.set(True)
            rin_signature.set(
                (input.left_structure(), input.right_structure(), cutoff)
            )
            compare_payload.set(None)
            say(
                f"Networks ready. Left {left_graph.number_of_nodes()} nodes / "
                f"{left_graph.number_of_edges()} edges, right "
                f"{right_graph.number_of_nodes()} nodes / "
                f"{right_graph.number_of_edges()} edges. Click Compare."
            )
        except Exception as error:  # noqa: BLE001
            rins_ready.set(False)
            fail("generate_rins", error)
            return
        trail.produced(
            "RINs generated for both structures",
            cutoff=cutoff, left=left_entry.key, right=right_entry.key,
        )
        LOGGER.info(
            "generate_rins completed cutoff=%s left=%s right=%s",
            cutoff,
            left_entry.key,
            right_entry.key,
            extra={"event": "generate_rins"},
        )

    @reactive.effect
    @reactive.event(input.compare)
    def _compare() -> None:
        trail.clicked("Compare / Align")
        left_graph, right_graph = graph_left.get(), graph_right.get()
        if left_graph is None or right_graph is None:
            say("Generate both networks first.")
            return

        label_left, label_right = labels.get()
        ptm = set(ptm_positions.get())
        variants = set(variant_positions.get())
        mode = input.mode()

        try:
            say("Comparing...")
            if mode == "same":
                result = diff_rins(left_graph, right_graph)
                structure_text = left_structure_text.get()
                compare_payload.set(
                    {
                        "mode": "same",
                        "summary": views.summary_html(result, label_left, label_right),
                        "contact_map": views.html_document(
                            views.contact_map_html(
                                result, label_left, label_right,
                                ptm_pos=ptm, var_pos=variants,
                            )
                        ),
                        "aligned": views.html_document(
                            views.aligned_network_html(
                                result, left_graph, right_graph, label_left, label_right,
                                ptm_pos=ptm, var_pos=variants,
                            )
                        ),
                        "force": views.html_document(
                            views.force_network_html(
                                result, left_graph, right_graph, label_left, label_right,
                                ptm_pos=ptm, var_pos=variants,
                            )
                        ),
                        # One document: the force graph and the NGL stage have to
                        # share a window for the two-way highlight bridge.
                        "linked": views.html_document(
                            views.linked_view_html(
                                result, left_graph, right_graph, label_left, label_right,
                                structure_text, left_structure_format.get(),
                                left_chain.get(), ptm_pos=ptm, var_pos=variants,
                            )
                        )
                        if structure_text
                        else "",
                    }
                )
                say(
                    f"Done. Jaccard {result['jaccard']:.3f} | conserved "
                    f"{len(result['conserved'])} | lost {len(result['lost'])} | "
                    f"gained {len(result['gained'])} | mutations "
                    f"{len(result['mutations'])}."
                )
            else:
                result = align_rins(left_graph, right_graph)
                compare_payload.set(
                    {
                        "mode": "align",
                        "summary": (
                            f"<h4>Graph alignment: {label_left} vs {label_right}</h4>"
                            f"<p>Jaccard: <b>{result['jaccard']:.3f}</b> | "
                            f"Conserved edges: <b>{len(result['conserved'])}</b> | "
                            f"Mapped nodes: <b>{len(result['mapping'])}</b></p>"
                        ),
                        "mapping": _mapping_table_html(result),
                    }
                )
                say(f"Alignment done. Jaccard {result['jaccard']:.3f}.")
        except Exception as error:  # noqa: BLE001
            fail("compare_rins", error)
            return
        trail.produced(
            f"networks compared ({mode})", jaccard=round(result["jaccard"], 4),
        )
        LOGGER.info(
            "compare_rins completed mode=%s jaccard=%.4f",
            mode,
            result["jaccard"],
            extra={"event": "compare_rins"},
        )

    @reactive.effect
    @reactive.event(input.fetch_ptm)
    def _fetch_ptm() -> None:
        trail.clicked("Fetch PTMs")
        value = accession.get()
        if not value:
            ptmvar_text.set("Fetch a UniProtKB accession first.")
            return
        try:
            ptmvar_text.set("Fetching PTMs...")
            scop3p = service.fetch_scop3p_ptm_positions(value)
            uniprot: set[int] = set()
            if input.include_uniprot_ptms():
                uniprot = service.fetch_uniprot_ptm_positions(value)
            combined = set(scop3p) | uniprot
            ptm_positions.set(frozenset(combined))

            lines = [
                f"PTMs - Scop3P: {len(scop3p)} | UniProt: {len(uniprot)} | "
                f"total unique: {len(combined)}"
            ]
            if not combined:
                lines.append(
                    "No PTMs found. Scop3P mainly covers human phosphoproteins."
                )
            lines.append("Re-run Compare to overlay PTM rings on the network.")
            ptmvar_text.set("\n".join(lines))
        except Exception as error:  # noqa: BLE001
            LOGGER.exception("fetch_ptm failed", extra={"event": "fetch_ptm"})
            trail.failed("fetch_ptm failed", error=type(error).__name__)
            ptmvar_text.set(f"PTM error: {error}")
            return
        LOGGER.info(
            "fetch_ptm completed accession=%s scop3p=%s uniprot=%s",
            value,
            len(scop3p),
            len(uniprot),
            extra={"event": "fetch_ptm"},
        )

    @reactive.effect
    @reactive.event(input.fetch_variants)
    def _fetch_variants() -> None:
        trail.clicked("Fetch disease-associated variants")
        value = accession.get()
        if not value:
            ptmvar_text.set("Fetch a UniProtKB accession first.")
            return
        try:
            ptmvar_text.set("Fetching disease-associated variants...")
            positions = service.fetch_uniprot_variant_positions(value, disease_only=True)
            variant_positions.set(frozenset(positions))
            lines = [f"Disease-associated variants: {len(positions)}"]
            if not positions:
                lines.append("No disease-associated variants found for this protein.")
            lines.append("Re-run Compare to overlay variant rings on the network.")
            ptmvar_text.set("\n".join(lines))
        except Exception as error:  # noqa: BLE001
            LOGGER.exception("fetch_variants failed", extra={"event": "fetch_variants"})
            trail.failed("fetch_variants failed", error=type(error).__name__)
            ptmvar_text.set(f"Variant error: {error}")
            return
        LOGGER.info(
            "fetch_variants completed accession=%s variants=%s",
            value,
            len(positions),
            extra={"event": "fetch_variants"},
        )

    @reactive.effect
    def _invalidate_stale_networks() -> None:
        """Disable Compare when the controls no longer describe the built networks.

        No notebook equivalent: there, moving the cutoff slider left Compare enabled
        over graphs built at the old cutoff, so the results silently disagreed with
        the UI. Modelled on ``apps/structure_viz/app.py``'s stale-input effect.
        """
        signature = (input.left_structure(), input.right_structure(), float(input.cutoff()))
        loaded = rin_signature.get()
        if loaded is None or signature == loaded:
            return
        clear_networks(
            "Structure selection or cutoff changed. Re-run Generate RINs."
        )

    @reactive.effect
    @reactive.event(input.mode)
    def _mode_changed() -> None:
        # A "same protein" payload must not linger while the UI says "align".
        if compare_payload.get() is not None:
            compare_payload.set(None)
            say("Mode changed. Re-run Compare.")

    # -- outputs ------------------------------------------------------------

    @render.text
    def status() -> str:
        return status_text.get()

    # Outputs inside a nav_panel that is not the initially-active one must opt out
    # of hidden-output suspension. Shiny decides suspension from the client-reported
    # ``.clientdata_output_<id>_hidden`` value, and Session._is_hidden() treats
    # "never reported" as hidden -- so an output in a tab that starts inactive is
    # suspended at load and is never woken when the user opens that tab. It sits at
    # "recalculating" forever, with no error anywhere. suspend_when_hidden=False
    # takes the output out of that mechanism entirely.
    @output(suspend_when_hidden=False)
    @render.text
    def ptmvar_status() -> str:
        return ptmvar_text.get()

    @render.ui
    def protein_info():
        if not protein_html.get():
            return _empty("Fetch a UniProtKB accession to see protein details.")
        return ui.HTML(protein_html.get())

    @render.ui
    def left_rin_stats():
        if not rin_left_html.get():
            return _empty("No network built yet.")
        return ui.HTML(rin_left_html.get())

    @render.ui
    def right_rin_stats():
        if not rin_right_html.get():
            return _empty("No network built yet.")
        return ui.HTML(rin_right_html.get())

    @render.ui
    def compare_actions():
        # Rendered rather than static so it can be disabled until both networks
        # exist, mirroring the notebook's ``cmp_btn.disabled`` flag.
        return gate(
            task_button("compare", "Compare / Align", class_="btn btn-success"),
            ready=rins_ready.get(),
            hint="Generate both networks first.",
        )

    @render.ui
    def view_selector():
        """The view radio buttons, disabled until there is something to switch between.

        Offering them before Compare has run invites a click that can only report
        "no comparison yet", so they stay inert until a payload exists. The choices are
        also trimmed to the modes that apply: a diff produces the four network views, an
        alignment produces the node mapping.
        """
        payload = compare_payload.get()
        if payload is None:
            # Shiny has no `disabled` for a radio group, so the wrapper carries a
            # class that the CSS makes inert. The server ignores the value anyway
            # while compare_payload is None.
            return ui.div(
                ui.input_radio_buttons(
                    "view",
                    "View",
                    DIFF_VIEW_CHOICES | ALIGN_VIEW_CHOICES,
                    selected="contact_map",
                    inline=True,
                ),
                ui.p(
                    "Available once Compare / Align has run.",
                    class_="scop3p-note ra-disabled-hint",
                ),
                class_="ra-view-disabled",
            )
        choices = DIFF_VIEW_CHOICES if payload["mode"] == "same" else ALIGN_VIEW_CHOICES
        current = input.view() if input.view() in choices else next(iter(choices))
        return ui.input_radio_buttons(
            "view", "View", choices, selected=current, inline=True
        )

    @output(suspend_when_hidden=False)
    @render.ui
    def summary_view():
        payload = compare_payload.get()
        if payload is None:
            # Silent: comparison_view carries the one "run Compare first" message, so
            # the card does not say the same thing twice.
            return None
        return ui.HTML(payload["summary"])

    @output(suspend_when_hidden=False)
    @render.ui
    def comparison_view():
        payload = compare_payload.get()
        key = input.view()

        if payload is None:
            return _empty("Generate both networks, then click Compare / Align.")

        if key == "mapping":
            if payload["mode"] != "align":
                return _empty(
                    "Node mapping is produced by 'Different proteins (align)' mode."
                )
            return ui.HTML(payload["mapping"])

        if payload["mode"] != "same":
            return _empty(
                "This view is only available in 'Same protein (diff)' mode, which "
                "compares matched residue positions. Use Node mapping for an "
                "alignment result."
            )

        document = payload.get(key) or ""
        if not document:
            return _empty(
                "The left structure's coordinates were not available, so the linked "
                "3D view cannot be drawn."
            )
        return _view_iframe(document, VIEW_HEIGHTS[key])


def _mapping_table_html(result: dict, limit: int = 20) -> str:
    """Top-scoring node correspondences from the Hungarian assignment."""
    rows = sorted(result["mapping"], key=lambda item: -item[2])[:limit]
    cells = "".join(
        "<tr>"
        f"<td style='padding:4px 10px;'>{first}</td>"
        f"<td style='padding:4px 10px;'>{second}</td>"
        f"<td style='padding:4px 10px;text-align:center;'>{score:.3f}</td>"
        "</tr>"
        for first, second, score in rows
    )
    return (
        f"<h5>Node mapping (top {len(rows)} of {len(result['mapping'])})</h5>"
        "<table style='border-collapse:collapse;font-size:13px;'>"
        "<tr style='background:#f5f5f5;'>"
        "<th style='padding:4px 10px;'>Network A</th>"
        "<th style='padding:4px 10px;'>Network B</th>"
        "<th style='padding:4px 10px;'>Score</th></tr>"
        f"{cells}</table>"
    )


content_ui = ui.div(
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
