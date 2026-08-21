"""Protein topology viewer, as a Shiny app.

The science lives in the ``topology`` package at ``notebooks/topology_viewer/``,
reached through :mod:`common.topology_bridge`. This module is the control layer only:
it replaces the ``ipywidgets`` front end in ``topology/app.py::make_app`` with Shiny
reactives and leaves ``build_view`` -- the whole 2D/3D pipeline and every line of
browser code -- untouched.

Two modes, deliberately kept apart, exactly as the package docstring requires:

    accession   fetch from AlphaFold DB or PDBe, with PTM and variant overlays
    file        an uploaded .pdb/.cif; no network, and no UniProt, PTM or variant
                pipeline at all

Conflating them is the bug the package was restructured to fix: an upload named
``P07949_relaxed_rank_1.pdb`` used to draw the prediction's topology beside
AlphaFold's *different* coordinates in 3D, with nothing to signal the mismatch. The
guard here is ``loaded_mode``, set by whichever handler actually loaded a structure --
never ``input.mode()``, because ``panel_conditional`` hides inputs without destroying
them, so a hidden accession field stays readable.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

from shiny import App, reactive, render, ui

_APPS_DIR = Path(__file__).resolve().parents[1]
if str(_APPS_DIR) not in sys.path:
    sys.path.append(str(_APPS_DIR))

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from common.busy import busy_indicators, task_button
from common.logging_utils import get_logger, new_trail
from common.structure_labels import ALPHAFOLD_OPTION_LABEL, structure_option_label  # noqa: E402
from common.topology_bridge import (  # noqa: E402
    TOPOLOGY_ERROR,
    annotations_module as ann,
    build_view,
    fetch_alphafold,
    load_structure,
    __topology_version__,
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


LOGGER = get_logger("scop3p.topology_viewer")

DEFAULT_ACCESSION = "P07949"
VIEWER_HEIGHT = int(os.getenv("SCOP3P_TOPOLOGY_HEIGHT", "1150"))
AFDB_CHOICE = "afdb"


def preferred_chain(refs, pdb_id: str, structure) -> str | None:
    """The chain this accession actually maps to, not merely the biggest one.

    A complex often carries a larger unrelated chain, and defaulting to it would draw
    a topology for the wrong protein. Lifted out of the widget closure so it can be
    unit-tested against a stub ``StructureRef``.
    """
    for ref in refs or []:
        if ref.pdb_id.upper() != pdb_id.upper():
            continue
        for chain in ref.chains:
            if chain in structure.residues_by_chain:
                return chain
    return None


def _chain_choices(structure) -> dict[str, str]:
    # Structure.chain_options() yields (label, value); Shiny wants value -> label.
    return {value: label for label, value in structure.chain_options()}


_INTRO = (
    "Draw a secondary-structure topology diagram beside the 3D structure. Fetch an "
    "AlphaFold DB model or a PDBe entry by UniProtKB accession, or upload a predicted "
    "structure and stay offline. Secondary structure is read from the file when it is "
    "present and derived from coordinates when it is not."
)


def _controls_card() -> ui.Tag:
    return scop3p_card(
        "Structure Source",
        ui.input_radio_buttons(
            "mode",
            "Source",
            {
                "accession": "Fetch from AlphaFold DB / PDBe",
                "file": "Upload a structure file",
            },
            selected="accession",
            inline=True,
        ),
        ui.panel_conditional(
            "input.mode === 'accession'",
            scop3p_field_row(
                ui.input_text(
                    "accession",
                    ACCESSION_LABEL,
                    value="",
                    placeholder=f"e.g. {DEFAULT_ACCESSION}",
                ),
                scop3p_example_button("load_example"),
            ),
            task_button(
                        "fetch_btn", "Fetch model", class_="btn btn-primary tv-block-btn"
            ),
            scop3p_structure_picker(
                "structure_choice", "Structure", {AFDB_CHOICE: ALPHAFOLD_OPTION_LABEL}
            ),
            ui.div(
                task_button(
                        "fetch_ptms", "Fetch PTMs", class_="btn btn-info"),
                task_button(
                        "fetch_variants", "Fetch variants", class_="btn btn-info"
                ),
                ui.input_action_button(
                    "clear_annotations", "Clear", class_="btn btn-danger"
                ),
                class_="tv-button-row",
            ),
            ui.input_checkbox("include_uniprot_ptms", "Include UniProt PTMs", value=True),
            ui.p(
                "Scop3P imports its mutations from UniProt, so UniProt is the primary "
                "source for disease variants.",
                class_="scop3p-note",
            ),
        ),
        ui.panel_conditional(
            "input.mode === 'file'",
            ui.input_file(
                "structure_upload",
                "Structure file",
                accept=[".pdb", ".ent", ".cif", ".mmcif"],
                multiple=False,
            ),
            ui.p(
                "File mode never contacts the network: the 3D panel is fed the bytes "
                "you uploaded, and no PTM or variant overlay is offered, because a "
                "local prediction has no reliable UniProt numbering.",
                class_="scop3p-note",
            ),
        ),
        ui.input_select("chain", "Chain", choices={}),
    )


app_ui = scop3p_shell(
    "Topology Viewer",
    _INTRO,
    ui.tags.style(
        """
.tv-block-btn {
  width: 100%;
  margin-bottom: 12px;
}
.tv-button-row {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  margin-bottom: 12px;
}
/* Session Status and Notes stack in the left header column. */
.tv-status-stack {
  display: flex;
  flex-direction: column;
  gap: 18px;
}
/* Controls beside the diagram. minmax(0, 1fr) rather than 1fr so the 1150px iframe
   cannot force the column wider than the grid and push the page into horizontal
   overflow. */
.tv-main-grid {
  display: grid;
  grid-template-columns: 360px minmax(0, 1fr);
  gap: 18px;
  align-items: stretch;
}
.tv-main-grid > .scop3p-card {
  height: 100%;
}
@media (max-width: 1200px) {
  .tv-main-grid { grid-template-columns: 1fr; }
}
"""
    ),
    ui.div(
        ui.div(
            scop3p_card(
                "Session Status",
                ui.output_text_verbatim("status"),
                extra_class="scop3p-status",
            ),
            scop3p_card("Notes", ui.output_ui("notes_panel")),
            class_="tv-status-stack",
        ),
        scop3p_card("Loaded Model", ui.output_ui("model_summary")),
        class_="scop3p-header-grid",
    ),
    # Controls and diagram side by side. The diagram is the reason the controls column
    # is kept narrow: build_view() returns an iframe containing its own 2D/3D
    # two-column grid, so whatever width this column takes is subtracted from both
    # panels inside it.
    ui.div(
        _controls_card(),
        scop3p_card("Topology", ui.output_ui("topology_view")),
        class_="tv-main-grid",
    ),
)


def server(input, output, session):
    # One trail per browser session: step numbers must not interleave across
    # sessions, and a module-level trail would be shared by every user.
    trail = new_trail()
    trail.opened("Topology Viewer")

    status_text = reactive.value(
        "Enter a UniProtKB accession and fetch a model, or upload a .pdb or .cif file."
    )
    summary_bits = reactive.value([])

    structure = reactive.value(None)
    source = reactive.value({})
    loaded_mode = reactive.value("")
    accession_loaded = reactive.value("")

    refs = reactive.value([])
    ptm_sites = reactive.value([])
    variant_sites = reactive.value([])
    notes = reactive.value([])
    numbering = reactive.value(None)
    numbering_source = reactive.value("")

    def say(message: str) -> None:
        status_text.set(message)

    def fail(event: str, error: Exception) -> None:
        LOGGER.exception("%s failed", event, extra={"event": event})
        trail.failed(f"{event} failed", error=type(error).__name__)
        status_text.set(f"{type(error).__name__}: {error}")

    def unavailable() -> bool:
        if TOPOLOGY_ERROR:
            status_text.set(TOPOLOGY_ERROR)
            return True
        return False

    def adopt(structure_value, source_value: dict, label: str, chain: str | None = None) -> None:
        """Take ownership of a freshly loaded structure and refresh the chain picker.

        The ipywidgets original needed a ``suspend_redraw`` flag here, because
        assigning ``.options`` then ``.value`` re-entered the observer and redrew from
        half-updated state -- and calling ``unobserve_all()`` to avoid that removed
        ipywidgets' own internal options observer, so the next ``.value`` assignment
        raised ``TraitError`` and the callback swallowed it, freezing the app.

        Neither problem exists here. ``ui.update_select`` carries choices and
        selection in one message with no re-validation of the stale value, and
        everything this function writes lands in a single reactive flush, so
        ``topology_view`` re-renders exactly once. Do not reintroduce the flag.
        """
        structure.set(structure_value)
        source.set(source_value)

        wanted = chain if chain in structure_value.residues_by_chain else None
        ui.update_select(
            "chain",
            choices=_chain_choices(structure_value),
            selected=wanted or structure_value.default_chain(),
        )

        bits = [
            label,
            f"{len(structure_value.chains)} chain(s)",
            f"SS from {structure_value.ss_source}",
        ]
        if numbering_source.get():
            bits.append(f"numbering: {numbering_source.get()}")
        if ptm_sites.get():
            bits.append(f"{len(ptm_sites.get())} PTMs")
        if variant_sites.get():
            bits.append(f"{len(variant_sites.get())} variants")
        summary_bits.set(bits)
        say(" · ".join(bits))

    def load_selected_structure() -> None:
        """Load whichever entry the structure picker points at.

        Numbering is resolved here rather than at draw time, because the choice of
        entry is exactly what determines it: AlphaFold models are already
        UniProt-numbered, while a PDB entry needs its SIFTS map or every annotation
        lands on the wrong residue.
        """
        accession = accession_loaded.get()
        choice = input.structure_choice()

        if choice == AFDB_CHOICE:
            _record, text, cif_url = fetch_alphafold(accession)
            structure_value = load_structure(text, f"AF-{accession}-F1.cif")
            structure_value.uniprot = accession
            numbering.set(None)  # identity: positions already match
            numbering_source.set("AlphaFold model is UniProt-numbered")
            adopt(
                structure_value,
                {
                    "kind": "afdb",
                    "url": cif_url,
                    "format": "mmcif",
                    "accession": accession,
                },
                f"AlphaFold model for {accession}",
            )
            return

        pdb_id = choice
        say(f"Loading PDB entry {pdb_id}.")
        text = ann.fetch_structure_file(pdb_id)
        structure_value = load_structure(text, f"{pdb_id}.cif")
        structure_value.uniprot = accession

        # Chains that genuinely map to this accession, per the file's own SIFTS
        # columns, rather than whichever chain happens to be largest.
        mapped_chains = structure_value.chains_for_accession(accession)
        preferred = (mapped_chains or [None])[0] or preferred_chain(
            refs.get(), pdb_id, structure_value
        )
        numbering.set(None)
        numbering_source.set("")

        if structure_value.has_sifts:
            # The updated mmCIF carries the UniProt correspondence alongside the
            # coordinates, so it cannot disagree with what is being drawn, and it
            # costs no extra request.
            mapping = structure_value.sifts_numbering(preferred or "", accession)
            if mapping:
                numbering.set(mapping)
                numbering_source.set("SIFTS columns in the mmCIF")

        if numbering.get() is None:
            try:
                numbering.set(ann.fetch_numbering(pdb_id, accession, preferred))
                numbering_source.set("PDBe mapping API")
            except Exception as error:  # noqa: BLE001 - reported to the user
                LOGGER.warning(
                    "numbering unavailable pdb_id=%s accession=%s error=%s",
                    pdb_id,
                    accession,
                    error,
                    extra={"event": "fetch_numbering"},
                )
                notes.set(
                    notes.get()
                    + [
                        f"No numbering map for {pdb_id} ({error}). Sites are hidden "
                        "rather than drawn at UniProt positions, which would place "
                        "them on the wrong residues."
                    ]
                )

        adopt(
            structure_value,
            {
                "kind": "pdbe",
                "url": ann.structure_file_url(pdb_id),
                "format": "mmcif",
                "accession": accession,
                "pdb_id": pdb_id,
            },
            f"PDB entry {pdb_id}",
            chain=preferred,
        )

    @reactive.effect
    @reactive.event(input.load_example)
    def _load_example() -> None:
        trail.clicked("Load example")
        ui.update_text("accession", value=DEFAULT_ACCESSION)
        say(f"Example accession {DEFAULT_ACCESSION} loaded. Click Fetch model.")

    @reactive.effect
    @reactive.event(input.fetch_btn)
    def _fetch() -> None:
        trail.clicked("Fetch model")
        if unavailable():
            return
        accession = (input.accession() or DEFAULT_ACCESSION).strip().upper()
        trail.entered(ACCESSION_LABEL, accession)
        try:
            accession_loaded.set(accession)
            loaded_mode.set("accession")
            notes.set([])
            numbering.set(None)
            numbering_source.set("")
            say(f"Looking up {accession}.")

            collected: list = []
            try:
                collected = ann.fetch_structures(accession)
            except Exception as error:  # noqa: BLE001 - degrade, do not abort
                LOGGER.warning(
                    "structure list unavailable accession=%s error=%s",
                    accession,
                    error,
                    extra={"event": "fetch_structures"},
                )
                notes.set([f"Structure list unavailable ({error})."])
            refs.set(collected)

            # Labels come from common.structure_labels so an entry reads the same here
            # as in Structure Visualisation and RIN Alignment. StructureRef.label() is
            # left alone deliberately: that dataclass is shared with the notebook and its
            # test suite, and this is a presentation concern belonging to the app.
            choices = {AFDB_CHOICE: ALPHAFOLD_OPTION_LABEL}
            for ref in collected:
                choices[ref.pdb_id] = structure_option_label(
                    ref.pdb_id,
                    method=ref.method,
                    resolution=ref.resolution,
                    chains=ref.chains,
                    coverage=ref.coverage,
                )
            ui.update_selectize("structure_choice", choices=choices, selected=AFDB_CHOICE)

            load_selected_structure()
        except Exception as error:  # noqa: BLE001
            fail("fetch_model", error)
            return
        LOGGER.info(
            "fetch_model completed accession=%s structures=%s",
            accession,
            len(refs.get()),
            extra={"event": "fetch_model"},
        )

    @reactive.effect
    @reactive.event(input.structure_choice)
    def _structure_changed() -> None:
        # Fires on initial render too; there is nothing to load until Fetch has run.
        if TOPOLOGY_ERROR or not accession_loaded.get():
            return
        try:
            load_selected_structure()
        except Exception as error:  # noqa: BLE001
            fail("structure_changed", error)

    @reactive.effect
    @reactive.event(input.structure_upload)
    def _upload() -> None:
        if unavailable():
            return
        uploaded = input.structure_upload()
        if not uploaded:
            return
        item = uploaded[0]
        name = item.get("name") or "uploaded"
        try:
            # Shiny hands us a file on disk, so no base64/bytes decoding is needed.
            text = Path(item["datapath"]).read_text(encoding="utf-8", errors="replace")
            say(f"Reading {name}.")
            structure_value = load_structure(text, name)

            # File mode stays offline and annotation-free: the 3D view is fed the
            # uploaded bytes, never a database model that merely shares a name.
            loaded_mode.set("file")
            accession_loaded.set("")
            refs.set([])
            ptm_sites.set([])
            variant_sites.set([])
            notes.set([])
            numbering.set(None)
            numbering_source.set("")
            adopt(
                structure_value,
                {"kind": "upload", "data": text, "format": structure_value.fmt},
                f"Uploaded {name}",
            )
        except Exception as error:  # noqa: BLE001
            fail("upload_structure", error)
            return
        LOGGER.info(
            "upload_structure completed name=%s chains=%s ss_source=%s",
            name,
            len(structure_value.chains),
            structure_value.ss_source,
            extra={"event": "upload_structure"},
        )

    @reactive.effect
    @reactive.event(input.fetch_ptms)
    def _fetch_ptms() -> None:
        trail.clicked("Fetch PTMs")
        if unavailable():
            return
        accession = accession_loaded.get()
        if not accession:
            say("Fetch a model first, so there is an accession to look up.")
            return
        try:
            say(f"Fetching PTMs for {accession}.")
            sites, site_notes = ann.fetch_ptms(
                accession, include_uniprot=input.include_uniprot_ptms()
            )
            ptm_sites.set(sites)
            notes.set(site_notes)
            sources = sorted({site.source for site in sites})
            say(
                f"{len(sites)} PTM sites"
                + (f" from {', '.join(sources)}" if sources else "")
            )
        except Exception as error:  # noqa: BLE001
            fail("fetch_ptms", error)
            return
        LOGGER.info(
            "fetch_ptms completed accession=%s sites=%s",
            accession,
            len(sites),
            extra={"event": "fetch_ptms"},
        )

    @reactive.effect
    @reactive.event(input.fetch_variants)
    def _fetch_variants() -> None:
        trail.clicked("Fetch disease-associated variants")
        if unavailable():
            return
        accession = accession_loaded.get()
        if not accession:
            say("Fetch a model first, so there is an accession to look up.")
            return
        try:
            say(f"Fetching disease variants for {accession}.")
            sites, site_notes = ann.fetch_variants(accession)
            variant_sites.set(sites)
            if site_notes:
                notes.set(notes.get() + site_notes)
            say(f"{len(sites)} disease variants.")
        except Exception as error:  # noqa: BLE001
            fail("fetch_variants", error)
            return
        LOGGER.info(
            "fetch_variants completed accession=%s sites=%s",
            accession,
            len(sites),
            extra={"event": "fetch_variants"},
        )

    @reactive.effect
    @reactive.event(input.clear_annotations)
    def _clear_annotations() -> None:
        trail.clicked("Clear annotations")
        ptm_sites.set([])
        variant_sites.set([])
        say("Annotations cleared.")

    @reactive.effect
    @reactive.event(input.mode)
    def _mode_hint() -> None:
        # Rewrites the hint only. It must never clear the loaded structure: the
        # ipywidgets version kept both control sets live at once, and switching the
        # radio is not a statement about what is currently loaded.
        if TOPOLOGY_ERROR:
            return
        if input.mode() == "file":
            say("Upload a .pdb or .cif file. No network requests are made in this mode.")
        else:
            say("Enter a UniProtKB accession and fetch a model.")

    @render.text
    def status() -> str:
        return status_text.get()

    @render.ui
    def model_summary():
        if TOPOLOGY_ERROR:
            return ui.p("Topology package unavailable.", class_="scop3p-note")
        bits = summary_bits.get()
        if not bits:
            return ui.p("No structure loaded yet.", class_="scop3p-note")
        return ui.TagList(
            ui.tags.ul(*[ui.tags.li(bit) for bit in bits]),
            ui.p(f"topology {__topology_version__}", class_="scop3p-note"),
        )

    @render.ui
    def notes_panel():
        collected = notes.get()
        if not collected:
            return ui.p(
                "Warnings about missing numbering maps or unavailable sources appear here.",
                class_="scop3p-note",
            )
        return ui.tags.ul(*[ui.tags.li(note) for note in collected])

    @render.ui
    def topology_view():
        if TOPOLOGY_ERROR:
            return scop3p_card(
                "Topology package not found",
                ui.pre(TOPOLOGY_ERROR),
                extra_class="scop3p-status",
            )
        structure_value = structure.get()
        if structure_value is None:
            return ui.p(
                "Fetch a model or upload a structure file to draw a topology.",
                class_="scop3p-note",
            )

        # Defensive read: the chain select may still carry the previous structure's
        # selection for the instant before its update message is applied.
        chain = input.chain()
        if chain not in structure_value.residues_by_chain:
            chain = structure_value.default_chain()

        sites = list(ptm_sites.get()) + list(variant_sites.get())
        show = bool(sites) and loaded_mode.get() == "accession"
        return ui.HTML(
            build_view(
                structure_value,
                chain,
                source.get(),
                height=VIEWER_HEIGHT,
                sites=sites if show else None,
                numbering=numbering.get() if show else None,
                accession=accession_loaded.get(),
                notes=notes.get(),
            )
        )


content_ui = ui.div(
    busy_indicators(),
    app_ui, scop3p_footer()
)

app = App(content_ui, server)
