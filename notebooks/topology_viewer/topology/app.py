"""The widget layer.

Two modes that share one renderer:

    accession   a UniProt accession, fetched from AlphaFold DB. Phase 2 adds
                the PDB entry dropdown and the SIFTS numbering map here; phase 3
                hangs PTMs and variants off that mapping.
    file        an uploaded .pdb or .cif. No network, and deliberately no
                UniProt, PTM or variant pipeline: topology and structure only.

The two are kept apart on purpose. The original code guessed an accession from
the uploaded filename and quietly pointed Mol* at AlphaFold, so a locally
predicted model was drawn as a topology while a *different* structure was shown
in 3D beside it, with nothing to signal the mismatch.
"""

from __future__ import annotations

import json
import traceback
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from . import annotations as ann
from . import logo as logo_module
from .elements import annotate_geometry, assign_sheets, build_elements, strand_contacts
from .io import Structure, load_structure
from .layout import build_layout
from .render import build_payload, render

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
USER_AGENT = "topology-viewer/1.0"
LAYOUT_MODES = ("sheet", "serpentine", "spatial")

# Width divided by height of the 2D panel the diagram is drawn into.
PANEL_ASPECT = 1.3


def _fetch(url: str, timeout: int = 60) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def fetch_alphafold(accession: str) -> Tuple[Dict[str, Any], str, str]:
    """Return (metadata, cif text, cif url) for a UniProt accession."""
    accession = accession.strip().upper()
    if not accession:
        raise ValueError("Enter a UniProt accession, for example P07949.")

    try:
        payload = json.loads(_fetch(AFDB_API.format(accession=accession)).decode("utf-8"))
    except urllib.error.HTTPError as error:
        if error.code == 404:
            raise ValueError(
                f"AlphaFold DB has no model for {accession}. Check the accession, "
                "or upload a structure file instead."
            ) from error
        raise ValueError(f"AlphaFold DB returned HTTP {error.code} for {accession}.") from error
    except urllib.error.URLError as error:
        raise ValueError(f"Could not reach AlphaFold DB: {error.reason}") from error

    record = payload[0] if isinstance(payload, list) and payload else payload
    if not isinstance(record, dict):
        raise ValueError(f"AlphaFold DB returned an unexpected response for {accession}.")

    cif_url = record.get("cifUrl")
    if not cif_url:
        entry = record.get("entryId") or f"AF-{accession}-F1"
        version = record.get("latestVersion") or 4
        cif_url = f"https://alphafold.ebi.ac.uk/files/{entry}-model_v{version}.cif"

    return record, _fetch(cif_url).decode("utf-8", errors="replace"), cif_url


def build_view(
    structure: Structure,
    chain: str,
    structure_source: Dict[str, Any],
    height: int = 1150,
    sites: Optional[List[Any]] = None,
    numbering: Optional[Dict[int, int]] = None,
    accession: str = "",
    notes: Optional[List[str]] = None,
) -> str:
    """Run the pipeline for one chain and return renderable HTML."""
    elements, residues = build_elements(structure, chain)
    annotate_geometry(elements, residues)
    contacts = strand_contacts(elements, residues)
    assign_sheets(elements, contacts)

    # The 2D panel is one half of a two-column grid, so the diagram is packed
    # to roughly that shape instead of being scaled down to fit a strip.
    layouts = {
        mode: build_layout(mode, elements, residues, contacts, target_aspect=PANEL_ASPECT)
        for mode in LAYOUT_MODES
    }

    # Sites are attached before the payload is built so each element carries its
    # own marks and counts, which is what the renderer draws from.
    annotation_block = None
    if sites:
        summary = ann.attach_sites(elements, sites, numbering)
        summary["accession"] = accession
        summary["notes"] = list(notes or [])
        summary["numbered"] = numbering is not None
        annotation_block = summary

    payload = build_payload(structure, chain, elements, residues, contacts, layouts)
    payload["structure_source"] = structure_source
    payload["annotations"] = annotation_block
    return render(payload, height=height)


def save_html(
    source: str,
    out_path: str = "topology_view.html",
    chain: Optional[str] = None,
) -> str:
    """Build a view and write it to a standalone HTML file.

    Bypasses Jupyter's output rendering completely, which makes it the fastest
    way to tell a pipeline problem apart from a display problem: if the file
    opens correctly in a browser, the Python side is fine and the issue is how
    the notebook front end is handling the output.

    ``source`` is either a UniProt accession or a path to a .pdb/.cif file.
    """
    from .render import standalone_document

    path = Path(source)
    if path.exists():
        structure = load_structure(path.read_text(encoding="utf-8", errors="replace"), path.name)
        structure_source = {
            "kind": "upload",
            "data": path.read_text(encoding="utf-8", errors="replace"),
            "format": structure.fmt,
        }
    else:
        record, text, cif_url = fetch_alphafold(source)
        structure = load_structure(text, f"AF-{source.upper()}-F1.cif")
        structure.uniprot = source.upper()
        structure_source = {
            "kind": "afdb", "url": cif_url, "format": "mmcif",
            "accession": source.upper(),
        }

    chain = chain or structure.default_chain()
    elements, residues = build_elements(structure, chain)
    annotate_geometry(elements, residues)
    contacts = strand_contacts(elements, residues)
    assign_sheets(elements, contacts)
    layouts = {
        mode: build_layout(mode, elements, residues, contacts, target_aspect=PANEL_ASPECT)
        for mode in LAYOUT_MODES
    }
    payload = build_payload(structure, chain, elements, residues, contacts, layouts)
    payload["structure_source"] = structure_source

    destination = Path(out_path).resolve()
    destination.write_text(standalone_document(payload), encoding="utf-8")
    return str(destination)


def diagnose(source: Optional[str] = None) -> None:
    """Print where the pipeline stands, without involving any widget."""
    import sys

    from . import __build__, __version__

    print("version     :", __version__, "-", __build__)
    print("module      :", __file__)
    print("python      :", sys.version.split()[0])

    try:
        import ipywidgets
        print("ipywidgets  :", ipywidgets.__version__)
    except ImportError:
        print("ipywidgets  : NOT INSTALLED  <- the app cannot render without it")

    try:
        import biotite
        print("biotite     :", biotite.__version__)
    except ImportError:
        print("biotite     : absent (built-in P-SEA will be used)")

    if source is None:
        print("\nPass a file path or accession to test the full pipeline, e.g.")
        print("    diagnose('fixtures/annotated.cif')")
        return

    print(f"\nBuilding {source} ...")
    try:
        written = save_html(source, "topology_diagnostic.html")
    except Exception as error:  # noqa: BLE001 - this function exists to report
        import traceback
        print("FAILED:", type(error).__name__, error)
        traceback.print_exc()
        return
    print("Wrote", written)
    print("Open that file directly in a browser. If it draws there, the Python")
    print("side is working and the problem is how the notebook renders output.")


def _set_options(picker: Any, options: List[Any], value: Any = None) -> None:
    """Replace a Dropdown's options without tripping over the old selection.

    Assigning ``options`` re-validates whatever ``value`` currently holds. If
    that value is absent from the new list -- which it almost always is, since
    the placeholder is "" -- ipywidgets raises
    ``TraitError: Invalid selection: value not found`` from inside the
    assignment, and the widget callback swallows it. Clearing the selection
    first makes the swap safe.
    """
    picker.value = None
    picker.options = options
    if value is not None:
        allowed = [
            item[1] if isinstance(item, tuple) else item for item in options
        ]
        picker.value = value if value in allowed else (allowed[0] if allowed else None)
    elif options:
        first = options[0]
        picker.value = first[1] if isinstance(first, tuple) else first


def make_app(default_accession: str = "P07949", viewer_height: int = 1150) -> Any:
    import ipywidgets as widgets
    from IPython.display import HTML, clear_output, display

    state: Dict[str, Any] = {
        "structure": None,
        "source": {},
        "mode": "accession",
        "suspend_redraw": False,
        "accession": "",
        "refs": [],
        "ptm_sites": [],
        "variant_sites": [],
        "notes": [],
        "numbering": None,
    }

    accession_input = widgets.Text(
        value="",
        placeholder=default_accession,
        description="UniProt",
        layout=widgets.Layout(width="230px"),
    )
    fetch_button = widgets.Button(
        description="Fetch model",
        button_style="primary",
        tooltip="Fetch the AlphaFold DB model for this accession",
        layout=widgets.Layout(width="140px"),
    )
    uploader = widgets.FileUpload(
        accept=".pdb,.ent,.cif,.mmcif",
        multiple=False,
        description="Upload structure",
        tooltip="A predicted or downloaded structure in PDB or mmCIF format",
        layout=widgets.Layout(width="180px"),
    )
    structure_picker = widgets.Dropdown(
        options=[("AlphaFold model", "afdb")],
        description="Structure",
        layout=widgets.Layout(width="430px"),
        style={"description_width": "70px"},
    )
    chain_picker = widgets.Dropdown(
        options=[],
        description="Chain",
        layout=widgets.Layout(width="250px"),
        style={"description_width": "50px"},
    )
    ptm_button = widgets.Button(
        description="Fetch PTMs",
        tooltip="Modifications from Scop3P",
        layout=widgets.Layout(width="130px"),
    )
    uniprot_ptm_toggle = widgets.Checkbox(
        value=True,
        description="Include UniProt PTMs",
        indent=False,
        layout=widgets.Layout(width="190px"),
    )
    variant_button = widgets.Button(
        description="Fetch variants",
        tooltip=(
            "Disease-associated variants from UniProt. Scop3P imports its "
            "mutations from UniProt, so UniProt is the primary source."
        ),
        layout=widgets.Layout(width="140px"),
    )
    clear_button = widgets.Button(
        description="Clear",
        tooltip="Remove all fetched annotations",
        layout=widgets.Layout(width="80px"),
    )
    annotation_bar = widgets.HBox(
        [ptm_button, uniprot_ptm_toggle, variant_button, clear_button],
        layout=widgets.Layout(display="none", flex_flow="row wrap",
                              gap="8px", align_items="center"),
    )
    structure_box = widgets.HBox(
        [structure_picker, chain_picker],
        layout=widgets.Layout(display="none", flex_flow="row wrap", gap="8px"),
    )
    chain_box = structure_box

    status = widgets.HTML()
    output = widgets.Output(layout=widgets.Layout(min_height="640px"))

    def say(message: str, tone: str = "muted") -> None:
        colours = {"muted": "#5f6b7a", "good": "#1f6f43", "warn": "#a15c00", "bad": "#b42318"}
        status.value = (
            f"<span style='color:{colours.get(tone, '#5f6b7a')};font-size:13px'>{message}</span>"
        )

    def guarded(handler):
        """Report failures instead of letting ipywidgets swallow them.

        Exceptions raised inside a widget callback are discarded by ipywidgets:
        no traceback reaches the notebook, and the app simply stops responding.
        Every handler goes through here so a failure is visible.
        """

        def wrapped(*args, **kwargs):
            try:
                return handler(*args, **kwargs)
            except Exception as error:  # noqa: BLE001 - the whole point
                detail = traceback.format_exc()
                say(f"{type(error).__name__}: {error}", "bad")
                with output:
                    clear_output(wait=True)
                    print(detail)

        return wrapped

    def redraw() -> None:
        if state.get("suspend_redraw"):
            return
        structure = state["structure"]
        if structure is None:
            return
        chain = chain_picker.value or structure.default_chain()
        collected = list(state["ptm_sites"]) + list(state["variant_sites"])
        show = bool(collected) and state["mode"] == "accession"
        html = build_view(
            structure,
            chain,
            state["source"],
            height=viewer_height,
            sites=collected if show else None,
            numbering=state["numbering"] if show else None,
            accession=state["accession"],
            notes=state["notes"],
        )
        with output:
            clear_output(wait=True)
            display(HTML(html))

    def adopt(
        structure: Structure,
        source: Dict[str, Any],
        label: str,
        preferred_chain: Optional[str] = None,
    ) -> None:
        state["structure"] = structure
        state["source"] = source

        options = structure.chain_options()
        offer_choice = len(options) > 1 or len(state["refs"]) > 0
        structure_box.layout.display = "flex" if offer_choice else "none"

        # Never call unobserve_all() on a Dropdown. It removes ipywidgets' own
        # internal observer, the one that rebuilds _options_values when options
        # change, so the very next assignment to .value raises
        # "Invalid selection: value not found". Widget callbacks swallow that
        # exception, leaving the app frozen on its last status message with no
        # traceback anywhere. A suppression flag does the same job safely.
        state["suspend_redraw"] = True
        try:
            wanted = preferred_chain if preferred_chain in structure.residues_by_chain else None
            _set_options(chain_picker, options, wanted or structure.default_chain())
        finally:
            state["suspend_redraw"] = False

        bits = [
            label,
            f"{len(structure.chains)} chain(s)",
            f"SS from {structure.ss_source}",
        ]
        if state["ptm_sites"]:
            bits.append(f"{len(state['ptm_sites'])} PTMs")
        if state["variant_sites"]:
            bits.append(f"{len(state['variant_sites'])} variants")
        message = " &middot; ".join(bits)
        for note in state["notes"]:
            message += f"<br><span style='color:#a15c00'>{note}</span>"
        say(message, "good")
        redraw()

    def load_selected_structure() -> None:
        """Load whichever structure the dropdown is pointing at.

        Numbering is resolved here rather than at draw time, because the choice
        of entry is exactly what determines it: AlphaFold models are already
        UniProt-numbered, while a PDB entry needs its SIFTS map or every mark
        lands on the wrong residue.
        """
        accession = state["accession"]
        choice = structure_picker.value

        if choice == "afdb":
            record, text, cif_url = fetch_alphafold(accession)
            structure = load_structure(text, f"AF-{accession}-F1.cif")
            structure.uniprot = accession
            state["numbering"] = None  # identity; positions already match
            adopt(
                structure,
                {"kind": "afdb", "url": cif_url, "format": "mmcif", "accession": accession},
                f"AlphaFold model for {accession}",
            )
            return

        pdb_id = choice
        say(f"Loading PDB entry {pdb_id}.")
        text = ann.fetch_structure_file(pdb_id)
        structure = load_structure(text, f"{pdb_id}.cif")
        structure.uniprot = accession

        # Chains that genuinely map to this accession, per the file's own SIFTS
        # columns, rather than whichever chain happens to be largest.
        mapped_chains = structure.chains_for_accession(accession)
        preferred = (mapped_chains or [None])[0] or _preferred_chain(pdb_id, structure)
        state["mapped_chains"] = mapped_chains
        state["numbering"] = None

        if structure.has_sifts:
            # The updated mmCIF carries the UniProt correspondence alongside the
            # coordinates, so it cannot disagree with what is being drawn, and
            # it costs no extra request.
            mapping = structure.sifts_numbering(preferred or "", accession)
            if mapping:
                state["numbering"] = mapping
                state["numbering_source"] = "SIFTS columns in the mmCIF"

        if state["numbering"] is None:
            try:
                state["numbering"] = ann.fetch_numbering(pdb_id, accession, preferred)
                state["numbering_source"] = "PDBe mapping API"
            except Exception as error:  # noqa: BLE001
                state["notes"] = state["notes"] + [
                    f"No numbering map for {pdb_id} ({error}). Sites are hidden "
                    "rather than drawn at UniProt positions, which would place "
                    "them on the wrong residues."
                ]

        adopt(
            structure,
            {
                "kind": "pdbe",
                "url": ann.structure_file_url(pdb_id),
                "format": "mmcif",
                "accession": accession,
                "pdb_id": pdb_id,
            },
            f"PDB entry {pdb_id}",
            preferred_chain=preferred,
        )

    def _preferred_chain(pdb_id: str, structure: Structure) -> Optional[str]:
        """Pick the chain this accession actually maps to, not just the biggest.

        A complex often has a larger unrelated chain; defaulting to it would show
        a topology for the wrong protein.
        """
        for ref in state["refs"]:
            if ref.pdb_id.upper() == pdb_id.upper():
                for chain in ref.chains:
                    if chain in structure.residues_by_chain:
                        return chain
        return None

    def on_fetch(_: Any = None) -> None:
        accession = (accession_input.value or default_accession).strip().upper()
        state["accession"] = accession
        state["mode"] = "accession"
        state["notes"] = []
        state["numbering"] = None

        say(f"Looking up {accession}.")
        refs: List[Any] = []
        try:
            refs = ann.fetch_structures(accession)
        except Exception as error:  # noqa: BLE001
            state["notes"].append(f"Structure list unavailable ({error}).")
        state["refs"] = refs

        annotation_bar.layout.display = "flex"

        options = [("AlphaFold model (full length)", "afdb")]
        options.extend((ref.label(), ref.pdb_id) for ref in refs)
        state["suspend_redraw"] = True
        try:
            _set_options(structure_picker, options, "afdb")
        finally:
            state["suspend_redraw"] = False

        load_selected_structure()

    def on_fetch_ptms(_: Any = None) -> None:
        accession = state["accession"]
        if not accession:
            say("Fetch a model first, so there is an accession to look up.", "warn")
            return
        say(f"Fetching PTMs for {accession}.")
        sites, notes = ann.fetch_ptms(accession, include_uniprot=uniprot_ptm_toggle.value)
        state["ptm_sites"] = sites
        state["notes"] = notes
        sources = sorted({site.source for site in sites})
        say(
            f"{len(sites)} PTM sites"
            + (f" from {', '.join(sources)}" if sources else "")
            + "".join(f"<br><span style='color:#a15c00'>{note}</span>" for note in notes),
            "good" if sites else "warn",
        )
        redraw()

    def on_fetch_variants(_: Any = None) -> None:
        accession = state["accession"]
        if not accession:
            say("Fetch a model first, so there is an accession to look up.", "warn")
            return
        say(f"Fetching disease variants for {accession}.")
        sites, notes = ann.fetch_variants(accession)
        state["variant_sites"] = sites
        say(
            f"{len(sites)} disease variants"
            + "".join(f"<br><span style='color:#a15c00'>{note}</span>" for note in notes),
            "good" if sites else "warn",
        )
        redraw()

    def on_clear(_: Any = None) -> None:
        state["ptm_sites"] = []
        state["variant_sites"] = []
        say("Annotations cleared.")
        redraw()

    def on_structure_change(_: Any = None) -> None:
        if state.get("suspend_redraw") or not state["accession"]:
            return
        load_selected_structure()

    def on_upload(change: Dict[str, Any]) -> None:
        value = change.get("new")
        if not value:
            return
        item = value[0] if isinstance(value, (list, tuple)) else list(value.values())[0]
        name = item.get("name") or item.get("metadata", {}).get("name") or "uploaded"
        raw = item.get("content")
        if raw is None:
            say("That upload arrived empty. Try selecting the file again.", "warn")
            return
        text = bytes(raw).decode("utf-8", errors="replace")

        say(f"Reading {name}.")
        try:
            structure = load_structure(text, name)
        except Exception as error:  # noqa: BLE001
            say(str(error), "bad")
            return

        # File mode stays offline: the 3D view is fed the uploaded bytes, never
        # a database model that merely shares a name with the file.
        state["mode"] = "file"
        adopt(
            structure,
            {"kind": "upload", "data": text, "format": structure.fmt},
            f"Uploaded {name}",
        )

    fetch_button.on_click(guarded(on_fetch))
    uploader.observe(guarded(on_upload), names="value")
    chain_picker.observe(guarded(lambda change: redraw()), names="value")
    structure_picker.observe(guarded(on_structure_change), names="value")
    ptm_button.on_click(guarded(on_fetch_ptms))
    variant_button.on_click(guarded(on_fetch_variants))
    clear_button.on_click(guarded(on_clear))

    header = widgets.HTML(
        f"""
        <div style="display:flex;align-items:center;gap:14px;
                    font-family:Inter,'Segoe UI',system-ui,sans-serif;margin-bottom:2px">
          <div style="flex:0 0 auto;line-height:0">{logo_module.mark(size=54)}</div>
          <div>
            <div style="font-size:20px;font-weight:650;color:#16202b">
              Protein topology viewer</div>
            <div style="font-size:13px;color:#5f6b7a;max-width:640px">
              Fetch an AlphaFold model by accession, or upload a predicted or downloaded
              structure. Secondary structure is read from the file when present and derived
              from coordinates when it is not.
            </div>
          </div>
        </div>
        """
    )
    controls = widgets.HBox(
        [accession_input, fetch_button, uploader],
        layout=widgets.Layout(flex_flow="row wrap", gap="8px", align_items="center"),
    )

    say("Enter a UniProt accession and fetch a model, or upload a .pdb or .cif file.")
    return widgets.VBox(
        [header, controls, structure_box, annotation_bar, status, output],
        layout=widgets.Layout(gap="10px", width="100%"),
    )
