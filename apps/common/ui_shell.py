from __future__ import annotations

import base64
from pathlib import Path

from shiny import ui


_IMAGE_DIR = Path(__file__).resolve().parents[1] / "assets" / "images"
_FOOTER_LOGOS = [
    ("CompOmics", "compomics.png", "60px"),
    ("Bio2Byte", "bio2byte.png", "60px"),
    ("IB2", "IB2.png", "60px"),
    ("VIB Data Core", "vib_Data_Core.png", "120px"),
    ("VIB", "vib.png", "60px"),
    ("UGent", "ugent.png", "60px"),
    ("VUB", "vub.png", "60px"),
    ("ELIXIR Belgium", "elixir-belgium.png", "60px"),
]


#: Multiplier applied to a logo's CSS display height before resizing. 2x keeps it crisp on
#: a high-DPI screen, which is where an under-sized raster is obvious.
_RETINA_SCALE = 2

#: Rendered images are inlined as data URIs, so every byte is paid on every page load --
#: there is no separate cacheable request to amortise them. The source logos are print
#: assets (elixir-belgium.png alone is 382 KB) displayed at 60-120px tall, so they were
#: costing 1.07 MB of HTML per load. Resizing to the height they are actually drawn at
#: removes about 95% of that. Keyed by (filename, target height) and computed once per
#: process, because a resize per page load would trade bandwidth for CPU.
_scaled_cache: dict[tuple[str, int], str] = {}


def _image_data_uri(filename: str) -> str | None:
    """The original file, inlined at full size. The fallback when resizing is unavailable."""
    path = _IMAGE_DIR / filename
    if not path.exists():
        return None
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{payload}"


def _css_pixels(value: str) -> int | None:
    """The number out of a CSS length like ``"60px"``."""
    text = str(value).strip().lower().removesuffix("px").strip()
    try:
        return int(float(text))
    except ValueError:
        return None


def _resized_png(filename: str, *, height: int | None = None, square: int | None = None):
    """PNG bytes for a resized copy, or None when Pillow is unavailable.

    Pillow is a transitive dependency here (it arrives through bokeh) rather than a declared
    one, so every caller has to cope with it being absent.
    """
    path = _IMAGE_DIR / filename
    if not path.exists():
        return None
    try:
        import io

        from PIL import Image

        with Image.open(path) as opened:
            image = opened.convert("RGBA")
            if square is not None:
                image.thumbnail((square, square), Image.LANCZOS)
                canvas = Image.new("RGBA", (square, square), (0, 0, 0, 0))
                canvas.paste(
                    image,
                    ((square - image.width) // 2, (square - image.height) // 2),
                )
                image = canvas
            elif height is not None and image.height > height:
                width = max(1, round(image.width * height / image.height))
                image = image.resize((width, height), Image.LANCZOS)
            buffer = io.BytesIO()
            image.save(buffer, format="PNG", optimize=True)
        return buffer.getvalue()
    except Exception:  # noqa: BLE001 - a logo must never break a page
        return None


def _scaled_image_data_uri(filename: str, display_height: str) -> str | None:
    """A logo inlined at the size it is actually displayed, not its source size."""
    target = _css_pixels(display_height)
    if target is None:
        return _image_data_uri(filename)

    target *= _RETINA_SCALE
    cache_key = (filename, target)
    cached = _scaled_cache.get(cache_key)
    if cached is not None:
        return cached or None

    payload = _resized_png(filename, height=target)
    if payload is None:
        fallback = _image_data_uri(filename)
        _scaled_cache[cache_key] = fallback or ""
        return fallback

    uri = f"data:image/png;base64,{base64.b64encode(payload).decode('ascii')}"
    _scaled_cache[cache_key] = uri
    return uri


#: The browser asks for /favicon.ico on every page and nothing answered, so each load
#: logged a 404. Declaring the icon in the head stops the request being made at all,
#: which is tidier than adding a route to every app plus the portal.
FAVICON_SOURCE = "scop3p.png"

#: Pixels per side for the generated icon. Browsers display a favicon at 16-32px, so the
#: 82 KB source is far larger than needed; 64 covers high-DPI tabs.
FAVICON_SIZE = 64

_favicon_cache: str | None = None


def _favicon_data_uri() -> str | None:
    """A small square PNG icon derived from the Scop3P logo, as a data URI.

    Squared rather than merely shrunk: the logo is wider than it is tall, and a browser fits
    a favicon into a square tab box, which would scale a wide image down to the height of
    its narrower side.

    Falls back to the full-size image if resizing is unavailable -- wasteful, but a working
    page, and a cosmetic icon must never take an app down.
    """
    global _favicon_cache
    if _favicon_cache is not None:
        return _favicon_cache or None

    payload = _resized_png(FAVICON_SOURCE, square=FAVICON_SIZE)
    if payload is None:
        fallback = _image_data_uri(FAVICON_SOURCE)
        _favicon_cache = fallback or ""
        return fallback

    _favicon_cache = f"data:image/png;base64,{base64.b64encode(payload).decode('ascii')}"
    return _favicon_cache


def favicon_tags() -> list[ui.Tag]:
    """Head tags declaring the tab icon, or nothing if the asset is missing."""
    href = _favicon_data_uri()
    if not href:
        return []
    return [
        ui.tags.link(rel="icon", type="image/png", href=href),
        # Safari and iOS look for this one specifically.
        ui.tags.link(rel="apple-touch-icon", href=href),
    ]


def _footer_logo_tags() -> list[ui.Tag]:
    tags: list[ui.Tag] = []
    for label, filename, height in _FOOTER_LOGOS:
        src = _scaled_image_data_uri(filename, height)
        if src is None:
            tags.append(ui.span(label, class_="scop3p-logo-fallback", title=f"Missing asset: {filename}"))
            continue
        tags.append(
            ui.tags.img(
                src=src,
                alt=label,
                title=label,
                class_="scop3p-footer-logo",
                style=f"height: {height}",
            )
        )
    return tags


#: The one name for a UniProt identifier field, used by every app. Previously each app
#: invented its own -- "ACC_ID (UniProt accession number)", "UniProt",
#: "UniProt accession (AlphaFold DB / PDBe)" -- which made the toolkit read as five
#: unrelated tools. Change it here and every app follows.
ACCESSION_LABEL = "UniProtKB accession"


def scop3p_structure_picker(input_id: str, label: str, choices: dict[str, str]):
    """A searchable single-select for choosing a structure.

    Searchable rather than a plain ``<select>`` because the option count follows the
    protein, not the interface: P04637 (p53) yields 629 selectable chains, and scanning
    that many entries by eye is not a real option. Typing "2IVT" or "1.0 A" narrows it
    immediately.

    Selectize is configured to behave like a picker rather than a tag editor: no free
    text, and no remove button, because clearing the only selection would leave the app
    with no structure at all. Every caller must update it with ``ui.update_selectize`` --
    ``ui.update_select`` targets a different client-side widget and fails silently, which
    looks exactly like an upstream returning nothing.
    """
    return ui.input_selectize(
        input_id,
        label,
        choices=choices,
        multiple=False,
        remove_button=False,
        options={"placeholder": "Type to filter structures"},
    )


def scop3p_field_row(*children: ui.TagChild, extra_class: str = "") -> ui.Tag:
    """Put an input and its buttons on one baseline.

    Shiny renders a text input as a label above a control with a margin below it,
    while an action button has no label. Dropped into adjacent grid columns the
    button therefore floats above the input it belongs to. This row pins every
    child to the bottom edge and drops the trailing margin, so controls line up.
    """
    classes = "scop3p-field-row"
    if extra_class:
        classes = f"{classes} {extra_class}"
    return ui.div(*children, class_=classes)


def scop3p_example_button(input_id: str, label: str = "Load example") -> ui.Tag:
    """A quiet button that fills the adjacent input with a worked example.

    Every accession field documents an example in its placeholder; this makes that
    example usable in one click instead of something to retype.
    """
    return ui.input_action_button(input_id, label, class_="scop3p-example-btn")


def scop3p_shell(app_name: str, intro: str, *children: ui.TagChild) -> ui.Tag:
    return ui.page_fluid(
        ui.head_content(*favicon_tags()),
        ui.tags.style(_SCOP3P_CSS),
        ui.div(
            ui.div(
                ui.h1("Scop3P-Toolkit", class_="scop3p-eyebrow"),
                ui.h2(app_name, class_="scop3p-title"),
                ui.p(intro, class_="scop3p-intro"),
                class_="scop3p-hero-copy",
            ),
            class_="scop3p-hero",
        ),
        ui.div(*children, class_="scop3p-shell"),
        title=f"Scop3P-Toolkit: {app_name}",
        style="min-height: 100dvh;"
    )


def scop3p_card(title: str, *children: ui.TagChild, extra_class: str = "") -> ui.Tag:
    classes = "scop3p-card"
    if extra_class:
        classes = f"{classes} {extra_class}"
    return ui.div(
        ui.h4(title, class_="scop3p-card-title"),
        *children,
        class_=classes,
    )


#: Footer disclaimer. These protocols are thin clients over other people's services --
#: UniProt, Scop3P, PDBe, the EBI Proteins API, RCSB and AlphaFold DB -- and those services
#: do fail transiently: dropped TLS handshakes, connections closed mid-response, truncated
#: JSON bodies. The app retries once and caches what it gets, but it cannot make an
#: unavailable service available, so the honest thing is to tell the user that retrying is
#: worth doing rather than let a network error read as a mistake on their part.
EXTERNAL_RESOURCES_NOTICE = (
    "These protocols query external online resources (UniProt, Scop3P, PDBe, the EBI "
    "Proteins API, RCSB PDB and AlphaFold DB) as you use them, so results depend on those "
    "services being reachable. If something fails unexpectedly, it is usually a temporary "
    "problem at the source rather than anything wrong with your input \u2014 please wait a "
    "moment and try the action again."
)


#: The preprint describing the toolkit. Kept as fields rather than one pre-formatted string
#: so the DOI can be linked, and so a test can check the rendered footer against them.
#: Transcribed from the BibTeX record, with its LaTeX escapes resolved (D{\'\i}az -> Díaz,
#: Adri{\'a}n -> Adrián).
CITATION = {
    "authors": (
        "Díaz A, Tichshenko N, Depoortere B, Andrade Buono R, De Geest P, "
        "Vranken WF, Martens L, Ramasamy P"
    ),
    "title": (
        "Scop3P-Toolkit: executable structure-aware workflows linking PTMs, "
        "peptides, and mutations to protein function"
    ),
    "venue": "bioRxiv",
    "year": "2026",
    "doi": "10.64898/2026.08.04.742789",
}

#: doi.org rather than the biorxiv URL in the record: a DOI keeps resolving if the preprint
#: is published in a journal, where the "early" biorxiv path would not.
CITATION_DOI_URL = f"https://doi.org/{CITATION['doi']}"


def citation_tags() -> list[ui.Tag]:
    """The "please cite" block, for the footer.

    Placed above the affiliation logos: whoever is looking for how to cite the work should
    reach it before the institutional marks, not after them.
    """
    return [
        ui.p("If you use Scop3P-Toolkit, please cite:", class_="scop3p-footer-head"),
        ui.p(
            f"{CITATION['authors']}. ",
            ui.em(f"{CITATION['title']}. "),
            f"{CITATION['venue']} ({CITATION['year']}). ",
            # target="_blank" with rel="noopener": without it the opened page can reach
            # back through window.opener.
            ui.a(
                f"doi:{CITATION['doi']}",
                href=CITATION_DOI_URL,
                target="_blank",
                rel="noopener noreferrer",
            ),
            class_="scop3p-footer-copy scop3p-footer-citation",
        ),
    ]


def scop3p_footer() -> ui.Tag:
    return ui.tags.footer(
        ui.div(
            ui.div(
                ui.h5("Scop3P-Toolkit", class_="scop3p-footer-head"),
                ui.p(
                    "Protein phosphorylation context across sequence, structure, proteomics, and variant evidence.",
                    class_="scop3p-footer-copy",
                ),
                # Every protocol reads live data from UniProt, Scop3P, PDBe, the EBI
                # Proteins API, RCSB and AlphaFold DB, so a failure here is usually an
                # upstream hiccup rather than bad input. Saying so in the footer means the
                # explanation is on screen when it is needed, instead of the user
                # concluding their accession is wrong.
                ui.p(
                    EXTERNAL_RESOURCES_NOTICE,
                    class_="scop3p-footer-copy scop3p-footer-notice",
                ),
                *citation_tags(),
                ui.div(*_footer_logo_tags(), class_="scop3p-footer-logos"),
                ui.p(
                    "Licensed under Apache 2.0.", 
                    class_="scop3p-footer-head"
                ),
            ),
            ui.div(
                ui.p(
                    "You can find data and documentation (in PDF) ",
                    ui.a(
                        "here.",
                        href="https://iomics.ugent.be/scop3p/documentation",
                        target="_blank",
                    ),
                    class_="scop3p-footer-link",
                ),
                ui.p(
                    "Official documentation of the Scop3P REST API: ",
                    ui.a(
                        "Swagger OpenAPI definition.",
                        href="https://iomics.ugent.be/scop3p/swagger-ui/index.html",
                        target="_blank",
                    ),
                    class_="scop3p-footer-link",
                ),
                ui.p(
                    "For any further questions, feedback or suggestions, please send an email to: ",
                    ui.a(
                        "scop3p.compomics@vib-ugent.be",
                        href="mailto:scop3p.compomics@vib-ugent.be",
                    ),
                    ", or ",
                    ui.a(
                        "bio2byte@vub.be",
                        href="mailto:bio2byte@vub.be",
                    ),
                    ".",
                    class_="scop3p-footer-link",
                ),
                ui.p(
                    "Scop3P is part of the ELIXIR Belgium infrastructure as a Core Service since 2020: ",
                    ui.a(
                        "Learn more",
                        href="https://www.elixir-belgium.org/services/scop3p/"
                    ),
                    ".",
                    class_="scop3p-footer-link",
                ),
                ui.p(
                    "Scop3P Toolkit has been developed in Belgium by CompOmics (UGENT) in collaboration with the Bio2Byte Lab (VUB) and the VIB Data Core team (VIB).",
                    class_="scop3p-footer-head"
                )
            ),
            class_="scop3p-footer-grid",
        ),
        class_="scop3p-footer",
    )


_SCOP3P_CSS = """
:root {
  --scop3p-ink: #1f2937;
  --scop3p-muted: #5f7084;
  --scop3p-line: #d8dee8;
  --scop3p-panel: rgba(255, 255, 255, 0.92);
  --scop3p-surface: #f4f1ea;
  --scop3p-accent: #1b5f94;
  --scop3p-accent-2: #d48a2c;
  --scop3p-shadow: 0 18px 42px rgba(19, 38, 58, 0.12);
}
body {
  background:
    radial-gradient(circle at top right, rgba(212, 138, 44, 0.12), transparent 24%),
    radial-gradient(circle at top left, rgba(27, 95, 148, 0.10), transparent 28%),
    linear-gradient(180deg, #f8f5ef 0%, #eef3f7 100%);
  color: var(--scop3p-ink);
}
.scop3p-shell {
  max-width: 1600px;
  margin: 0 auto;
  padding: 0 10px 32px;
}
.scop3p-hero {
  max-width: 1600px;
  margin: 0 auto;
  padding: 20px 10px 18px;
}
.scop3p-hero-copy {
  background: linear-gradient(135deg, #10263c 0%, #21547b 58%, #2f6a95 100%);
  color: #fff;
  border-radius: 20px;
  padding: 24px 26px;
  box-shadow: 0 22px 46px rgba(16, 38, 60, 0.22);
}
.scop3p-eyebrow {
  margin: 0 0 8px;
  font-size: 0.9rem;
  letter-spacing: 0.18em;
  text-transform: uppercase;
  color: rgba(255,255,255,0.68);
}
.scop3p-title {
  margin: 0;
  font-size: 2rem;
  font-weight: 700;
  letter-spacing: -0.03em;
}
.scop3p-intro {
  margin: 10px 0 0;
  max-width: 75ch;
  color: rgba(255,255,255,0.86);
}
.scop3p-card {
  background: var(--scop3p-panel);
  border: 1px solid var(--scop3p-line);
  border-radius: 16px;
  padding: 16px;
  box-shadow: var(--scop3p-shadow);
  backdrop-filter: blur(4px);
}
.scop3p-card-title {
  margin: 0 0 12px;
  font-size: 1.1rem;
}
.scop3p-status pre {
  margin: 0;
  white-space: pre-wrap;
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 13px;
}
.scop3p-note {
  color: var(--scop3p-muted);
  font-size: 0.95rem;
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
/* Columns stretch to a common height, so a results card is never a short box
   floating beside a tall controls card. Grid items stretch by default; the point is
   that we do NOT set align-items:start here. */
.scop3p-two-col {
  display: grid;
  grid-template-columns: 420px minmax(0, 1fr);
  gap: 18px;
  align-items: stretch;
}
/* A card that is a grid or flex child fills the space it is given. */
.scop3p-two-col > .scop3p-card,
.scop3p-header-grid > .scop3p-card {
  height: 100%;
}
.scop3p-field-row {
  display: flex;
  flex-wrap: wrap;
  align-items: flex-end;
  gap: 10px;
  margin-bottom: 14px;
}
.scop3p-field-row > .shiny-input-container {
  flex: 1 1 200px;
  margin-bottom: 0;
}
/* Buttons in a field row size to their label and stay on the input's baseline.
   width/flex are pinned because app stylesheets legitimately set `.btn { width: 100% }`
   for their own button grids, and that would otherwise stretch these across the row. */
.scop3p-field-row > .btn,
.scop3p-field-row .shiny-input-container > .btn {
  margin-bottom: 0;
  white-space: nowrap;
  width: auto;
  flex: 0 0 auto;
}
.scop3p-example-btn {
  background: transparent;
  border: 1px dashed var(--scop3p-line);
  color: var(--scop3p-accent);
  font-size: 0.85rem;
  font-weight: 600;
  padding: 7px 12px;
  border-radius: 9px;
  white-space: nowrap;
}
.scop3p-example-btn:hover,
.scop3p-example-btn:focus {
  border-style: solid;
  border-color: var(--scop3p-accent);
  background: rgba(27, 95, 148, 0.07);
  color: var(--scop3p-accent);
}
.scop3p-header-grid {
  display: grid;
  grid-template-columns: 1.4fr 0.8fr;
  gap: 18px;
  margin-bottom: 18px;
}
.nav-tabs {
  margin-top: 10px;
  border-bottom-color: var(--scop3p-line);
}
.nav-tabs .nav-link {
  color: #274560;
  font-weight: 600;
}
.nav-tabs .nav-link.active {
  background: rgba(255,255,255,0.9);
  border-color: var(--scop3p-line) var(--scop3p-line) transparent;
}
.shiny-input-container {
  width: 100%;
}
.btn-primary {
  background: linear-gradient(135deg, #1f6fb2 0%, #3696c9 100%);
  border-color: #1f6fb2;
}
.btn-warning {
  background: linear-gradient(135deg, #d48a2c 0%, #ebb64d 100%);
  border-color: #d48a2c;
  color: #1f2328;
}
.btn-info {
  background: linear-gradient(135deg, #5ca9c7 0%, #72c4dc 100%);
  border-color: #5ca9c7;
  color: #13222d;
}
.btn-success {
  background: linear-gradient(135deg, #3f8d6a 0%, #56ad82 100%);
  border-color: #3f8d6a;
}
.btn-danger {
  background: linear-gradient(135deg, #b24a4a 0%, #d66c55 100%);
  border-color: #b24a4a;
}
.scop3p-footer {
  margin-top: 28px;
  width: 100vw;
  margin-left: calc(50% - 50vw);
  margin-right: calc(50% - 50vw);
  padding: 24px 10px 0px 10px;
}
.scop3p-footer-grid {
  max-width: 100%;
  display: grid;
  grid-template-columns: 1.4fr 0.8fr;
  gap: 18px;
  background: rgba(16, 38, 60, 0.96);
  color: #f7fafc;
  padding: 20px 22px;
  box-shadow: 0 20px 42px rgba(16, 38, 60, 0.24);
}
.scop3p-footer-head {
  margin: 0 0 8px;
  font-weight: 700;
}
.scop3p-footer-copy,
.scop3p-footer-link {
  margin: 0 0 8px;
  color: rgba(247,250,252,0.82);
}
.scop3p-footer-notice {
  font-size: 0.9em;
  line-height: 1.5;
  color: rgba(247,250,252,0.68);
}
.scop3p-footer-citation {
  font-size: 0.9em;
  line-height: 1.5;
  margin-bottom: 14px;
}
.scop3p-footer-logos {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 14px;
  margin-top: 12px;
  margin-bottom: 24px;
}
.scop3p-footer-logo {
  height: 60px;
  width: auto;
  max-width: 150px;
  object-fit: contain;
  opacity: 0.92;
}
.scop3p-logo-fallback {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  min-height: 34px;
  padding: 0 12px;
  border: 1px solid rgba(255,255,255,0.18);
  border-radius: 999px;
  color: rgba(247,250,252,0.82);
  font-size: 0.85rem;
  font-weight: 600;
  letter-spacing: 0.02em;
}
.scop3p-footer a {
  color: #87d7ff;
  text-decoration: none;
}
.scop3p-footer a:hover {
  text-decoration: underline;
}
@media (max-width: 1000px) {
  .scop3p-header-grid,
  .scop3p-two-col,
  .scop3p-footer-grid {
    grid-template-columns: 1fr;
  }
}
"""
