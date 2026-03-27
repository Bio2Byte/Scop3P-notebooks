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


def _image_data_uri(filename: str) -> str | None:
    path = _IMAGE_DIR / filename
    if not path.exists():
        return None
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{payload}"


def _footer_logo_tags() -> list[ui.Tag]:
    tags: list[ui.Tag] = []
    for label, filename, height in _FOOTER_LOGOS:
        src = _image_data_uri(filename)
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


def scop3p_shell(app_name: str, intro: str, *children: ui.TagChild) -> ui.Tag:
    return ui.page_fluid(
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


def scop3p_footer() -> ui.Tag:
    return ui.tags.footer(
        ui.div(
            ui.div(
                ui.h5("Scop3P-Toolkit", class_="scop3p-footer-head"),
                ui.p(
                    "Protein phosphorylation context across sequence, structure, proteomics, and variant evidence.",
                    class_="scop3p-footer-copy",
                ),
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
.scop3p-two-col {
  display: grid;
  grid-template-columns: 420px minmax(0, 1fr);
  gap: 18px;
  align-items: start;
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
