from __future__ import annotations

from shiny import ui


def scop3p_shell(app_name: str, intro: str, *children: ui.TagChild) -> ui.Tag:
    return ui.page_fluid(
        ui.tags.style(_SCOP3P_CSS),
        ui.div(
            ui.div(
                ui.h1("Scop3P", class_="scop3p-eyebrow"),
                ui.h2(app_name, class_="scop3p-title"),
                ui.p(intro, class_="scop3p-intro"),
                class_="scop3p-hero-copy",
            ),
            class_="scop3p-hero",
        ),
        ui.div(*children, class_="scop3p-shell"),
        scop3p_footer(),
        title=f"Scop3P: {app_name}",
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
                ui.h5("Scop3P", class_="scop3p-footer-head"),
                ui.p(
                    "Protein phosphorylation context across sequence, structure, proteomics, and variant evidence.",
                    class_="scop3p-footer-copy",
                ),
            ),
            ui.div(
                ui.p("Licensed under CC BY 4.0", class_="scop3p-footer-head"),
                ui.p(
                    ui.a(
                        "Data and documentation",
                        href="https://iomics.ugent.be/scop3p/documentation",
                        target="_blank",
                    ),
                    class_="scop3p-footer-link",
                ),
                ui.p(
                    ui.a(
                        "Scop3P API",
                        href="https://iomics.ugent.be/scop3p/api",
                        target="_blank",
                    ),
                    class_="scop3p-footer-link",
                ),
                ui.p(
                    ui.a(
                        "scop3p.compomics@vib-ugent.be",
                        href="mailto:scop3p.compomics@vib-ugent.be",
                    ),
                    class_="scop3p-footer-link",
                ),
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
  padding: 24px 10px 34px;
}
.scop3p-footer-grid {
  max-width: 1600px;
  margin: 0 auto;
  display: grid;
  grid-template-columns: 1.4fr 0.8fr;
  gap: 18px;
  background: rgba(16, 38, 60, 0.96);
  color: #f7fafc;
  border-radius: 18px;
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
