import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from shiny import App, ui

app_ui = ui.page_fluid(
    ui.h2("Structure Visualisation (Pending Conversion)"),
    ui.p("This app entrypoint is scaffolded. Conversion from the Voilà notebook is planned next."),
)


def server(input, output, session):
    return None


app = App(app_ui, server)
