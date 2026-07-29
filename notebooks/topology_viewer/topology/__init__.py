"""Protein secondary-structure topology viewer."""

__version__ = "1.6.0-phase3"
__build__ = "logo"

from . import annotations
from . import logo
from .annotations import probe as probe_apis
from .app import build_view, diagnose, make_app, save_html
from .elements import annotate_geometry, assign_sheets, build_elements, strand_contacts
from .io import Structure, load_structure, sniff_format
from .layout import build_layout
from .render import build_payload, render, standalone_document

__all__ = [
    "make_app", "build_view", "save_html", "diagnose",
    "load_structure", "sniff_format", "Structure",
    "build_elements", "strand_contacts", "annotate_geometry", "assign_sheets",
    "build_layout", "build_payload", "render", "standalone_document",
    "annotations", "logo", "probe_apis", "__version__", "__build__",
]
