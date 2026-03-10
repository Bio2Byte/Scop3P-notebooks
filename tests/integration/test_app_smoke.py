from __future__ import annotations

from shiny import App

from apps.peptide_mapper.app import app as peptide_mapper_app
from apps.structure_viz.app import app as structure_viz_app


def test_peptide_mapper_app_constructs() -> None:
    assert isinstance(peptide_mapper_app, App)


def test_structure_viz_app_constructs() -> None:
    assert isinstance(structure_viz_app, App)
