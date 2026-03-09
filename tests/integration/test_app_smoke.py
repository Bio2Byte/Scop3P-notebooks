from __future__ import annotations

from shiny import App

from apps.peptide_mapper.app import app


def test_peptide_mapper_app_constructs() -> None:
    assert isinstance(app, App)
