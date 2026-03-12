from __future__ import annotations

from shiny import App
from starlette.testclient import TestClient

from mutation_effect.app import app as mutation_effect_app
from peptide_mapper.app import app as peptide_mapper_app
from portal.main import app as portal_app
from structure_viz.app import app as structure_viz_app


def test_peptide_mapper_app_constructs() -> None:
    assert isinstance(peptide_mapper_app, App)


def test_structure_viz_app_constructs() -> None:
    assert isinstance(structure_viz_app, App)


def test_mutation_effect_app_constructs() -> None:
    assert isinstance(mutation_effect_app, App)


def test_portal_app_constructs() -> None:
    assert callable(portal_app)


def test_portal_root_selector_and_cookie() -> None:
    client = TestClient(portal_app)

    default_response = client.get("/")
    assert default_response.status_code == 200
    assert "App selector" in default_response.text
    assert "Peptide Mapper" in default_response.text
    assert "scop3p_app=peptide-mapper" in default_response.headers["set-cookie"]

    selected_response = client.get("/?app=structure-viz")
    assert selected_response.status_code == 200
    assert "Structure Visualisation" in selected_response.text
    assert "scop3p_app=structure-viz" in selected_response.headers["set-cookie"]
