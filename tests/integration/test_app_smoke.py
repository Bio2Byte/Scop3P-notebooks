from __future__ import annotations

from pathlib import Path

from common.logging_utils import configure_logging, get_log_file_path, get_metadata_path
import pytest
from shiny import App
from starlette.testclient import TestClient

from help.app import PROTOCOLS, app as help_app
from mutation_effect.app import app as mutation_effect_app
from peptide_mapper.app import app as peptide_mapper_app
from portal.main import APP_OPTIONS, app as portal_app
from rinalign.app import app as rinalign_app
from structure_viz.app import app as structure_viz_app
from topology_viewer.app import app as topology_viewer_app


def test_peptide_mapper_app_constructs() -> None:
    assert isinstance(peptide_mapper_app, App)


def test_structure_viz_app_constructs() -> None:
    assert isinstance(structure_viz_app, App)


def test_mutation_effect_app_constructs() -> None:
    assert isinstance(mutation_effect_app, App)


def test_topology_viewer_app_constructs() -> None:
    assert isinstance(topology_viewer_app, App)


def test_rinalign_app_constructs() -> None:
    assert isinstance(rinalign_app, App)


def test_help_app_constructs() -> None:
    assert isinstance(help_app, App)


def test_help_documents_every_protocol_app() -> None:
    """The Help page has to keep up with the navbar.

    Every app in APP_OPTIONS except Help itself must have a card, and every Help card
    must point at a real app key -- otherwise Help either hides a tool or offers a
    dead link.
    """
    documented = {protocol.key for protocol in PROTOCOLS}
    registered = set(APP_OPTIONS) - {"help"}
    assert documented == registered, (
        f"undocumented: {sorted(registered - documented)}, "
        f"stale in Help: {sorted(documented - registered)}"
    )


def test_help_cards_are_substantive() -> None:
    """Guards against a protocol being added with placeholder text."""
    for protocol in PROTOCOLS:
        assert protocol.question.endswith("?"), protocol.key
        assert len(protocol.mission) > 80, protocol.key
        assert len(protocol.scope) >= 3, protocol.key
        assert len(protocol.use_cases) >= 3, protocol.key
        assert protocol.spec.startswith("docs/use-cases/"), protocol.key


def test_help_spec_links_point_at_real_files() -> None:
    root = Path(__file__).resolve().parents[2]
    for protocol in PROTOCOLS:
        assert (root / protocol.spec).is_file(), f"{protocol.key}: missing {protocol.spec}"


def test_portal_app_constructs() -> None:
    assert callable(portal_app)


def test_portal_root_selector_and_cookie() -> None:
    client = TestClient(portal_app)

    default_response = client.get("/")
    assert default_response.status_code == 200
    assert "Scop3P-Toolkit" in default_response.text
    assert "Peptide Mapper" in default_response.text
    assert "scop3p_app=peptide-mapper" in default_response.headers["set-cookie"]

    selected_response = client.get("/?app=structure-viz")
    assert selected_response.status_code == 200
    assert "Structure Visualisation" in selected_response.text
    assert "scop3p_app=structure-viz" in selected_response.headers["set-cookie"]

    # Every registered app must appear in the injected navbar, whichever app is
    # being served. This is what catches an APP_OPTIONS entry that was forgotten.
    for label, _icon, _app in APP_OPTIONS.values():
        assert label in default_response.text, f"{label} is missing from the navbar"


@pytest.mark.parametrize("key", sorted(APP_OPTIONS))
def test_portal_serves_every_registered_app(key: str) -> None:
    label = APP_OPTIONS[key][0]
    client = TestClient(portal_app)
    response = client.get(f"/?app={key}")
    assert response.status_code == 200
    assert label in response.text
    assert f"scop3p_app={key}" in response.headers["set-cookie"]


def test_new_apps_are_registered_in_the_portal() -> None:
    """They are unreachable in the published all-in-one image otherwise."""
    assert APP_OPTIONS["topology-viewer"][2] is topology_viewer_app
    assert APP_OPTIONS["rinalign"][2] is rinalign_app


def test_topology_bridge_resolves_in_a_checkout() -> None:
    """The Topology Viewer app is useless if the bridge cannot find the package.

    Fails loudly if notebooks/topology_viewer moves, or if apps/common is relocated
    and the bridge's parents[2] arithmetic stops pointing at the repository root.
    """
    from common.topology_bridge import TOPOLOGY_ERROR, __topology_version__, build_view

    assert TOPOLOGY_ERROR is None, TOPOLOGY_ERROR
    assert callable(build_view)
    assert __topology_version__


def test_portal_logging_metadata_configured() -> None:
    configure_logging()
    assert get_log_file_path() is not None
    metadata_path = get_metadata_path()
    assert metadata_path is not None
    assert metadata_path.exists()


def test_every_app_footer_carries_the_external_resources_notice() -> None:
    """The protocols are thin clients over other people's services.

    Those services fail transiently -- dropped handshakes, truncated bodies -- and without
    saying so a network error reads as "my accession must be wrong". The notice belongs on
    every protocol, not just the one that happened to get it, so it is asserted per app.
    """
    from common.ui_shell import EXTERNAL_RESOURCES_NOTICE

    # A distinctive fragment, so rewording the sentence does not silently drop the check.
    fragment = "external online resources"
    assert fragment in EXTERNAL_RESOURCES_NOTICE

    for key, (_label, _icon, app) in APP_OPTIONS.items():
        rendered = str(app.ui)
        assert fragment in rendered, f"{key} does not show the external-resources notice"
        assert "try the action again" in rendered, f"{key} does not tell the user to retry"


def test_the_notice_names_the_services_it_depends_on() -> None:
    """Naming them lets a user check a status page rather than guess."""
    from common.ui_shell import EXTERNAL_RESOURCES_NOTICE

    for service in ("UniProt", "Scop3P", "PDBe", "RCSB", "AlphaFold"):
        assert service in EXTERNAL_RESOURCES_NOTICE, f"{service} is not named"


def test_the_notice_follows_the_tagline() -> None:
    """Placement was requested specifically: after the one-line description."""
    from common.ui_shell import EXTERNAL_RESOURCES_NOTICE, scop3p_footer

    rendered = str(scop3p_footer())
    tagline = "Protein phosphorylation context across sequence"
    assert tagline in rendered
    assert rendered.index(tagline) < rendered.index("external online resources")
    assert EXTERNAL_RESOURCES_NOTICE.split(" — ")[0][:40] in rendered.replace("&#8212;", "—")


def test_every_app_declares_a_favicon() -> None:
    """Without it the browser asks for /favicon.ico on every load and gets a 404.

    Declaring the icon stops the request being made at all, which beats adding a route to
    six apps plus the portal. Asserted per app so one shell change cannot leave some behind.
    """
    for key, (_label, _icon, app) in APP_OPTIONS.items():
        rendered = str(app.ui)
        assert 'rel="icon"' in rendered, f"{key} declares no favicon"
        assert "data:image/png;base64," in rendered, f"{key}'s favicon has no payload"


def test_the_favicon_is_small_enough_to_inline_on_every_page() -> None:
    """It is a data URI, so its size is paid on every page load.

    The source logo is 82 KB, which would be ~110 KB of base64 per load for something a
    browser draws at 16-32px.
    """
    from common.ui_shell import _favicon_data_uri

    uri = _favicon_data_uri()
    assert uri is not None
    assert len(uri) < 20_000, f"favicon data URI is {len(uri)/1024:.0f} KB; downscale it"


def test_the_favicon_is_square() -> None:
    """A browser fits the icon to a square box, so a wide logo ends up tiny."""
    import base64
    import io

    from common.ui_shell import FAVICON_SIZE, _favicon_data_uri

    PIL = pytest.importorskip("PIL.Image")
    raw = base64.b64decode(_favicon_data_uri().split(",", 1)[1])
    with PIL.open(io.BytesIO(raw)) as image:
        assert image.size == (FAVICON_SIZE, FAVICON_SIZE)


def test_a_missing_asset_does_not_break_the_page(monkeypatch, tmp_path) -> None:
    """A cosmetic icon must never take an app down."""
    from common import ui_shell

    monkeypatch.setattr(ui_shell, "_IMAGE_DIR", tmp_path)
    monkeypatch.setattr(ui_shell, "_favicon_cache", None)
    assert ui_shell.favicon_tags() == []
    assert str(ui_shell.scop3p_shell("X", "y"))  # still renders


def test_pillow_being_absent_does_not_break_the_page(monkeypatch) -> None:
    """Pillow is transitive (via bokeh), not declared, so it may not be there.

    The fallback is the full-size image: wasteful, but a working page.
    """
    import builtins

    from common import ui_shell

    real_import = builtins.__import__

    def no_pil(name, *args, **kwargs):
        if name.startswith("PIL"):
            raise ImportError("no Pillow")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_pil)
    monkeypatch.setattr(ui_shell, "_favicon_cache", None)
    uri = ui_shell._favicon_data_uri()
    monkeypatch.undo()
    monkeypatch.setattr(ui_shell, "_favicon_cache", None)

    assert uri is not None and uri.startswith("data:image/png;base64,")


def test_every_app_shows_the_citation() -> None:
    """A published tool has to tell people how to cite it, from wherever they are."""
    from common.ui_shell import CITATION

    for key, (_label, _icon, app) in APP_OPTIONS.items():
        rendered = str(app.ui)
        assert "please cite" in rendered, f"{key} shows no citation"
        assert CITATION["doi"] in rendered, f"{key} shows no DOI"


def test_the_citation_carries_every_field_from_the_record() -> None:
    from common.ui_shell import CITATION, scop3p_footer

    rendered = str(scop3p_footer())
    for field in ("title", "venue", "year", "doi"):
        assert CITATION[field] in rendered, f"the footer omits the {field}"
    # First and last author, so a truncated list is caught.
    assert "Díaz A" in rendered
    assert "Ramasamy P" in rendered


def test_the_doi_is_a_link_that_opens_in_a_new_tab() -> None:
    """Requested explicitly, and it is the right default: losing the app to a navigation
    would discard whatever the user had loaded."""
    from common.ui_shell import CITATION_DOI_URL, scop3p_footer

    rendered = str(scop3p_footer())
    assert f'href="{CITATION_DOI_URL}"' in rendered
    assert 'target="_blank"' in rendered
    assert "noopener" in rendered, (
        "target=_blank without rel=noopener lets the opened page reach window.opener"
    )


def test_the_doi_link_uses_doi_org_not_a_versioned_preprint_path() -> None:
    """A DOI keeps resolving if the preprint is published; an "early" biorxiv path may not."""
    from common.ui_shell import CITATION_DOI_URL

    assert CITATION_DOI_URL.startswith("https://doi.org/")
    assert "biorxiv.org" not in CITATION_DOI_URL


def test_the_citation_sits_above_the_affiliation_logos() -> None:
    """Where it was asked for: someone looking for how to cite should reach it before the
    institutional marks."""
    from common.ui_shell import scop3p_footer

    rendered = str(scop3p_footer())
    assert rendered.index("please cite") < rendered.index("scop3p-footer-logos")


def test_the_latex_escapes_from_the_bibtex_are_resolved() -> None:
    """The record spells the first author D{\\'\\i}az; a footer must not."""
    from common.ui_shell import CITATION

    joined = " ".join(CITATION.values())
    for artefact in ("{", "}", "\\'", "\\`", "\\textendash"):
        assert artefact not in joined, f"unresolved LaTeX {artefact!r} in the citation"


def test_the_navbar_links_to_the_preprint() -> None:
    """The navbar is the one thing on screen in every protocol, so it is where a reader
    who has not scrolled to the footer will look."""
    from common.ui_shell import CITATION_DOI_URL

    with TestClient(portal_app) as client:
        body = client.get("/?app=structure-viz").text

    assert "read pre-print" in body, "the navbar does not offer the preprint"
    assert f'href="{CITATION_DOI_URL}"' in body, "the navbar link is not the DOI"
    assert "Tools for exploring and extending Scop3P" in body, "the subtitle was lost"


def test_the_navbar_preprint_link_opens_in_a_new_tab() -> None:
    """Navigating away would discard whatever the user had loaded in the app."""
    with TestClient(portal_app) as client:
        body = client.get("/").text

    marker = body.index("read pre-print")
    anchor = body.rfind("<a", 0, marker)
    tag = body[anchor:marker]
    assert 'target="_blank"' in tag
    assert "noopener" in tag, "target=_blank without rel=noopener exposes window.opener"


def test_the_navbar_and_the_footer_cite_the_same_doi() -> None:
    """Two links to the same paper must not drift apart."""
    from common.ui_shell import CITATION, CITATION_DOI_URL

    with TestClient(portal_app) as client:
        body = client.get("/").text

    assert body.count(CITATION_DOI_URL) >= 2, "navbar and footer do not share the DOI URL"
    assert CITATION["doi"] in body
