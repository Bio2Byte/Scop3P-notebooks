"""Browser libraries served from the app rather than a public CDN.

Three properties matter here, and each would fail silently:

* **Compression.** Shiny sends static files raw. Vendoring molstar.js without gzip puts
  5.16 MB on the wire against 1.45 MB from the CDN -- 3.5x *worse*, and invisible on
  localhost where everything looks instant.
* **Portability of exports.** The apps display and download the same HTML string. A local
  ``/vendor/...`` URL is right for display and wrong for a file someone shares.
* **Falling back.** No vendor directory is the normal state for ``shiny run``, and must
  resolve to the pinned CDN URL rather than a broken link.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from common import vendor
from common.vendor import (
    ASSETS,
    URL_PREFIX,
    asset_url,
    is_vendored,
    rewrite_cdn_urls,
    static_assets,
    to_portable,
    vendor_report,
)


@pytest.fixture
def vendored(tmp_path, monkeypatch):
    """A vendor directory containing every declared asset."""
    for asset in ASSETS:
        (tmp_path / asset.filename).write_text(f"/* {asset.filename} */")
    monkeypatch.setenv("SCOP3P_VENDOR_DIR", str(tmp_path))
    return tmp_path


@pytest.fixture
def not_vendored(tmp_path, monkeypatch):
    """No vendor directory, as when running from a checkout."""
    monkeypatch.setenv("SCOP3P_VENDOR_DIR", str(tmp_path / "absent"))
    monkeypatch.setattr(vendor, "DEFAULT_VENDOR_DIR", tmp_path / "absent")
    monkeypatch.setattr(vendor, "_repo_vendor_dir", lambda: tmp_path / "also-absent")
    return None


# --------------------------------------------------------------------------------------
# The manifest
# --------------------------------------------------------------------------------------


def test_every_asset_is_pinned_to_an_exact_version() -> None:
    """``@latest`` would make behaviour depend on when the image was built."""
    for asset in ASSETS:
        assert "@latest" not in asset.url, f"{asset.key} is unpinned"
        assert any(character.isdigit() for character in asset.filename), (
            f"{asset.filename} carries no version"
        )


def test_every_asset_has_an_integrity_hash() -> None:
    """A CDN serving different bytes for a pinned version is otherwise invisible."""
    for asset in ASSETS:
        assert len(asset.sha256) == 64, f"{asset.key} has no sha256"
        int(asset.sha256, 16)


def test_asset_keys_and_filenames_are_unique() -> None:
    assert len({a.key for a in ASSETS}) == len(ASSETS)
    assert len({a.filename for a in ASSETS}) == len(ASSETS)


def test_the_filename_records_the_version_from_the_url() -> None:
    """So a bumped URL with a stale filename cannot serve the old file under a new pin."""
    import re

    for asset in ASSETS:
        match = re.search(r"@([0-9]+\.[0-9]+\.[0-9]+)", asset.url)
        if match:
            assert match.group(1) in asset.filename, (
                f"{asset.filename} does not carry the version in {asset.url}"
            )


# --------------------------------------------------------------------------------------
# Resolution
# --------------------------------------------------------------------------------------


def test_local_urls_are_used_when_the_files_are_present(vendored) -> None:
    for asset in ASSETS:
        assert asset_url(asset.key) == f"{URL_PREFIX}/{asset.filename}"
    assert set(vendor_report().values()) == {"local"}


def test_the_pinned_cdn_url_is_used_when_nothing_is_vendored(not_vendored) -> None:
    """A developer running `shiny run` has no vendor directory; that must still work."""
    for asset in ASSETS:
        assert asset_url(asset.key) == asset.url
    assert set(vendor_report().values()) == {"cdn"}


def test_an_empty_file_does_not_count_as_vendored(tmp_path, monkeypatch) -> None:
    """A truncated download would otherwise produce a 0-byte asset served as real."""
    asset = ASSETS[0]
    (tmp_path / asset.filename).write_bytes(b"")
    monkeypatch.setenv("SCOP3P_VENDOR_DIR", str(tmp_path))
    assert not is_vendored(asset)
    assert asset_url(asset.key) == asset.url


def test_static_assets_mounts_the_directory(vendored) -> None:
    mapping = static_assets()
    assert mapping == {URL_PREFIX: vendored}


def test_static_assets_is_empty_without_a_directory(not_vendored) -> None:
    """An empty mapping is what makes App() construction safe with nothing vendored."""
    assert static_assets() == {}


# --------------------------------------------------------------------------------------
# Rewriting, in both directions
# --------------------------------------------------------------------------------------


def test_cdn_urls_are_retargeted_to_local(vendored) -> None:
    """For HTML built by the topology package, which cannot import this module."""
    asset = ASSETS[0]
    html = f'<script src="{asset.url}"></script>'
    assert rewrite_cdn_urls(html) == f'<script src="{URL_PREFIX}/{asset.filename}"></script>'


def test_nothing_is_retargeted_when_nothing_is_vendored(not_vendored) -> None:
    asset = ASSETS[0]
    html = f'<script src="{asset.url}"></script>'
    assert rewrite_cdn_urls(html) == html


def test_only_exact_pinned_urls_are_retargeted(vendored) -> None:
    """A loose match could silently repoint an unrelated script."""
    other = "https://unpkg.com/ngl@9.9.9/dist/ngl.js"
    assert rewrite_cdn_urls(f'src="{other}"') == f'src="{other}"'


def test_exports_are_made_portable_again(vendored) -> None:
    """The round trip that keeps a downloaded file working on someone else's machine."""
    asset = ASSETS[0]
    displayed = rewrite_cdn_urls(f'<script src="{asset.url}"></script>')
    assert URL_PREFIX in displayed
    exported = to_portable(displayed)
    assert exported == f'<script src="{asset.url}"></script>'
    assert URL_PREFIX not in exported


def test_a_downloaded_viewer_never_points_at_a_local_path(vendored, tmp_path) -> None:
    """The concrete failure: /vendor/... resolves only inside the container.

    Exercised through ``export_html`` rather than by calling ``to_portable`` here -- a
    correct helper that the export path does not call is still a broken export, and testing
    the helper alone passes either way.
    """
    from common.viewer import NGLViewerBuilder

    html = NGLViewerBuilder.build_html(
        accession="P1",
        pdb_path=_a_pdb(),
        union_ranges=[],
        intersection_positions=[],
        modification_positions=[],
    )
    assert URL_PREFIX in html, "display should use the local copy"

    written = tmp_path / "session.html"
    NGLViewerBuilder.export_html(written, html)
    exported = written.read_text(encoding="utf-8")
    assert URL_PREFIX not in exported, "the exported file keeps a container-only path"
    assert "unpkg.com" in exported, "the exported file has no usable script source"


def test_the_download_handler_makes_its_payload_portable() -> None:
    """The other export route: a download stream rather than a written file."""
    import ast

    source = (
        Path(__file__).resolve().parents[2] / "apps" / "peptide_mapper" / "app.py"
    ).read_text(encoding="utf-8")
    tree = ast.parse(source)
    handler = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == "download_html"
    )
    body = ast.get_source_segment(source, handler) or ""
    assert "to_portable" in body, (
        "the HTML download yields the displayed payload, which points at /vendor/..."
    )


def _a_pdb() -> Path:
    import tempfile

    path = Path(tempfile.mkdtemp()) / "x.pdb"
    path.write_text("ATOM      1  CA  ALA A   1       0.0   0.0   0.0  1.00  0.00\n")
    return path


# --------------------------------------------------------------------------------------
# Compression
# --------------------------------------------------------------------------------------


def test_compression_is_enabled_on_every_app() -> None:
    """Without it, vendoring is a pessimisation rather than an optimisation."""
    import ast

    apps_dir = Path(__file__).resolve().parents[2] / "apps"
    for app_file in sorted(apps_dir.glob("*/app.py")):
        tree = ast.parse(app_file.read_text(encoding="utf-8"))
        calls = [
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        ]
        assert "enable_compression" in calls, (
            f"{app_file.parent.name} does not enable compression; Shiny serves static "
            "files raw, so molstar.js would go out at 5 MB instead of 1.45 MB"
        )
        assert "static_assets" in calls, f"{app_file.parent.name} serves no vendored assets"


def test_the_portal_compresses_and_does_not_rewrite_compressed_bodies() -> None:
    """The portal injects a navbar into the app's HTML.

    Handing it a gzipped body meant injecting into gzip bytes, which the browser rejected
    with "incorrect header check". It strips Accept-Encoding inbound and compresses its own
    output instead.
    """
    source = (
        Path(__file__).resolve().parents[2] / "apps" / "portal" / "main.py"
    ).read_text(encoding="utf-8")
    assert "_without_accept_encoding" in source
    assert "GZipMiddleware" in source
    assert source.index("GZipMiddleware") > 0


def test_compression_leaves_websockets_alone() -> None:
    """Shiny's reactivity is a websocket; compressing it would break every interaction."""
    from starlette.middleware.gzip import GZipMiddleware

    import inspect

    source = inspect.getsource(GZipMiddleware.__call__)
    assert '"http"' in source, "GZipMiddleware no longer gates on the http scope type"


# --------------------------------------------------------------------------------------
# Loading the manifest on its own
# --------------------------------------------------------------------------------------
# The image build copies apps/common/vendor.py alone to /tmp so the fetch script can read
# the manifest without the rest of the repo or the science stack. Nothing else exercises
# that, which is how a module-level Path(__file__).parents[2] shipped and failed the build
# with IndexError while merely importing the module.


def _load_manifest_from(path: Path):
    """Load the manifest the way scripts/fetch-vendor-assets.py does."""
    import importlib.util
    import sys

    spec = importlib.util.spec_from_file_location("_manifest_probe", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(spec.name, None)
    return module


def test_the_manifest_does_not_walk_parent_directories_at_import_time() -> None:
    """The rule the build failure taught, encoded directly.

    ``_REPO_VENDOR_DIR = Path(__file__).resolve().parents[2] / "vendor"`` at module level
    raised IndexError when the build copied this file to ``/tmp/vendor.py``, which has two
    parents rather than three -- and it failed while *importing*, so nothing about vendoring
    even ran.

    Checked structurally rather than by loading from a shallow path, because pytest's
    tmp_path is itself several levels deep and so cannot reproduce the shallow case: a
    depth-based test passes while the bug is present. Inside a function is fine; the module
    may be imported from anywhere, and only a call has to cope.
    """
    import ast

    tree = ast.parse(Path(vendor.__file__).read_text(encoding="utf-8"))
    offenders = []
    for node in tree.body:  # module level only
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        for child in ast.walk(node):
            if (
                isinstance(child, ast.Subscript)
                and isinstance(child.value, ast.Attribute)
                and child.value.attr == "parents"
            ):
                offenders.append(child.lineno)
    assert not offenders, (
        f"module-level .parents[...] at line(s) {offenders}: this module is copied to "
        "/tmp during the image build, where that index does not exist"
    )


def test_a_short_parent_chain_is_handled_rather_than_raising(monkeypatch) -> None:
    """The behaviour that the structural rule protects.

    Simulated by pointing the module's __file__ at a two-parent path, since a real one
    cannot be created inside a test sandbox.
    """
    monkeypatch.setattr(vendor, "__file__", "/tmp/vendor.py")
    assert vendor._repo_vendor_dir() is None


def test_the_manifest_loads_from_a_copied_location(tmp_path) -> None:
    """What the build actually does: this file alone, away from the repo."""
    source = Path(vendor.__file__).read_text(encoding="utf-8")
    target = tmp_path / "vendor.py"
    target.write_text(source, encoding="utf-8")

    module = _load_manifest_from(target)
    assert module.ASSETS, "the manifest loaded but declares no assets"
    assert module.DEFAULT_VENDOR_DIR


def test_the_manifest_still_resolves_a_directory_when_shallow(tmp_path, monkeypatch) -> None:
    """Losing the checkout-local candidate must not break resolution outright."""
    source = Path(vendor.__file__).read_text(encoding="utf-8")
    target = tmp_path / "vendor.py"
    target.write_text(source, encoding="utf-8")
    module = _load_manifest_from(target)

    # No override and no container path: it should simply report nothing vendored.
    monkeypatch.delenv("SCOP3P_VENDOR_DIR", raising=False)
    monkeypatch.setattr(module, "DEFAULT_VENDOR_DIR", tmp_path / "absent")
    assert module.vendor_dir() is None
    assert module.asset_url(module.ASSETS[0].key) == module.ASSETS[0].url


def test_the_fetch_script_needs_no_project_dependencies() -> None:
    """It runs during the build, before the app's imports are guaranteed to work.

    Importing through the ``common`` package would execute apps/common/__init__.py, which
    pulls in pandas and the rest of the science stack.
    """
    import ast

    path = Path(__file__).resolve().parents[2] / "scripts" / "fetch-vendor-assets.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))

    # Checked against the parsed imports, not the text: the script's own docstring says
    # "deliberately not ``from common.vendor import ...``", and a substring search finds
    # that and fails on the explanation for the rule it is enforcing.
    project_imports = [
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("common")
    ]
    project_imports += [
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
        if alias.name.startswith("common")
    ]
    assert not project_imports, (
        f"the fetch script imports {project_imports} through the package, which executes "
        "apps/common/__init__.py and drags in pandas"
    )
    assert "spec_from_file_location" in path.read_text(encoding="utf-8")
