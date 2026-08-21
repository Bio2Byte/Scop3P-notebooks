"""Browser libraries served from the app instead of a public CDN.

Every 3D view, network diagram and icon in the toolkit is drawn by a library fetched at
runtime. Measured from a warm connection, that was 2.07 MB across **three** hosts:

    Mol* 4.18.0        1.45 MB    2.61 s   cdn.jsdelivr.net
    NGL 2.3.1          355 KB     0.23 s   unpkg.com
    vis-network        156 KB     0.34 s   cdn.jsdelivr.net
    D3 7.9.0            93 KB     0.19 s   cdnjs.cloudflare.com
    Font Awesome 6.5.2  22 KB     0.15 s   cdnjs.cloudflare.com

Three hosts means three DNS + TCP + TLS handshakes, about 0.25 s of connection setup before
any library byte arrives -- more on a slow link, and paid on every cold load. Serving them
from the app removes the handshakes entirely and turns each fetch into a local read.

The stronger reason is not speed. Until now, a bad five minutes at unpkg or jsdelivr broke
every 3D view in the toolkit, with no more recourse than we had when UniProt was dropping
handshakes. Vendoring also makes a figure reproducible: the exact library ships in the
image, alongside the versions already recorded in ``metadata.yml``.

Two deliberate limits:

* **Falling back to the CDN is normal, not an error.** A developer running ``shiny run``
  has no vendor directory, and a half-built image should degrade rather than break. Every
  URL resolves locally when the file is present and to the pinned CDN URL when it is not.
* **Exported HTML keeps the CDN.** The RIN tab hands the user a standalone file, and pyvis
  builds it with ``cdn_resources = "remote"`` on purpose. Rewriting those to ``/vendor/...``
  would produce a file that only works inside the container -- broken the moment it is
  shared, which is the entire point of an export.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

#: URL prefix the apps mount the vendored directory at.
URL_PREFIX = "/vendor"

#: Where the image puts the downloaded files. Overridable for local experiments.
DEFAULT_VENDOR_DIR = Path("/opt/scop3p/vendor")

def _repo_vendor_dir() -> Path | None:
    """A checkout-local ``vendor/``, when this file sits at ``apps/common/vendor.py``.

    Computed on demand and guarded rather than at import time. The image build copies this
    module alone to ``/tmp/vendor.py`` so the fetch script can read the manifest without the
    rest of the repo, and there ``parents[2]`` does not exist -- which raised IndexError
    while merely *importing* the module, failing the build in a place that had nothing to do
    with vendoring.
    """
    parents = Path(__file__).resolve().parents
    if len(parents) < 3:
        return None
    return parents[2] / "vendor"


@dataclass(frozen=True, slots=True)
class VendoredAsset:
    """One pinned browser library.

    ``sha256`` is recorded so the build fails if a CDN ever serves different bytes for a
    pinned version. That has happened to other people and is invisible without a check.
    """

    key: str
    filename: str
    url: str
    sha256: str = ""

    @property
    def local_url(self) -> str:
        return f"{URL_PREFIX}/{self.filename}"


#: The manifest. One place, so a version bump cannot leave two copies disagreeing -- four
#: NGL references used to say ``@latest`` while two pinned 2.3.1, which meant two viewers
#: in one app could load different builds of the same library.
ASSETS: tuple[VendoredAsset, ...] = (
    VendoredAsset(
        key="ngl",
        filename="ngl-2.3.1.js",
        url="https://unpkg.com/ngl@2.3.1/dist/ngl.js",
        sha256="0e8fea984b0e306d948d675f30e10f5a275ab5b4ce2135191a6787ec1b29dc5d",
    ),
    VendoredAsset(
        key="molstar-js",
        filename="molstar-4.18.0.js",
        url="https://cdn.jsdelivr.net/npm/molstar@4.18.0/build/viewer/molstar.js",
        sha256="0dba8aea4c75a6816bdf900ba52e379e103715e9efe4e533cd9e8d304110a27d",
    ),
    VendoredAsset(
        key="molstar-css",
        filename="molstar-4.18.0.css",
        url="https://cdn.jsdelivr.net/npm/molstar@4.18.0/build/viewer/molstar.css",
        sha256="bea8d630b6bf8b4ad005459343a9611712deb0597eb42957e52aef5d7d594dcc",
    ),
    # Font Awesome is deliberately NOT vendored. Its CSS references eight webfont files by
    # relative path ("../webfonts/fa-solid-900.woff2"), so shipping the stylesheet alone
    # would resolve those against /vendor/ and every icon would silently fall back to a
    # blank box. Doing it properly means vendoring the fonts and rewriting the CSS, for a
    # 22 KB stylesheet on a host (cdnjs) that D3 already made us depend on -- so it is left
    # on the CDN and this comment records why, rather than someone "finishing the job".
    VendoredAsset(
        key="d3",
        filename="d3-7.9.0.min.js",
        url="https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js",
        sha256="f2094bbf6141b359722c4fe454eb6c4b0f0e42cc10cc7af921fc158fceb86539",
    ),
)

_BY_KEY = {asset.key: asset for asset in ASSETS}


def vendor_dir() -> Path | None:
    """The directory holding the downloaded libraries, or None if there is not one."""
    override = os.getenv("SCOP3P_VENDOR_DIR")
    candidates = [Path(override)] if override else []
    candidates.append(DEFAULT_VENDOR_DIR)
    repo_dir = _repo_vendor_dir()
    if repo_dir is not None:
        candidates.append(repo_dir)
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    return None


def is_vendored(asset: VendoredAsset) -> bool:
    directory = vendor_dir()
    if directory is None:
        return False
    path = directory / asset.filename
    return path.is_file() and path.stat().st_size > 0


def asset_url(key: str) -> str:
    """Where the browser should fetch this library from.

    The local copy when it exists, the pinned CDN URL otherwise. Callers never need to know
    which, and the pinned URL means a fallback still loads the version that was tested.
    """
    asset = _BY_KEY[key]
    return asset.local_url if is_vendored(asset) else asset.url


def static_assets() -> dict[str, Path]:
    """The ``static_assets`` mapping for ``shiny.App``, empty when nothing is vendored.

    Every app mounts the same prefix, so ``/vendor/...`` resolves whichever app the portal
    happens to be serving -- no portal-level route is needed.
    """
    directory = vendor_dir()
    if directory is None:
        return {}
    return {URL_PREFIX: directory}


def rewrite_cdn_urls(html: str) -> str:
    """Point any pinned CDN URL in a fragment of HTML at the local copy.

    For HTML built somewhere that cannot import this module. The topology viewer's JS lives
    in ``notebooks/topology_viewer/topology``, which is shared with the Voila notebook and
    its test suite and must not depend on ``apps/`` -- so its Mol* and NGL URLs are
    rewritten here, at the boundary where the Shiny app embeds them, rather than by changing
    a package other people also run.

    Only exact pinned URLs are replaced, so this cannot silently retarget something else.
    """
    for asset in ASSETS:
        if is_vendored(asset):
            html = html.replace(asset.url, asset.local_url)
    return html


def to_portable(html: str) -> str:
    """Point local ``/vendor/...`` URLs back at the pinned CDN, for HTML that leaves here.

    The apps display and *export* the same HTML string -- Peptide Mapper stores one payload
    and both renders it and offers it as a download. A local URL is right for the first and
    wrong for the second: ``/vendor/ngl-2.3.1.js`` resolves only inside the container, so an
    exported file would work for the person who made it and be broken for everyone they
    send it to, which defeats the purpose of an export.

    Applied on the way out, so display keeps the fast local copy and the file people share
    keeps working.
    """
    for asset in ASSETS:
        html = html.replace(asset.local_url, asset.url)
    return html


def vendor_report() -> dict[str, str]:
    """Per-asset resolution, for a startup log line or a test."""
    return {
        asset.key: "local" if is_vendored(asset) else "cdn" for asset in ASSETS
    }


# ---------------------------------------------------------------------------
# Compression
# ---------------------------------------------------------------------------
# Serving these files locally is only faster if they are compressed. Shiny's static handler
# sends them raw: measured, molstar.js went out as 5,159,386 bytes with no Content-Encoding
# at all, against 1.45 MB gzipped from the CDN. Vendoring without this would put 3.5x *more*
# bytes on the wire for any remote user -- the opposite of the intent, and invisible on
# localhost where it all looks instant.
#
# Compressing also covers the page itself, which is worth as much again: the HTML is ~380 KB
# because the footer logos are inlined as data URIs, and base64 compresses well.

#: Below this, compression costs more than it saves.
COMPRESSION_MINIMUM_BYTES = 1024


def enable_compression(app) -> None:  # noqa: ANN001 - a shiny.App, kept loose for testing
    """Add gzip to a Shiny app's response path.

    Applied to the Starlette app underneath rather than by wrapping the Shiny app, because
    ``shiny run`` expects to be handed a ``shiny.App`` and a wrapped ASGI callable would not
    give it one. The middleware stack is built lazily on first request, so adding to it after
    construction is fine.

    GZipMiddleware only touches ``scope["type"] == "http"``, so the websocket Shiny uses for
    its reactive updates passes through untouched -- worth stating because breaking that
    would disable every interaction rather than merely slowing it.
    """
    from starlette.middleware.gzip import GZipMiddleware

    inner = getattr(app, "starlette_app", None)
    if inner is None or not hasattr(inner, "add_middleware"):
        return
    inner.add_middleware(GZipMiddleware, minimum_size=COMPRESSION_MINIMUM_BYTES)
