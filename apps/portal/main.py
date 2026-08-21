from __future__ import annotations

import base64
from http.cookies import SimpleCookie
from pathlib import Path
from typing import Any
from urllib.parse import parse_qsl, urlencode

from starlette.datastructures import Headers, MutableHeaders

from starlette.middleware.gzip import GZipMiddleware

from common import vendor
from common.logging_utils import get_logger
from common.ui_shell import CITATION_DOI_URL
from help.app import app as help_app
from mutation_effect.app import app as mutation_effect_app
from peptide_mapper.app import app as peptide_mapper_app
from rinalign.app import app as rinalign_app
from structure_viz.app import app as structure_viz_app
from topology_viewer.app import app as topology_viewer_app


LOGGER = get_logger("scop3p.portal")
# Dict order is navbar order. The key is the ?app= value and the cookie value.
APP_OPTIONS = {
    "peptide-mapper": ("Peptide Mapper", "fa-solid fa-map-pin", peptide_mapper_app),
    "structure-viz": ("Structure Visualisation", "fa-solid fa-cube", structure_viz_app),
    "topology-viewer": ("Topology Viewer", "fa-solid fa-diagram-project", topology_viewer_app),
    "mutation-effect": ("Mutation Effect", "fa-solid fa-bolt", mutation_effect_app),
    "rinalign": ("RIN Alignment", "fa-solid fa-circle-nodes", rinalign_app),
    # Last on purpose: help sits at the end of the navbar, where help belongs.
    "help": ("Help?", "fa-solid fa-circle-question", help_app),
}
DEFAULT_APP_KEY = "peptide-mapper"
COOKIE_NAME = "scop3p_app"
_LOGO_PATH = Path(__file__).resolve().parents[1] / "assets" / "images" / "scop3p-nobg.png"


def _normalize_app_key(value: str | None) -> str:
    return value if value in APP_OPTIONS else DEFAULT_APP_KEY


def _get_selected_app_key(scope: dict[str, Any]) -> str:
    query_items = dict(parse_qsl(scope.get("query_string", b"").decode("utf-8"), keep_blank_values=True))
    if "app" in query_items:
        selected = _normalize_app_key(query_items["app"])
        LOGGER.info(
            "portal navbar clicked requested_app=%s selected_app=%s",
            query_items["app"],
            selected,
            extra={"event": "navbar_click"},
        )
        return selected

    headers = Headers(scope=scope)
    cookie_header = headers.get("cookie", "")
    cookie = SimpleCookie()
    cookie.load(cookie_header)
    if COOKIE_NAME in cookie:
        return _normalize_app_key(cookie[COOKIE_NAME].value)
    return DEFAULT_APP_KEY


def _strip_selector_query(scope: dict[str, Any]) -> dict[str, Any]:
    query_items = parse_qsl(scope.get("query_string", b"").decode("utf-8"), keep_blank_values=True)
    filtered = [(key, value) for key, value in query_items if key != "app"]
    if len(filtered) == len(query_items):
        return scope

    new_scope = dict(scope)
    new_scope["query_string"] = urlencode(filtered, doseq=True).encode("utf-8")
    return new_scope


def _logo_data_uri() -> str:
    payload = base64.b64encode(_LOGO_PATH.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{payload}"


def _selector_navbar(selected_key: str) -> str:
    links = []
    for key, (label, icon_class, _) in APP_OPTIONS.items():
        active_class = " toolkit-link-active" if key == selected_key else ""
        links.append(
            f'<a class="toolkit-link{active_class}" href="/?app={key}">'
            f'<i class="{icon_class}" aria-hidden="true"></i>'
            f'<span>{label}</span>'
            f"</a>"
        )
    logo_src = _logo_data_uri()

    return f"""
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.2/css/all.min.css">
<style>
  body {{
    padding-top: 78px !important;
  }}
  .toolkit-navbar {{
    position: fixed;
    inset: 0 0 auto 0;
    z-index: 9999;
    background: rgba(16, 38, 60, 0.96);
    color: #fff;
    border-bottom: 1px solid rgba(255,255,255,0.08);
    backdrop-filter: blur(10px);
    box-shadow: 0 10px 24px rgba(16, 38, 60, 0.16);
  }}
  .toolkit-navbar-inner {{
    /* Same width as the content below, from the property ui_shell defines. The literal
       fallback matters: an app served without the portal never loads this CSS, but the
       portal's own selector page does, and it has no ui_shell styles to read from. */
    max-width: var(--scop3p-max-width, 98vw);
    margin: 0 auto;
    padding: 14px 18px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 18px;
  }}
  .toolkit-brand {{
    display: flex;
    align-items: center;
    gap: 12px;
  }}
  .toolkit-brand-copy {{
    display: flex;
    flex-direction: column;
    gap: 2px;
  }}
  .toolkit-brand-logo {{
    display: block;
    height: 42px;
    width: auto;
    object-fit: contain;
  }}
  .toolkit-brand strong {{
    font-size: 1.1rem;
    letter-spacing: 0.02em;
  }}
  .toolkit-brand span {{
    color: rgba(255,255,255,0.72);
    font-size: 0.88rem;
  }}
  .toolkit-preprint {{
    color: rgba(255,255,255,0.92);
    text-decoration: underline;
    text-decoration-color: rgba(255,255,255,0.35);
    text-underline-offset: 2px;
  }}
  .toolkit-preprint:hover {{
    text-decoration-color: rgba(255,255,255,0.9);
  }}
  .toolkit-links {{
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
  }}
  .toolkit-link {{
    display: inline-flex;
    align-items: center;
    gap: 10px;
    text-decoration: none;
    color: rgba(255,255,255,0.84);
    padding: 9px 14px;
    border-radius: 999px;
    border: 1px solid rgba(255,255,255,0.12);
    background: rgba(255,255,255,0.06);
    font-weight: 600;
  }}
  .toolkit-link-active {{
    color: #11263c;
    background: #ffffff;
    border-color: #ffffff;
  }}
  .toolkit-link i {{
    font-size: 0.95rem;
  }}
  @media (max-width: 1000px) {{
    .toolkit-navbar-inner {{
      flex-direction: column;
      align-items: flex-start;
    }}
    body {{
      padding-top: 124px !important;
    }}
  }}
</style>
<div class="toolkit-navbar">
  <div class="toolkit-navbar-inner">
    <div class="toolkit-brand">
      <img class="toolkit-brand-logo" src="{logo_src}" alt="Scop3P logo" />
      <div class="toolkit-brand-copy">
        <strong>Scop3P-Toolkit</strong>
        <span>Tools for exploring and extending Scop3P (<a class="toolkit-preprint" href="{CITATION_DOI_URL}" target="_blank" rel="noopener noreferrer">read pre-print</a>)</span>
      </div>
    </div>
    <nav class="toolkit-links" aria-label="Tools for exploring and extending Scop3P">
      {''.join(links)}
    </nav>
  </div>
</div>
"""


def _inject_navbar(html_text: str, selected_key: str) -> str:
    navbar = _selector_navbar(selected_key)
    if "<body>" in html_text:
        return html_text.replace("<body>", f"<body>{navbar}", 1)
    if "<body " in html_text:
        marker = html_text.find(">")
        if marker != -1:
            return html_text[: marker + 1] + navbar + html_text[marker + 1 :]
    return navbar + html_text


def _set_selection_cookie(headers: MutableHeaders, selected_key: str) -> None:
    headers.append(
        "set-cookie",
        f"{COOKIE_NAME}={selected_key}; Path=/; HttpOnly; SameSite=Lax",
    )


class SingleRootPortal:
    def __init__(self) -> None:
        self.apps = {key: app for key, (_, _, app) in APP_OPTIONS.items()}
        LOGGER.info("portal initialized apps=%s", ",".join(self.apps), extra={"event": "portal_startup"})

    async def __call__(self, scope, receive, send) -> None:  # noqa: ANN001
        scope_type = scope["type"]
        if scope_type == "lifespan":
            await self._handle_lifespan(receive, send)
            return

        selected_key = _get_selected_app_key(scope)
        selected_app = self.apps[selected_key]
        delegated_scope = _strip_selector_query(scope)
        LOGGER.info(
            "portal dispatch scope=%s selected_app=%s path=%s",
            scope_type,
            selected_key,
            scope.get("path", "-"),
            extra={"event": "portal_dispatch"},
        )

        if scope_type == "http":
            await self._dispatch_http(selected_app, delegated_scope, receive, send, selected_key)
            return

        await selected_app(delegated_scope, receive, send)

    async def _dispatch_http(self, app, scope, receive, send, selected_key: str) -> None:  # noqa: ANN001
        path = scope.get("path", "")
        method = scope.get("method", "GET").upper()
        should_inject = path == "/" and method == "GET"

        # This portal rewrites the app's HTML to inject the navbar, so it must never be
        # handed a compressed body -- injecting into gzip bytes produces something the
        # browser cannot decode ("incorrect header check"). Dropping Accept-Encoding on the
        # way in keeps the app's response plain; the portal is itself wrapped in
        # GZipMiddleware, so what reaches the client is still compressed.
        scope = _without_accept_encoding(scope)

        if not should_inject:
            async def passthrough(message):  # noqa: ANN001
                if message["type"] == "http.response.start":
                    headers = MutableHeaders(raw=list(message["headers"]))
                    _set_selection_cookie(headers, selected_key)
                    message = {**message, "headers": headers.raw}
                await send(message)

            await app(scope, receive, passthrough)
            return

        response_start: dict[str, Any] | None = None
        body_parts: list[bytes] = []

        async def capture(message):  # noqa: ANN001
            nonlocal response_start
            if message["type"] == "http.response.start":
                response_start = message
                return

            if message["type"] != "http.response.body":
                await send(message)
                return

            body_parts.append(message.get("body", b""))
            if message.get("more_body", False):
                return

            if response_start is None:
                raise RuntimeError("Missing http.response.start for portal response.")

            headers = MutableHeaders(raw=list(response_start["headers"]))
            _set_selection_cookie(headers, selected_key)
            body = b"".join(body_parts)
            if "text/html" in headers.get("content-type", "").lower():
                injected = _inject_navbar(body.decode("utf-8", errors="ignore"), selected_key)
                body = injected.encode("utf-8")
                headers["content-length"] = str(len(body))

            await send({**response_start, "headers": headers.raw})
            await send({"type": "http.response.body", "body": body})

        await app(scope, receive, capture)

    async def _handle_lifespan(self, receive, send) -> None:  # noqa: ANN001
        while True:
            message = await receive()
            if message["type"] == "lifespan.startup":
                LOGGER.info("portal lifespan startup", extra={"event": "portal_lifespan"})
                await send({"type": "lifespan.startup.complete"})
            elif message["type"] == "lifespan.shutdown":
                LOGGER.info("portal lifespan shutdown", extra={"event": "portal_lifespan"})
                await send({"type": "lifespan.shutdown.complete"})
                return


def _without_accept_encoding(scope: dict) -> dict:
    """A copy of the scope with the Accept-Encoding request header removed."""
    headers = [
        (name, value)
        for name, value in scope.get("headers", [])
        if name.lower() != b"accept-encoding"
    ]
    return {**scope, "headers": headers}


# Compression is applied here, outside the navbar injection, so it acts on the final body.
# It covers both the HTML (~380 KB, and base64-heavy, so it compresses well) and the
# vendored browser libraries proxied through from the selected app.
app = GZipMiddleware(
    SingleRootPortal(), minimum_size=vendor.COMPRESSION_MINIMUM_BYTES
)
