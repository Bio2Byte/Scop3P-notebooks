from __future__ import annotations

import base64
from http.cookies import SimpleCookie
from pathlib import Path
from typing import Any
from urllib.parse import parse_qsl, urlencode

from starlette.datastructures import Headers, MutableHeaders

from mutation_effect.app import app as mutation_effect_app
from peptide_mapper.app import app as peptide_mapper_app
from structure_viz.app import app as structure_viz_app


APP_OPTIONS = {
    "peptide-mapper": ("Peptide Mapper", "fa-solid fa-map-pin", peptide_mapper_app),
    "structure-viz": ("Structure Visualisation", "fa-solid fa-cube", structure_viz_app),
    "mutation-effect": ("Mutation Effect", "fa-solid fa-bolt", mutation_effect_app),
}
DEFAULT_APP_KEY = "peptide-mapper"
COOKIE_NAME = "scop3p_app"
_LOGO_PATH = Path(__file__).resolve().parents[1] / "assets" / "images" / "scop3p-nobg.png"


def _normalize_app_key(value: str | None) -> str:
    return value if value in APP_OPTIONS else DEFAULT_APP_KEY


def _get_selected_app_key(scope: dict[str, Any]) -> str:
    query_items = dict(parse_qsl(scope.get("query_string", b"").decode("utf-8"), keep_blank_values=True))
    if "app" in query_items:
        return _normalize_app_key(query_items["app"])

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
        active_class = " suite-link-active" if key == selected_key else ""
        links.append(
            f'<a class="suite-link{active_class}" href="/?app={key}">'
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
  .suite-navbar {{
    position: fixed;
    inset: 0 0 auto 0;
    z-index: 9999;
    background: rgba(16, 38, 60, 0.96);
    color: #fff;
    border-bottom: 1px solid rgba(255,255,255,0.08);
    backdrop-filter: blur(10px);
    box-shadow: 0 10px 24px rgba(16, 38, 60, 0.16);
  }}
  .suite-navbar-inner {{
    max-width: 1800px;
    margin: 0 auto;
    padding: 14px 18px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 18px;
  }}
  .suite-brand {{
    display: flex;
    align-items: center;
    gap: 12px;
  }}
  .suite-brand-copy {{
    display: flex;
    flex-direction: column;
    gap: 2px;
  }}
  .suite-brand-logo {{
    display: block;
    height: 42px;
    width: auto;
    object-fit: contain;
  }}
  .suite-brand strong {{
    font-size: 1.1rem;
    letter-spacing: 0.02em;
  }}
  .suite-brand span {{
    color: rgba(255,255,255,0.72);
    font-size: 0.88rem;
  }}
  .suite-links {{
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
  }}
  .suite-link {{
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
  .suite-link-active {{
    color: #11263c;
    background: #ffffff;
    border-color: #ffffff;
  }}
  .suite-link i {{
    font-size: 0.95rem;
  }}
  @media (max-width: 1000px) {{
    .suite-navbar-inner {{
      flex-direction: column;
      align-items: flex-start;
    }}
    body {{
      padding-top: 124px !important;
    }}
  }}
</style>
<div class="suite-navbar">
  <div class="suite-navbar-inner">
    <div class="suite-brand">
      <img class="suite-brand-logo" src="{logo_src}" alt="Scop3P logo" />
      <div class="suite-brand-copy">
        <strong>Scop3P-Toolkit</strong>
        <span>Tools for exploring and extending Scop3P</span>
      </div>
    </div>
    <nav class="suite-links" aria-label="Tools for exploring and extending Scop3P">
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

    async def __call__(self, scope, receive, send) -> None:  # noqa: ANN001
        scope_type = scope["type"]
        if scope_type == "lifespan":
            await self._handle_lifespan(receive, send)
            return

        selected_key = _get_selected_app_key(scope)
        selected_app = self.apps[selected_key]
        delegated_scope = _strip_selector_query(scope)

        if scope_type == "http":
            await self._dispatch_http(selected_app, delegated_scope, receive, send, selected_key)
            return

        await selected_app(delegated_scope, receive, send)

    async def _dispatch_http(self, app, scope, receive, send, selected_key: str) -> None:  # noqa: ANN001
        path = scope.get("path", "")
        method = scope.get("method", "GET").upper()
        should_inject = path == "/" and method == "GET"

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
                await send({"type": "lifespan.startup.complete"})
            elif message["type"] == "lifespan.shutdown":
                await send({"type": "lifespan.shutdown.complete"})
                return


app = SingleRootPortal()
