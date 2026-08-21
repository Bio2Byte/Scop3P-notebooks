"""Images inlined into the page, and the size budget that implies.

Every image in the shell is a base64 data URI, so its bytes are paid on *every* page load --
there is no separate cacheable request to amortise them. The partner logos are print assets
(ELIXIR Belgium is 2656x1752, 498 KB as base64) displayed at 60-120px tall, and together they
were putting 1.07 MB of HTML in front of every user on every load of every protocol.

The budget assertions below are the point of this file. They are deliberately expressed as
limits rather than exact sizes, so replacing a logo with another print-resolution file fails
here rather than quietly costing a megabyte again.
"""

from __future__ import annotations

import base64
import io

import pytest

from common.ui_shell import (
    _FOOTER_LOGOS,
    _IMAGE_DIR,
    _RETINA_SCALE,
    _css_pixels,
    _image_data_uri,
    _scaled_image_data_uri,
    scop3p_footer,
)

Image = pytest.importorskip("PIL.Image")


#: The whole footer, logos included, must stay under this. It was 1094 KB.
FOOTER_BUDGET_KB = 200

#: No single inlined logo may exceed this.
LOGO_BUDGET_KB = 40


def _kb(text: str | None) -> float:
    return len(text or "") / 1024


# --------------------------------------------------------------------------------------
# The budget
# --------------------------------------------------------------------------------------


def test_the_footer_stays_within_its_page_weight_budget() -> None:
    size = _kb(str(scop3p_footer()))
    assert size < FOOTER_BUDGET_KB, (
        f"the footer is {size:.0f} KB of HTML on every page load; it was reduced from "
        "1094 KB by resizing the logos to the height they are displayed at"
    )


@pytest.mark.parametrize(
    "label, filename, height", _FOOTER_LOGOS, ids=[row[0] for row in _FOOTER_LOGOS]
)
def test_no_single_logo_blows_the_budget(label, filename, height) -> None:
    size = _kb(_scaled_image_data_uri(filename, height))
    assert size < LOGO_BUDGET_KB, f"{label} inlines at {size:.0f} KB"


def test_resizing_actually_saves_something_substantial() -> None:
    """Guards against the resize silently becoming a no-op."""
    before = sum(_kb(_image_data_uri(f)) for _l, f, _h in _FOOTER_LOGOS)
    after = sum(_kb(_scaled_image_data_uri(f, h)) for _l, f, h in _FOOTER_LOGOS)
    assert after < before * 0.3, (
        f"resizing saved only {100 * (1 - after / before):.0f}% "
        f"({before:.0f} KB -> {after:.0f} KB)"
    )


# --------------------------------------------------------------------------------------
# Correctness of the resize
# --------------------------------------------------------------------------------------


def _decode(uri: str):
    return Image.open(io.BytesIO(base64.b64decode(uri.split(",", 1)[1])))


@pytest.mark.parametrize(
    "label, filename, height", _FOOTER_LOGOS, ids=[row[0] for row in _FOOTER_LOGOS]
)
def test_a_logo_is_rendered_at_the_height_it_is_displayed(label, filename, height) -> None:
    """Smaller than displayed would look soft; larger is wasted bytes."""
    target = _css_pixels(height) * _RETINA_SCALE
    with _decode(_scaled_image_data_uri(filename, height)) as image:
        assert image.height <= target, f"{label} is taller than it is drawn"
        with Image.open(_IMAGE_DIR / filename) as source:
            # Only shrunk, never enlarged: upscaling adds bytes and no detail.
            assert image.height == min(target, source.height)


@pytest.mark.parametrize(
    "label, filename, height", _FOOTER_LOGOS, ids=[row[0] for row in _FOOTER_LOGOS]
)
def test_the_aspect_ratio_survives(label, filename, height) -> None:
    """A squashed partner logo is worse than a heavy one."""
    with Image.open(_IMAGE_DIR / filename) as source:
        source_ratio = source.width / source.height
    with _decode(_scaled_image_data_uri(filename, height)) as image:
        assert abs(image.width / image.height - source_ratio) < 0.02, f"{label} distorted"


def test_transparency_is_preserved() -> None:
    """These are logos on a dark footer; losing alpha would show a white box."""
    with _decode(_scaled_image_data_uri("bio2byte.png", "60px")) as image:
        assert image.mode == "RGBA"


def test_every_logo_still_renders() -> None:
    rendered = str(scop3p_footer())
    for label, _filename, _height in _FOOTER_LOGOS:
        assert f'alt="{label}"' in rendered, f"{label} is missing from the footer"
    assert "scop3p-logo-fallback" not in rendered, "a logo fell back to a text placeholder"


# --------------------------------------------------------------------------------------
# Caching and degradation
# --------------------------------------------------------------------------------------


def test_the_resize_happens_once_per_process() -> None:
    """Resizing per page load would trade bandwidth for CPU on every request."""
    from common import ui_shell

    ui_shell._scaled_cache.clear()
    first = _scaled_image_data_uri("bio2byte.png", "60px")
    assert ui_shell._scaled_cache, "nothing was cached"
    calls = {"n": 0}
    real = ui_shell._resized_png

    def counting(*args, **kwargs):
        calls["n"] += 1
        return real(*args, **kwargs)

    ui_shell._resized_png = counting
    try:
        second = _scaled_image_data_uri("bio2byte.png", "60px")
    finally:
        ui_shell._resized_png = real
    assert second == first
    assert calls["n"] == 0, "the second call re-encoded the image"


def test_a_missing_logo_degrades_to_a_label(monkeypatch, tmp_path) -> None:
    from common import ui_shell

    monkeypatch.setattr(ui_shell, "_IMAGE_DIR", tmp_path)
    monkeypatch.setattr(ui_shell, "_scaled_cache", {})
    rendered = str(ui_shell.scop3p_footer())
    assert "scop3p-logo-fallback" in rendered
    assert "CompOmics" in rendered


def test_without_pillow_the_page_still_renders(monkeypatch) -> None:
    """Pillow is transitive (via bokeh), so it may simply not be installed."""
    import builtins

    from common import ui_shell

    real_import = builtins.__import__

    def no_pil(name, *args, **kwargs):
        if name.startswith("PIL"):
            raise ImportError("no Pillow")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(ui_shell, "_scaled_cache", {})
    monkeypatch.setattr(builtins, "__import__", no_pil)
    try:
        uri = _scaled_image_data_uri("bio2byte.png", "60px")
    finally:
        monkeypatch.undo()
        ui_shell._scaled_cache.clear()
    assert uri is not None and uri.startswith("data:image/png;base64,")


@pytest.mark.parametrize(
    "value, expected",
    [
        ("60px", 60),
        ("120px", 120),
        ("60", 60),
        (" 60 px ", 60),   # tolerant of stray whitespace, which is the useful behaviour
        ("auto", None),
        ("", None),
        ("3rem", None),    # a unit this cannot convert must not guess a pixel count
    ],
)
def test_css_pixel_parsing(value, expected) -> None:
    assert _css_pixels(value) == expected


def test_an_unparseable_height_falls_back_to_the_full_image() -> None:
    """Better a heavy logo than a wrongly-sized one, and never a missing one."""
    assert _scaled_image_data_uri("bio2byte.png", "auto") == _image_data_uri("bio2byte.png")
