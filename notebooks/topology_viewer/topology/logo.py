"""The app's logo.

Built from the same grammar the viewer draws with -- strand arrows, a helix
capsule, connector loops, N and C termini, one annotated site -- rather than a
generic molecule glyph. It reads as a miniature of the actual output, and it
reuses the palette, so the mark and the tool stay in step if either changes.

The topology is real: two antiparallel strands joined by a loop below, then a
loop above into a helix, N to C. Nothing here is arranged only for looks.
"""

from __future__ import annotations

# Same values the renderer uses.
HELIX = "#d9606b"
STRAND_A = "#5fb9e0"
STRAND_B = "#e8912d"
CONNECTOR = "#6b7887"
SITE = "#F0C808"
INK = "#16202b"
MUTED = "#5f6b7a"
PANEL = "#f4f7fa"
LINE = "#e4ebf2"

_SHAFT = 7.0
_HEAD_W = 15.0
_HEAD_L = 16.0


def _arrow_down(x: float, top: float, bottom: float) -> str:
    points = [
        (x - _SHAFT, top),
        (x - _SHAFT, bottom - _HEAD_L),
        (x - _HEAD_W, bottom - _HEAD_L),
        (x, bottom),
        (x + _HEAD_W, bottom - _HEAD_L),
        (x + _SHAFT, bottom - _HEAD_L),
        (x + _SHAFT, top),
    ]
    return " ".join(f"{px:g},{py:g}" for px, py in points)


def _arrow_up(x: float, top: float, bottom: float) -> str:
    points = [
        (x + _SHAFT, bottom),
        (x + _SHAFT, top + _HEAD_L),
        (x + _HEAD_W, top + _HEAD_L),
        (x, top),
        (x - _HEAD_W, top + _HEAD_L),
        (x - _SHAFT, top + _HEAD_L),
        (x - _SHAFT, bottom),
    ]
    return " ".join(f"{px:g},{py:g}" for px, py in points)


def _glyph(background: bool = True, rounded: bool = True, offset_x: float = 0.0) -> str:
    """The drawing itself, without an <svg> wrapper.

    Kept separate so the icon and the lockup share one definition instead of one
    being carved out of the other's markup.
    """
    plate = ""
    if background:
        radius = 26 if rounded else 0
        plate = (
            f'<rect x="{2 + offset_x:g}" y="2" width="124" height="124" '
            f'rx="{radius}" ry="{radius}" '
            f'fill="{PANEL}" stroke="{LINE}" stroke-width="2"/>'
        )

    return f"""{plate}
  <g transform="translate({offset_x:g} 0)" stroke-linejoin="round">
    <path d="M 26 98 L 26 110 L 64 110 L 64 98" fill="none"
          stroke="{CONNECTOR}" stroke-width="3"/>
    <path d="M 64 30 L 64 18 L 102 18 L 102 30" fill="none"
          stroke="{CONNECTOR}" stroke-width="3"/>
    <polygon points="{_arrow_down(26, 30, 98)}"
             fill="{STRAND_A}" stroke="{INK}" stroke-width="2.4"/>
    <polygon points="{_arrow_up(64, 30, 98)}"
             fill="{STRAND_B}" stroke="{INK}" stroke-width="2.4"/>
    <rect x="93" y="30" width="18" height="68" rx="9" ry="9"
          fill="{HELIX}" stroke="{INK}" stroke-width="2.4"/>
    <circle cx="102" cy="58" r="6" fill="{SITE}" stroke="#ffffff" stroke-width="2"/>
    <text x="26" y="22" font-family="Inter, 'Segoe UI', system-ui, sans-serif"
          font-size="13" font-weight="700" fill="{MUTED}" text-anchor="middle">N</text>
    <text x="102" y="116" font-family="Inter, 'Segoe UI', system-ui, sans-serif"
          font-size="13" font-weight="700" fill="{MUTED}" text-anchor="middle">C</text>
  </g>"""


def mark(size: int = 128, background: bool = True, rounded: bool = True) -> str:
    """The square icon on its own."""
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 128 128" \
width="{size}" height="{size}" role="img" aria-label="Protein topology viewer">
  <title>Protein topology viewer</title>
  {_glyph(background, rounded)}
</svg>"""


def lockup(height: int = 56, subtitle: str = "topology viewer") -> str:
    """Icon plus wordmark, for a page header."""
    total_width = 470
    width = int(height * total_width / 128.0)

    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_width} 128" \
height="{height}" width="{width}" role="img" aria-label="Protein topology viewer">
  <title>Protein topology viewer</title>
  {_glyph(True, True)}
  <text x="154" y="60" font-family="Inter, 'Segoe UI', system-ui, sans-serif"
        font-size="34" font-weight="650" fill="{INK}">Protein</text>
  <text x="154" y="100" font-family="Inter, 'Segoe UI', system-ui, sans-serif"
        font-size="34" font-weight="650" fill="{HELIX}">{subtitle}</text>
</svg>"""


def favicon() -> str:
    """A stripped-back mark for small sizes.

    The full logo carries termini labels, a site marker and outlines that turn
    into noise below about 32px, so the small version drops them rather than
    scaling them down into mush.
    """
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 128 128" \
width="32" height="32" role="img" aria-label="Topology viewer">
  <rect width="128" height="128" rx="26" ry="26" fill="{PANEL}"/>
  <path d="M 30 96 L 30 108 L 64 108 L 64 96" fill="none"
        stroke="{CONNECTOR}" stroke-width="5"/>
  <path d="M 64 32 L 64 20 L 98 20 L 98 32" fill="none"
        stroke="{CONNECTOR}" stroke-width="5"/>
  <polygon points="{_arrow_down(30, 32, 96)}" fill="{STRAND_A}"/>
  <polygon points="{_arrow_up(64, 32, 96)}" fill="{STRAND_B}"/>
  <rect x="89" y="32" width="18" height="64" rx="9" ry="9" fill="{HELIX}"/>
</svg>"""
