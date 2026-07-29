"""Render the layout to standalone SVG so it can be inspected as an image."""
import sys
from pathlib import Path
from topology.io import load_structure
from topology.elements import build_elements, strand_contacts, annotate_geometry
from topology.layout import build_layout

HELIX = "#d9606b"
SHEET = ["#e8912d", "#5fb9e0", "#e0563f", "#6fbb6f", "#9186cf", "#cf9a63", "#3fa39b", "#c76fae"]
LONE = "#8fa2b8"

def rounded(points, radius=12):
    if len(points) < 4: return ""
    parts = [f"M {points[0]} {points[1]}"]
    for i in range(2, len(points)-2, 2):
        px,py = points[i-2],points[i-1]; cx,cy = points[i],points[i+1]; nx,ny = points[i+2],points[i+3]
        li = ((cx-px)**2+(cy-py)**2)**.5; lo = ((nx-cx)**2+(ny-cy)**2)**.5
        r = min(radius, li/2 if li else 0, lo/2 if lo else 0)
        if r <= 0.5: parts.append(f"L {cx} {cy}"); continue
        t1, t2 = r/li, r/lo
        parts.append(f"L {cx-(cx-px)*t1} {cy-(cy-py)*t1}")
        parts.append(f"Q {cx} {cy} {cx+(nx-cx)*t2} {cy+(ny-cy)*t2}")
    parts.append(f"L {points[-2]} {points[-1]}")
    return " ".join(parts)

def svg_for(L, title):
    x0,y0,x1,y1 = L["extents"]; pad=52
    w,h = x1-x0+2*pad, y1-y0+2*pad
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0-pad} {y0-pad} {w} {h}" width="{w:.0f}" height="{h:.0f}">',
           f'<rect x="{x0-pad}" y="{y0-pad}" width="{w}" height="{h}" fill="#ffffff"/>',
           f'<text x="{x0-pad+14}" y="{y0-pad+26}" font-family="Inter,sans-serif" font-size="15" font-weight="600" fill="#16202b">{title}</text>']
    for sg in L.get("segments", []):
        out.append(f'<rect x="{sg["x"]}" y="{sg["y"]}" width="{sg["width"]}" height="{sg["height"]}" rx="10" fill="#f4f7fa" stroke="#e4ebf2"/>')
        out.append(f'<text x="{sg["x"]+12}" y="{sg["y"]+15}" font-family="Inter,sans-serif" font-size="10" font-weight="600" fill="#93a1b3">{sg["label"]}</text>')
    for c in L["connectors"]:
        out.append(f'<path d="{rounded(c["path"], 4)}" fill="none" stroke="#6b7887" stroke-width="1.4" stroke-linejoin="miter"/>')
    for t in L["termini"]:
        out.append(f'<path d="M {t["anchor"][0]} {t["anchor"][1]} L {t["x"]} {t["y"]}" stroke="#93a1b3" stroke-width="1.8" fill="none"/>')
        dy = -12 if t["type"]=="N" else 12
        out.append(f'<text x="{t["x"]}" y="{t["y"]+dy}" font-family="Inter,sans-serif" font-size="12" font-weight="700" fill="#5f6b7a" text-anchor="middle" dominant-baseline="middle">{t["type"]}</text>')
    for e in L["elements"]:
        meta = META.get(e["id"], {})
        if e["kind"] == "helix":
            f = HELIX
        elif meta.get("sheet_index"):
            f = SHEET[(meta["sheet_index"] - 1) % len(SHEET)]
        else:
            f = LONE
        if e["kind"]=="strand":
            pts = " ".join(f'{e["path"][i]},{e["path"][i+1]}' for i in range(0,len(e["path"]),2))
            out.append(f'<polygon points="{pts}" fill="{f}" stroke="#16202b" stroke-width="1.4"/>')
        else:
            out.append(f'<rect x="{e["x"]-16}" y="{e["y"]}" width="32" height="{e["h"]}" rx="15" fill="{f}" stroke="#16202b" stroke-width="1.4"/>')
        if e["kind"] == "strand" and meta.get("sheet"):
            by = e["y"] + e["h"]/2
            out.append(f'<rect x="{e["x"]-13}" y="{by-9}" width="26" height="18" rx="5" fill="rgba(255,255,255,0.92)" stroke="rgba(0,0,0,0.18)" stroke-width="0.8"/>')
            out.append(f'<text x="{e["x"]}" y="{by}" font-family="Inter,sans-serif" font-size="10" font-weight="700" fill="#1d2733" text-anchor="middle" dominant-baseline="middle">{meta["sheet"]}</text>')
        out.append(f'<text x="{e["x"]}" y="{e["y"]-11}" font-family="Inter,sans-serif" font-size="11.5" font-weight="700" fill="#16202b" text-anchor="middle" dominant-baseline="middle">{e["id"]}</text>')
        tn = e["start"] if e["direction"]>0 else e["stop"]
        bn = e["stop"] if e["direction"]>0 else e["start"]
        out.append(f'<text x="{e["x"]-20}" y="{e["y"]+6}" font-family="Inter,sans-serif" font-size="9.5" fill="#5f6b7a" text-anchor="end" dominant-baseline="middle">{tn}</text>')
        out.append(f'<text x="{e["x"]-20}" y="{e["y"]+e["h"]-2}" font-family="Inter,sans-serif" font-size="9.5" fill="#5f6b7a" text-anchor="end" dominant-baseline="middle">{bn}</text>')
    out.append("</svg>")
    return "\n".join(out)

path = sys.argv[1] if len(sys.argv)>1 else "fixtures/annotated.cif"
st = load_structure(Path(path).read_text(), Path(path).name)
ch = st.default_chain()
els,res = build_elements(st,ch); annotate_geometry(els,res); cts = strand_contacts(els,res)
from topology.elements import assign_sheets
assign_sheets(els, cts)
META = {e["id"]: e for e in els}
for mode in ["sheet","serpentine","spatial"]:
    L = build_layout(mode, els, res, cts)
    Path(f"preview_{mode}.svg").write_text(svg_for(L, f"{Path(path).name}  chain {ch}  [{mode}]"))
    print(mode, L["extents"])
