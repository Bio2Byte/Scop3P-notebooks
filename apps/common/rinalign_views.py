"""RINAlign view builders.

Ported verbatim from ``notebooks/RINAlign_align_and compare_networks.ipynb`` cell 6.
The JavaScript is copied character-for-character and must stay that way: every view
is one large f-string with hand-doubled ``{{``/``}}`` around JS object literals and
function bodies, so reformatting silently corrupts the browser code. There is a
regression test for exactly that in ``tests/unit/test_rinalign_views.py``.

Two changes only:

``_voila_iframe`` is gone. It existed because ``display(HTML(...))`` does not execute
``<script>`` tags under Voila, and it hand-escaped the document into a ``srcdoc``
attribute. Shiny's ``ui.tags.iframe(srcdoc=...)`` escapes attribute values itself, so
the manual ``html.escape`` would double-escape. :func:`html_document` supplies the
same document shell without it, and the app wraps the result.

Element ids are derived with md5 rather than ``hash()``. ``hash()`` of a ``str`` is
salted per interpreter via ``PYTHONHASHSEED``, so the notebook's ids were
self-consistent within one render but different after every restart, which made these
builders impossible to snapshot-test.

Each script-bearing view must get its **own** iframe, and ``linked_view_html`` must
get exactly **one**. See the module tests and ``docs/use-cases/rinalign.md``.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any, Iterable

from common.vendor import asset_url


def _stable_uid(seed: Iterable[Any]) -> str:
    """Six hex characters derived from the data, stable across interpreter restarts."""
    return hashlib.md5(repr(list(seed)).encode()).hexdigest()[:6]


def html_document(inner_html: str) -> str:
    """Wrap a view fragment in the document shell the notebook's iframe used.

    The escaping ``_voila_iframe`` did is deliberately absent: the caller passes this
    to ``ui.tags.iframe(srcdoc=...)``, which escapes it once, correctly.
    """
    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        "<style>body{margin:0;font-family:sans-serif;}</style></head>"
        "<body>" + inner_html + "</body></html>"
    )


# ── Visualizations ────────────────────────────────────────────

def summary_html(diff, lL, lR):
    nc,nl,ng,nm = len(diff['conserved']),len(diff['lost']),len(diff['gained']),len(diff['mutations'])
    jc = diff['jaccard']
    h = f"<div style='font-family:sans-serif;'><h3 style='margin:12px 0 8px;'>{lL} vs {lR}</h3>"
    h += "<div style='display:flex;gap:10px;flex-wrap:wrap;margin:10px 0;'>"
    _desc = {
        'Conserved': 'Contacts (residue-residue pairs) present in BOTH structures.',
        'Lost': f'Contacts present in {lL} but absent in {lR}.',
        'Gained': f'Contacts present in {lR} but absent in {lL}.',
        'Mutations': 'Positions where the residue IDENTITY differs between the two structures (structure A vs B). NOT a PTM or a disease variant.',
    }
    for bg,fg,lbl,val in [('#e8f5e9','#1b5e20','Conserved',nc),('#ffebee','#b71c1c','Lost',nl),('#e3f2fd','#0d47a1','Gained',ng),('#fff8e1','#e65100','Mutations',nm)]:
        h += f"<div title=\"{_desc.get(lbl,'')}\" style='background:{bg};padding:10px 18px;border-radius:8px;cursor:help;'><div style='font-size:11px;color:{fg};'>{lbl} <span style='opacity:.55;'>&#9432;</span></div><div style='font-size:26px;font-weight:600;color:{fg};'>{val}</div></div>"
    h += f"<div title='Overlap of the two contact sets: shared / union. 1.0 = identical contacts.' style='background:#f3e5f5;padding:10px 18px;border-radius:8px;cursor:help;'><div style='font-size:11px;color:#6a1b9a;'>Jaccard <span style='opacity:.55;'>&#9432;</span></div><div style='font-size:26px;font-weight:600;color:#4a148c;'>{jc:.3f}</div></div>"
    h += "</div>"
    h += "<div style='font-size:11px;color:#777;margin:-2px 0 10px;line-height:1.5;max-width:760px;'><b>Conserved / Lost / Gained</b> count <b>contacts (edges)</b> between residues; <b>Mutations</b> count <b>residues</b> whose identity differs between the two structures (A vs B) &mdash; not PTMs or disease variants. Hover any card for the full definition.</div>"
    if diff['mutations']:
        h += f"<h4>Mutations</h4><table style='border-collapse:collapse;font-size:13px;'><tr style='background:#f5f5f5;'><th style='padding:4px 12px;text-align:center;'>Pos</th><th style='padding:4px 12px;text-align:center;'>{lL}</th><th style='padding:4px 12px;text-align:center;'>{lR}</th></tr>"
        for m in diff['mutations']:
            h += f"<tr><td style='padding:4px 12px;text-align:center;'>{m['position']}</td><td style='padding:4px 12px;text-align:center;'>{m['left']}</td><td style='padding:4px 12px;text-align:center;color:#e65100;font-weight:600;'>{m['right']}</td></tr>"
        h += "</table>"
    changed = [r for r in diff['residue_impact'] if r['net_change']!=0 or r['is_mutation']][:15]
    if changed:
        h += f"<h4>Residues with changes</h4><table style='border-collapse:collapse;font-size:12px;width:100%;'>"
        h += (f"<tr style='background:#f5f5f5;'>"
              f"<th style='text-align:center;padding:4px 10px;'>Pos</th>"
              f"<th style='text-align:center;padding:4px 10px;'>{lL}</th>"
              f"<th style='text-align:center;padding:4px 10px;'>{lR}</th>"
              f"<th style='text-align:center;padding:4px 10px;' title='Number of contacts (degree) at this residue in {lL}'>Deg L</th>"
              f"<th style='text-align:center;padding:4px 10px;' title='Number of contacts (degree) at this residue in {lR}'>Deg R</th>"
              f"<th style='text-align:center;padding:4px 10px;' title='Contacts at this residue in {lL} but not {lR}'>Lost</th>"
              f"<th style='text-align:center;padding:4px 10px;' title='Contacts at this residue in {lR} but not {lL}'>Gained</th></tr>")
        for r in changed:
            rb = 'background:#fff8e1;' if r['is_mutation'] else ''
            h += f"<tr style='{rb}'><td style='text-align:center;padding:3px 10px;'>{r['position']}</td><td style='text-align:center;padding:3px 10px;'>{r['resname_L']}</td><td style='text-align:center;padding:3px 10px;'>{r['resname_R']}</td><td style='text-align:center;padding:3px 10px;'>{r['deg_L']}</td><td style='text-align:center;padding:3px 10px;'>{r['deg_R']}</td><td style='text-align:center;padding:3px 10px;color:#c62828;'>{r['lost'] or ''}</td><td style='text-align:center;padding:3px 10px;color:#1565c0;'>{r['gained'] or ''}</td></tr>"
        h += "</table>"
    h += "</div>"
    return h


def contact_map_html(diff, lL, lR, ptm_pos=None, var_pos=None):
    """Interactive contact map with zoom, pan, hover, and click."""
    positions = sorted(diff['matched_pos'])
    if not positions: return "<p>No overlapping residues.</p>"
    pos_js = json.dumps(positions)
    cons_js = json.dumps([list(e) for e in diff['conserved']])
    lost_js = json.dumps([list(e) for e in diff['lost']])
    gain_js = json.dumps([list(e) for e in diff['gained']])
    # KNOWN GAP, carried over from the notebook unchanged.
    #
    # ``ptm_pos`` and ``var_pos`` are accepted but not drawn. The notebook computed
    # ``ptm_js``, ``var_js`` and a mutation list here and then never interpolated any
    # of them into the script below, so the contact map has no overlay layer at all --
    # only the aligned, force and linked views mark PTMs and variants. The parameters
    # are kept so every view shares one call signature, and so wiring the overlay up
    # later is a change in one place. Adding it means writing new JavaScript, which is
    # a feature rather than part of this migration; see docs/use-cases/rinalign.md.
    #
    # A dead ``rn_map`` loop also stood here, assigning the literal 'UNK' to every
    # position behind an ``if hasattr(diff, '__getitem__'): pass``, and was never
    # read. Residue names come from ``ri_map`` below, which works.
    ri_map = {}
    for r in diff.get('residue_impact',[]):
        ri_map[r['position']] = r.get('resname_L','?')
    rn_js = json.dumps(ri_map)
    uid = 'cm' + _stable_uid(positions[:3])

    return f'''<div style="font-family:sans-serif;max-width:700px;">
    <div style="display:flex;gap:5px;flex-wrap:wrap;margin:0 0 6px;align-items:center;">
      <button onclick="{uid}T('cons',this)" class="cmbtn" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #81C784;background:#e8f5e9;color:#2e7d32;">Conserved</button>
      <button onclick="{uid}T('lost',this)" class="cmbtn" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #E57373;background:#ffebee;color:#c62828;">Lost</button>
      <button onclick="{uid}T('gain',this)" class="cmbtn" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #64B5F6;background:#e3f2fd;color:#1565c0;">Gained</button>
      <span style="margin-left:auto;display:flex;gap:3px;align-items:center;font-size:11px;color:#888;">
        <button onclick="{uid}Z(1.4)" style="width:22px;height:22px;font-size:13px;border-radius:3px;cursor:pointer;border:1px solid #ccc;background:#fff;">+</button>
        <button onclick="{uid}Z(1/1.4)" style="width:22px;height:22px;font-size:13px;border-radius:3px;cursor:pointer;border:1px solid #ccc;background:#fff;">-</button>
        <button onclick="{uid}R()" style="padding:1px 6px;font-size:10px;border-radius:3px;cursor:pointer;border:1px solid #ccc;background:#fff;">Reset</button>
        <span id="{uid}zl">1.0x</span>
      </span>
    </div>
    <div style="position:relative;">
      <canvas id="{uid}c" width="640" height="640" style="display:block;border-radius:6px;border:1px solid #e0e0e0;cursor:crosshair;width:100%;"></canvas>
      <div id="{uid}tip" style="position:absolute;pointer-events:none;background:#fff;border:1px solid #ccc;border-radius:4px;padding:3px 7px;font-size:11px;display:none;z-index:2;white-space:nowrap;"></div>
    </div>
    <div id="{uid}nfo" style="font-size:12px;color:#666;margin:4px 0;min-height:16px;">Scroll to zoom, drag to pan. Hover for details. {lL} (lower) vs {lR} (upper).</div>
    </div>
    <script>
    (function(){{
      var P={pos_js},C={cons_js},L={lost_js},G={gain_js},RN={rn_js};
      var N=P.length,pi={{}};P.forEach(function(p,i){{pi[p]=i;}});
      var sh={{cons:1,lost:1,gain:1}},zm=1,px=0,py=0,dg=false,dsx,dsy,dpx,dpy;
      var cv=document.getElementById('{uid}c'),cx=cv.getContext('2d'),W=cv.width,H=cv.height,mg=55,ms=W-mg*2;
      var cc={{cons:'#81C784',lost:'#E57373',gain:'#64B5F6',mut:'#FF9800'}};
      function tc(i){{return mg+i*(ms/N)*zm+px;}}
      function fc(v){{return Math.floor((v-mg-px)/((ms/N)*zm));}}
      function fcy(v){{return Math.floor((v-mg-py)/((ms/N)*zm));}}
      function dr(){{
        cx.clearRect(0,0,W,H);cx.fillStyle='#fff';cx.fillRect(0,0,W,H);
        var cl=ms/N*zm;cx.save();cx.beginPath();cx.rect(mg,mg,ms,ms);cx.clip();
        for(var i=0;i<N;i++){{var x=tc(i),y=mg+i*cl+py;if(x+cl<mg||x>W-mg||y+cl<mg||y>H-mg)continue;cx.fillStyle='#e8e8e4';cx.fillRect(x,y,Math.max(cl,1),Math.max(cl,1));}}
        function dE(edges,col,tri){{cx.fillStyle=col;edges.forEach(function(e){{var ii=pi[e[0]],jj=pi[e[1]];if(ii===undefined||jj===undefined)return;if(tri!=='lower'){{cx.fillRect(tc(jj),mg+ii*cl+py,Math.max(cl,1),Math.max(cl,1));}}if(tri!=='upper'){{cx.fillRect(tc(ii),mg+jj*cl+py,Math.max(cl,1),Math.max(cl,1));}}}});}}
        if(sh.cons)dE(C,cc.cons,'both');if(sh.lost)dE(L,cc.lost,'lower');if(sh.gain)dE(G,cc.gain,'upper');
        cx.restore();cx.fillStyle='#fff';cx.fillRect(0,0,W,mg);cx.fillRect(0,0,mg,H);cx.fillRect(W-mg,0,mg,H);cx.fillRect(0,H-mg,W,mg);
        cx.strokeStyle='#ccc';cx.lineWidth=0.5;cx.strokeRect(mg,mg,ms,ms);
        var ls=Math.max(1,Math.floor(N/(20*zm)));cx.fillStyle='#666';cx.font='10px monospace';cx.textAlign='center';
        for(var i=0;i<N;i+=ls){{var x=tc(i)+cl/2,y=mg+i*cl+py+cl/2;if(x>mg&&x<W-mg)cx.fillText(P[i],x,mg-3);if(y>mg&&y<H-mg){{cx.textAlign='right';cx.fillText(P[i],mg-3,y+3);cx.textAlign='center';}}}}
        cx.fillStyle='#888';cx.font='10px sans-serif';cx.textAlign='center';cx.fillText('{lR} (upper)',W/2,14);cx.save();cx.translate(12,H/2);cx.rotate(-Math.PI/2);cx.fillText('{lL} (lower)',0,0);cx.restore();
        var lx=mg+4,ly=H-mg+14;cx.font='9px sans-serif';
        [[cc.cons,'Conserved'],[cc.lost,'Lost'],[cc.gain,'Gained']].forEach(function(it){{cx.fillStyle=it[0];cx.fillRect(lx,ly-6,7,7);cx.fillStyle='#666';cx.textAlign='left';cx.fillText(it[1],lx+10,ly);lx+=cx.measureText(it[1]).width+22;}});
      }}
      window['{uid}T']=function(k,b){{sh[k]=sh[k]?0:1;b.style.opacity=sh[k]?1:.3;dr();}};
      window['{uid}Z']=function(f){{var nz=Math.max(.5,Math.min(20,zm*f)),cx2=W/2,cy2=H/2;px=cx2-(cx2-px)*(nz/zm);py=cy2-(cy2-py)*(nz/zm);zm=nz;document.getElementById('{uid}zl').textContent=zm.toFixed(1)+'x';dr();}};
      window['{uid}R']=function(){{zm=1;px=0;py=0;document.getElementById('{uid}zl').textContent='1.0x';dr();}};
      cv.addEventListener('wheel',function(e){{e.preventDefault();var r=cv.getBoundingClientRect(),mx=(e.clientX-r.left)*(W/r.width),my=(e.clientY-r.top)*(H/r.height),f=e.deltaY<0?1.15:1/1.15,nz=Math.max(.5,Math.min(20,zm*f));px=mx-(mx-px)*(nz/zm);py=my-(my-py)*(nz/zm);zm=nz;document.getElementById('{uid}zl').textContent=zm.toFixed(1)+'x';dr();}},{{passive:false}});
      cv.addEventListener('mousedown',function(e){{dg=true;dsx=e.clientX;dsy=e.clientY;dpx=px;dpy=py;cv.style.cursor='grabbing';}});
      window.addEventListener('mousemove',function(e){{
        if(dg){{var r=cv.getBoundingClientRect(),sc=W/r.width;px=dpx+(e.clientX-dsx)*sc;py=dpy+(e.clientY-dsy)*sc;dr();return;}}
        var r=cv.getBoundingClientRect(),mx=(e.clientX-r.left)*(W/r.width),my=(e.clientY-r.top)*(H/r.height);
        var ci=fc(mx),ri=fcy(my),tip=document.getElementById('{uid}tip');
        if(ci>=0&&ci<N&&ri>=0&&ri<N&&mx>mg&&mx<W-mg&&my>mg&&my<H-mg){{
          var pC=P[ci],pR2=P[ri],rnC=RN[pC]||'?',rnR=RN[pR2]||'?';
          var lab=rnR+pR2+' vs '+rnC+pC,type='';var k=Math.min(pR2,pC)+','+Math.max(pR2,pC);
          C.forEach(function(x){{if(x[0]+','+x[1]===k)type=' [conserved]';}});L.forEach(function(x){{if(x[0]+','+x[1]===k)type=' [lost]';}});G.forEach(function(x){{if(x[0]+','+x[1]===k)type=' [gained]';}});
          tip.textContent=lab+type;tip.style.display='block';var tr=cv.getBoundingClientRect();tip.style.left=(e.clientX-tr.left+10)+'px';tip.style.top=(e.clientY-tr.top-22)+'px';
          dr();var cl=ms/N*zm;cx.save();cx.beginPath();cx.rect(mg,mg,ms,ms);cx.clip();cx.strokeStyle='rgba(0,0,0,0.08)';cx.lineWidth=0.5;
          var hx=tc(ci)+cl/2,hy=mg+ri*cl+py+cl/2;cx.beginPath();cx.moveTo(hx,mg);cx.lineTo(hx,H-mg);cx.stroke();cx.beginPath();cx.moveTo(mg,hy);cx.lineTo(W-mg,hy);cx.stroke();cx.restore();
        }}else{{tip.style.display='none';}}
      }});
      window.addEventListener('mouseup',function(){{dg=false;cv.style.cursor='crosshair';}});
      cv.addEventListener('click',function(e){{
        var r=cv.getBoundingClientRect(),mx=(e.clientX-r.left)*(W/r.width),my=(e.clientY-r.top)*(H/r.height);
        var ci=fc(mx),ri=fcy(my);if(ci>=0&&ci<N&&ri>=0&&ri<N){{
          var pC=P[ci],pR2=P[ri],k=Math.min(pR2,pC)+','+Math.max(pR2,pC),t='no contact';
          C.forEach(function(x){{if(x[0]+','+x[1]===k)t='conserved';}});L.forEach(function(x){{if(x[0]+','+x[1]===k)t='lost';}});G.forEach(function(x){{if(x[0]+','+x[1]===k)t='gained';}});
          document.getElementById('{uid}nfo').innerHTML='<b>'+(RN[pR2]||'?')+pR2+'</b> vs <b>'+(RN[pC]||'?')+pC+'</b> — '+t;
        }}
      }});
      dr();
    }})();
    </script>'''


def aligned_network_html(diff, GL, GR, lL, lR, ptm_pos=None, var_pos=None):
    """Interactive aligned network overlay with range slider."""
    all_pos = sorted(set(
        [GL.nodes[n].get('position') for n in GL if GL.nodes[n].get('position') is not None] +
        [GR.nodes[n].get('position') for n in GR if GR.nodes[n].get('position') is not None]
    ))
    pLeft = diff.get('pos_to_left',{})
    pRight = diff.get('pos_to_right',{})
    ptm_set = set(int(x) for x in (ptm_pos or []))
    var_set = set(int(x) for x in (var_pos or []))
    nodes_js = []
    for p in all_pos:
        inA,inB = p in pLeft, p in pRight
        aaA = GL.nodes[pLeft[p]].get('one_letter','?') if inA else '?'
        aaB = GR.nodes[pRight[p]].get('one_letter','?') if inB else '?'
        rnA = GL.nodes[pLeft[p]].get('resname','?') if inA else '?'
        rnB = GR.nodes[pRight[p]].get('resname','?') if inB else '?'
        im = any(m['position']==p for m in diff['mutations'])
        nodes_js.append({'p':p,'aa':aaA if inA else aaB,'rn':rnA,'rnB':rnB if im else None,'inA':inA,'inB':inB,'mut':im,'ptm':(p in ptm_set),'var':(p in var_set)})
    edges_js = {'cons':[list(e) for e in diff['conserved']],'lost':[list(e) for e in diff['lost']],'gain':[list(e) for e in diff['gained']],'onlyA':[list(e) for e in diff.get('onlyA_edges',[])]}
    nc,nl,ng,noa,nm = len(diff['conserved']),len(diff['lost']),len(diff['gained']),len(diff.get('onlyA_edges',[])),len(diff['mutations'])
    mn = min(all_pos) if all_pos else 0
    mx = max(all_pos) if all_pos else 100
    nj = json.dumps(nodes_js)
    ej = json.dumps(edges_js)

    return f'''<div style="font-family:sans-serif;max-width:750px;">
    <h4 style="margin:18px 0 8px;">Aligned network overlay</h4>
    <div style="display:flex;gap:5px;flex-wrap:wrap;margin:0 0 6px;">
      <button id="bc" onclick="tg('cons',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #81C784;background:#e8f5e9;color:#2e7d32;">Conserved ({nc})</button>
      <button id="bl" onclick="tg('lost',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #E57373;background:#ffebee;color:#c62828;">Lost ({nl})</button>
      <button id="bg" onclick="tg('gain',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #64B5F6;background:#e3f2fd;color:#1565c0;">Gained ({ng})</button>
      <button id="bm" title="Residues whose identity differs between the two structures (A vs B) — not a PTM/variant" onclick="tg('mut',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #FFB74D;background:#fff8e1;color:#e65100;">Mutations ({nm})</button>
      <button onclick="tg('ptm',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #E58FD8;background:#FBEAF8;color:#8A1878;">PTMs</button>
      <button onclick="tg('var',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #7FCFCF;background:#E6F7F7;color:#0A6E6E;">Variants</button>
      <button id="bo" onclick="tg('onlyA',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #F48FB1;background:#fce4ec;color:#880e4f;">Only {lL} ({noa})</button>
    </div>
    <div style="display:flex;gap:5px;flex-wrap:wrap;margin:0 0 6px;align-items:center;">
      <span style="font-size:11px;color:#888;">Layout:</span>
      <button onclick="sV('arc')" id="va" style="padding:2px 7px;font-size:10px;border:1px solid #ccc;border-radius:3px;cursor:pointer;background:#e0e0e0;margin-left:10px;">Arc</button>
      <button onclick="sV('linear')" id="vl" style="padding:2px 7px;font-size:10px;border:1px solid #ccc;border-radius:3px;cursor:pointer;">Linear</button>
      <button onclick="sV('circular')" id="vc" style="padding:2px 7px;font-size:10px;border:1px solid #ccc;border-radius:3px;cursor:pointer;">Circular</button>
    </div>
    <div style="display:flex;gap:6px;align-items:center;font-size:12px;color:#666;margin:4px 0;">
      <label>Range:</label>
      <input type="range" id="rlo" min="{mn}" max="{mx}" value="{mn}" style="flex:1;" oninput="uR()">
      <span id="rlov">{mn}</span><span>-</span>
      <input type="range" id="rhi" min="{mn}" max="{mx}" value="{mx}" style="flex:1;" oninput="uR()">
      <span id="rhiv">{mx}</span>
      <button onclick="rR()" style="padding:2px 6px;font-size:10px;border:1px solid #ccc;border-radius:3px;cursor:pointer;">Reset</button>
    </div>
    <div id="ng"></div>
    <div id="ni" style="font-size:12px;color:#666;min-height:32px;margin:4px 0;line-height:1.5;">Click a residue to inspect.</div>
    </div>
    <script>
    (function(){{var N={nj},E={ej},S={{cons:1,lost:1,gain:1,mut:1,onlyA:1,ptm:1,var:1}},vw='arc',sel=null,lo={mn},hi={mx};
    var C={{cons:'#81C784',lost:'#E57373',gain:'#64B5F6',mut:'#FFB74D',onlyA:'#F48FB1',nd:'#B0BEC5'}};
    window.tg=function(k,b){{S[k]=S[k]?0:1;b.style.opacity=S[k]?1:.3;dr();}};
    window.sV=function(v){{vw=v;document.getElementById('va').style.background=v=='arc'?'#e0e0e0':'transparent';document.getElementById('vl').style.background=v=='linear'?'#e0e0e0':'transparent';document.getElementById('vc').style.background=v=='circular'?'#e0e0e0':'transparent';dr();}};
    window.uR=function(){{lo=+document.getElementById('rlo').value;hi=+document.getElementById('rhi').value;if(lo>hi){{var t=lo;lo=hi;hi=t;}}document.getElementById('rlov').textContent=lo;document.getElementById('rhiv').textContent=hi;dr();}};
    window.rR=function(){{lo={mn};hi={mx};document.getElementById('rlo').value={mn};document.getElementById('rhi').value={mx};document.getElementById('rlov').textContent={mn};document.getElementById('rhiv').textContent={mx};dr();}};
    window.cn=function(p){{sel=p;var r=N.find(function(x){{return x.p==p;}});var info=document.getElementById('ni');if(!r)return;var h='<b>'+r.rn+(r.mut?' &rarr; '+r.rnB:'')+'</b> pos '+p;if(r.mut)h+=' <span style="color:#e65100">(mut)</span>';function pr(el){{return el.filter(function(e){{return e[0]==p||e[1]==p;}}).map(function(e){{var q=e[0]==p?e[1]:e[0];var rr=N.find(function(x){{return x.p==q;}});return rr?rr.aa+'<sub>'+q+'</sub>':q;}});}}
    var a=pr(E.cons),b=pr(E.lost),c=pr(E.gain);if(a.length)h+='<br><span style="color:#2e7d32">Conserved: </span>'+a.join(', ');if(b.length)h+='<br><span style="color:#c62828">Lost: </span>'+b.join(', ');if(c.length)h+='<br><span style="color:#1565c0">Gained: </span>'+c.join(', ');info.innerHTML=h;dr();}};
    function dr(){{var fn=N.filter(function(r){{return r.p>=lo&&r.p<=hi;}}),n=fn.length;if(!n){{document.getElementById('ng').innerHTML='<p style="color:#999;">No residues in range.</p>';return;}}
    var pi={{}};fn.forEach(function(r,i){{pi[r.p]=i;}});var W=700,H=vw=='linear'?Math.max(180,Math.min(300,n*2)):Math.max(280,Math.min(420,n*4));var pos=[];
    if(vw=='arc')fn.forEach(function(r,i){{var t=n>1?i/(n-1):.5,a=Math.PI*.15+t*Math.PI*.7;pos.push({{x:W/2-Math.cos(a)*(W-80)*.48+40,y:H-45-Math.sin(a)*(H-70)*.85}});}});
    else if(vw=='linear'){{var g=(W-60)/Math.max(1,n-1);fn.forEach(function(r,i){{pos.push({{x:30+i*g,y:H/2}});}});}}
    else{{var cx=W/2,cy=H/2,rad=Math.min(W,H)/2-40;fn.forEach(function(r,i){{var a=-Math.PI/2+2*Math.PI*i/n;pos.push({{x:cx+rad*Math.cos(a),y:cy+rad*Math.sin(a)}});}});}}
    var svg='<svg width="100%" viewBox="0 0 '+W+' '+H+'" xmlns="http://www.w3.org/2000/svg">';
    function dE(el,col,dash,sw,op){{el.forEach(function(e){{var ii=pi[e[0]],jj=pi[e[1]];if(ii===undefined||jj===undefined)return;var p1=pos[ii],p2=pos[jj],d=dash?'stroke-dasharray="3 2"':'';if(vw=='linear'){{var mx2=(p1.x+p2.x)/2,dist=Math.abs(p2.x-p1.x);svg+='<path d="M'+p1.x+','+p1.y+' Q'+mx2+','+(p1.y-dist*.35)+' '+p2.x+','+p2.y+'" fill="none" stroke="'+col+'" stroke-width="'+sw+'" '+d+' opacity="'+op+'"/>';}}else{{svg+='<line x1="'+p1.x+'" y1="'+p1.y+'" x2="'+p2.x+'" y2="'+p2.y+'" stroke="'+col+'" stroke-width="'+sw+'" '+d+' opacity="'+op+'"/>';}}}});}}
    if(S.onlyA)dE(E.onlyA,C.onlyA,1,.9,.45);if(S.cons)dE(E.cons,C.cons,0,1.2,.55);if(S.lost)dE(E.lost,C.lost,1,2.2,.85);if(S.gain)dE(E.gain,C.gain,0,2.2,.85);
    var ls=Math.max(1,Math.floor(n/28));fn.forEach(function(r,i){{var p=pos[i];if(!r.inB&&!S.onlyA)return;var rd=r.inB?(n>100?5:n>50?7:10):(n>100?3:5);var f=r.inB?C.nd:'transparent',st=C.nd,sw=.5;if(r.mut&&S.mut){{f=C.mut;st='#e65100';sw=1.5;}}if(!r.inB){{st=C.onlyA;sw=1;}}if(sel==r.p){{st='#333';sw=2.5;}}var op=r.inB?1:.5;svg+='<g style="cursor:pointer" onclick="cn('+r.p+')">';svg+='<circle cx="'+p.x+'" cy="'+p.y+'" r="'+rd+'" fill="'+f+'" stroke="'+st+'" stroke-width="'+sw+'" opacity="'+op+'"/>';if(r.ptm&&S.ptm)svg+='<circle cx="'+p.x+'" cy="'+p.y+'" r="'+(rd+3)+'" fill="none" stroke="#C724B1" stroke-width="2"/>';if(r['var']&&S['var'])svg+='<circle cx="'+p.x+'" cy="'+p.y+'" r="'+(rd+5)+'" fill="none" stroke="#0AA5A5" stroke-width="2" stroke-dasharray="2 2"/>';if(n<=55)svg+='<text x="'+p.x+'" y="'+p.y+'" text-anchor="middle" dominant-baseline="central" font-size="'+(n>40?6:8)+'" font-weight="bold" fill="'+(r.mut&&S.mut?'#4a2800':'#fff')+'" opacity="'+op+'">'+r.aa+'</text>';if(i%ls==0||r.mut||sel==r.p)svg+='<text x="'+p.x+'" y="'+(p.y+rd+9)+'" text-anchor="middle" font-size="'+(n>80?5:7)+'" fill="#999">'+r.p+'</text>';svg+='</g>';}});
    svg+='</svg>';document.getElementById('ng').innerHTML=svg;}}
    dr();}})();
    </script>'''



def force_network_html(diff, GL, GR, lL, lR, ptm_pos=None, var_pos=None):
    """D3 force-directed network with proper edge rendering and node coloring."""
    all_pos = sorted(set(
        [GL.nodes[n].get('position') for n in GL if GL.nodes[n].get('position') is not None] +
        [GR.nodes[n].get('position') for n in GR if GR.nodes[n].get('position') is not None]
    ))
    pLeft = diff.get('pos_to_left',{})
    pRight = diff.get('pos_to_right',{})
    ptm_set = set(int(x) for x in (ptm_pos or []))
    var_set = set(int(x) for x in (var_pos or []))
    nodes_js = []
    for p in all_pos:
        inA,inB = p in pLeft, p in pRight
        aaA = GL.nodes[pLeft[p]].get('one_letter','?') if inA else '?'
        rnA = GL.nodes[pLeft[p]].get('resname','?') if inA else '?'
        rnB = GR.nodes[pRight[p]].get('resname','?') if inB else '?'
        im = any(m['position']==p for m in diff['mutations'])
        cat = 'both' if (inA and inB) else ('onlyA' if inA else 'onlyB')
        degA = GL.degree(pLeft[p]) if inA else 0
        degB = GR.degree(pRight[p]) if inB else 0
        lostH = len([e for e in diff['lost'] if p in e])
        gainH = len([e for e in diff['gained'] if p in e])
        nodes_js.append({'id':p,'aa':aaA if inA else (GR.nodes[pRight[p]].get('one_letter','?') if inB else '?'),
                         'rn':rnA if inA else rnB,'rnB':rnB if im else None,
                         'mut':im,'ptm':(p in ptm_set),'var':(p in var_set),'cat':cat,'deg':max(degA,degB),'change':abs(gainH-lostH)})

    nP = sum(1 for n in nodes_js if n['ptm'])
    nV = sum(1 for n in nodes_js if n['var'])

    # Build ONE combined edge array — D3 forceLink will mutate this in place
    all_edges_js = []
    for e in diff['conserved']:
        all_edges_js.append({'source':e[0],'target':e[1],'type':'cons'})
    for e in diff['lost']:
        all_edges_js.append({'source':e[0],'target':e[1],'type':'lost'})
    for e in diff['gained']:
        all_edges_js.append({'source':e[0],'target':e[1],'type':'gain'})

    nc,nl,ng,nm = len(diff['conserved']),len(diff['lost']),len(diff['gained']),len(diff['mutations'])
    nA = len(diff.get('only_left_pos',[]))
    nB = len(diff.get('only_right_pos',[]))
    uid = 'fn' + _stable_uid(all_pos[:3])
    nj = json.dumps(nodes_js)
    aej = json.dumps(all_edges_js)

    return f'''<div style="font-family:sans-serif;max-width:750px;">
    <h4 style="margin:18px 0 8px;">Force-directed network layout</h4>
    <div style="display:flex;gap:5px;flex-wrap:wrap;margin:0 0 6px;align-items:center;">
      <button onclick="{uid}T('cons',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #5DCAA5;background:#E1F5EE;color:#085041;">Conserved ({nc})</button>
      <button onclick="{uid}T('lost',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #F09595;background:#FCEBEB;color:#791F1F;">Lost ({nl})</button>
      <button onclick="{uid}T('gain',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #85B7EB;background:#E6F1FB;color:#0C447C;">Gained ({ng})</button>
      <span style="width:1px;height:16px;background:#ddd;margin:0 3px;display:inline-block;"></span>
      <button onclick="{uid}T('both',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #AFA9EC;background:#EEEDFE;color:#3C3489;">Both models</button>
      <button onclick="{uid}T('onlyA',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #85B7EB;background:#E6F1FB;color:#0C447C;">Only {lL} ({nA})</button>
      <button onclick="{uid}T('onlyB',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #F0997B;background:#FAECE7;color:#712B13;">Only {lR} ({nB})</button>
      <button title="Residues whose identity differs between the two structures (A vs B) — not a PTM or disease variant" onclick="{uid}T('mut',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #FAC775;background:#FAEEDA;color:#633806;">Mutations ({nm})</button>
      <button onclick="{uid}T('ptm',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #E58FD8;background:#FBEAF8;color:#8A1878;">PTMs ({nP})</button>
      <button onclick="{uid}T('var',this)" style="padding:3px 9px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid #7FCFCF;background:#E6F7F7;color:#0A6E6E;">Variants ({nV})</button>
    </div>
    <div style="display:flex;gap:8px;align-items:center;font-size:11px;color:#888;margin:4px 0;flex-wrap:wrap;">
      <label>Repulsion</label><input type="range" id="{uid}ch" min="-500" max="-20" value="-150" step="10" style="width:80px;" oninput="{uid}S('c',this.value)">
      <label>Link dist</label><input type="range" id="{uid}ld" min="15" max="100" value="40" step="5" style="width:80px;" oninput="{uid}S('d',this.value)">
      <label>Node size</label>
      <select id="{uid}sz" onchange="{uid}R()" style="font-size:11px;padding:1px 4px;width:auto;min-width:80px;display:inline-block;"><option value="degree">Degree</option><option value="change">Contact change</option><option value="uniform">Uniform</option></select>
      <label style="display:inline-flex;align-items:center;gap:2px;"><input type="checkbox" id="{uid}lb" checked onchange="{uid}R()"> Labels</label>
      <label style="display:inline-flex;align-items:center;gap:2px;"><input type="checkbox" id="{uid}pn" onchange="{uid}R()"> Pos #</label>
      <span style="margin-left:auto;display:inline-flex;gap:3px;">
        <button onclick="{uid}H()" style="padding:2px 7px;font-size:10px;border-radius:4px;cursor:pointer;border:1px solid #ccc;background:#fff;">Re-layout</button>
        <button onclick="{uid}F()" id="{uid}fb" style="padding:2px 7px;font-size:10px;border-radius:4px;cursor:pointer;border:1px solid #ccc;background:#fff;">Freeze</button>
      </span>
    </div>
    <div id="{uid}g" style="width:100%;min-height:520px;"><p id="{uid}ld2" style="color:#888;font-size:13px;padding:20px;">Loading force layout...</p></div>
    <div style="display:flex;gap:12px;flex-wrap:wrap;font-size:11px;color:#888;padding:6px 10px;background:#f5f5f3;border-radius:6px;margin:4px 0;">
      <span style="font-weight:600;">Nodes:</span>
      <span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:10px;height:10px;border-radius:50%;background:#7F77DD;display:inline-block;"></span> Both</span>
      <span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:10px;height:10px;border-radius:50%;background:#378ADD;display:inline-block;"></span> Only {lL}</span>
      <span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:10px;height:10px;border-radius:50%;background:#D85A30;display:inline-block;"></span> Only {lR}</span>
      <span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:10px;height:10px;border-radius:50%;background:#EF9F27;border:2px solid #854F0B;display:inline-block;"></span> Mutation</span>
      <span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:10px;height:10px;border-radius:50%;background:#fff;border:2px solid #C724B1;display:inline-block;"></span> PTM</span>
      <span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:10px;height:10px;border-radius:50%;background:#fff;border:2px dashed #0AA5A5;display:inline-block;"></span> Variant</span>
      <span style="font-weight:600;margin-left:8px;">Edges:</span>
      <span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:16px;border-top:2px solid #1D9E75;display:inline-block;"></span> Conserved</span>
      <span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:16px;border-top:2px dashed #E24B4A;display:inline-block;"></span> Lost</span>
      <span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:16px;border-top:2px solid #378ADD;display:inline-block;"></span> Gained</span>
    </div>
    <div style="font-size:11px;color:#555;background:#fbfbf9;border:1px solid #eee;border-radius:6px;padding:6px 10px;margin:4px 0;line-height:1.5;"><b>What the categories mean &mdash;</b> <b>Conserved:</b> contact present in <i>both</i> models. <b>Lost:</b> contact in {lL} but not {lR}. <b>Gained:</b> contact in {lR} but not {lL}. <b>Only {lL} / Only {lR} / Both:</b> residue present in that model only, or in both. <b>PTM / Variant rings:</b> residue is an annotated PTM (magenta) or disease variant (teal). Use the buttons above to show/hide each set.</div><div id="{uid}nfo" style="font-size:12px;color:#666;min-height:32px;margin:4px 0;line-height:1.6;">Drag nodes to rearrange. Scroll to zoom. Click a residue to inspect.</div>
    </div>
    <script>
    (function(){{
      function _run(){{
        var el=document.getElementById('{uid}ld2');if(el)el.remove();
        var nodes={nj};
        var allEdges={aej};
        var nMap={{}};nodes.forEach(function(n){{nMap[n.id]=n;}});
        var sh={{cons:1,lost:1,gain:1,both:1,onlyA:1,onlyB:1,mut:1,ptm:1,var:1}},frozen=false,sel=null;
        var NC={{both:'#7F77DD',bothS:'#534AB7',onlyA:'#378ADD',onlyAS:'#185FA5',onlyB:'#D85A30',onlyBS:'#993C1D',mut:'#EF9F27',mutS:'#854F0B',sel:'#333'}};
        var EC={{cons:'#1D9E75',lost:'#E24B4A',gain:'#378ADD'}};
        var W=680,H=520;
        var container=document.getElementById('{uid}g');
        container.innerHTML='';
        var svg=d3.select(container).append('svg').attr('width','100%').attr('viewBox','0 0 '+W+' '+H);
        var gp=svg.append('g');
        var __zoom=d3.zoom().scaleExtent([0.2,10]).on('zoom',function(ev){{gp.attr('transform',ev.transform);}});svg.call(__zoom);

        function nR(d){{var m=document.getElementById('{uid}sz').value;if(m==='degree')return Math.max(5,Math.min(20,3+d.deg*1.2));if(m==='change')return Math.max(5,Math.min(20,5+d.change*3));return 9;}}
        function nFill(d){{if(d.mut&&sh.mut)return NC.mut;if(d.cat==='onlyA')return NC.onlyA;if(d.cat==='onlyB')return NC.onlyB;return NC.both;}}
        function nStroke(d){{if(sel===d.id)return NC.sel;if(d.mut&&sh.mut)return NC.mutS;if(d.cat==='onlyA')return NC.onlyAS;if(d.cat==='onlyB')return NC.onlyBS;return NC.bothS;}}

        var sim=d3.forceSimulation(nodes)
          .force('link',d3.forceLink(allEdges).id(function(d){{return d.id;}}).distance(40).strength(0.35))
          .force('charge',d3.forceManyBody().strength(-150))
          .force('center',d3.forceCenter(W/2,H/2))
          .force('collision',d3.forceCollide().radius(function(d){{return nR(d)+3;}}))
          .force('x',d3.forceX(W/2).strength(0.03))
          .force('y',d3.forceY(H/2).strength(0.03))
          .alphaDecay(0.018)
          .on('tick',tick);

        var lkG=gp.append('g');
        var riG=gp.append('g');
        var ptmG=gp.append('g');
        var varG=gp.append('g');
        var ndG=gp.append('g');
        var lbG=gp.append('g');
        var psG=gp.append('g');

        function tick(){{
          var ve=allEdges.filter(function(e){{return sh[e.type];}});
          var vn=nodes.filter(function(n){{if(n.cat==='onlyA'&&!sh.onlyA)return false;if(n.cat==='onlyB'&&!sh.onlyB)return false;if(n.cat==='both'&&!sh.both)return false;return true;}});

          var ln=lkG.selectAll('line').data(ve,function(d){{return d.source.id+'-'+d.target.id+'-'+d.type;}});
          ln.exit().remove();
          ln=ln.enter().append('line').merge(ln);
          ln.attr('x1',function(d){{return d.source.x;}})
            .attr('y1',function(d){{return d.source.y;}})
            .attr('x2',function(d){{return d.target.x;}})
            .attr('y2',function(d){{return d.target.y;}})
            .attr('stroke',function(d){{return EC[d.type];}})
            .attr('stroke-width',function(d){{return d.type==='cons'?0.8:2.5;}})
            .attr('stroke-dasharray',function(d){{return d.type==='lost'?'5 3':null;}})
            .attr('opacity',function(d){{return d.type==='cons'?0.25:0.65;}});

          var ri=riG.selectAll('circle').data(vn.filter(function(d){{return d.mut&&sh.mut;}}),function(d){{return 'r'+d.id;}});
          ri.exit().remove();
          ri=ri.enter().append('circle').merge(ri);
          ri.attr('cx',function(d){{return d.x;}}).attr('cy',function(d){{return d.y;}})
            .attr('r',function(d){{return nR(d)+4;}})
            .attr('fill','none').attr('stroke','#EF9F27').attr('stroke-width',2.5)
            .attr('stroke-dasharray','3 2').attr('opacity',0.8);

          var pr=ptmG.selectAll('circle').data(vn.filter(function(d){{return d.ptm&&sh.ptm;}}),function(d){{return 'p'+d.id;}});
          pr.exit().remove();
          pr=pr.enter().append('circle').merge(pr);
          pr.attr('cx',function(d){{return d.x;}}).attr('cy',function(d){{return d.y;}})
            .attr('r',function(d){{return nR(d)+7;}})
            .attr('fill','none').attr('stroke','#C724B1').attr('stroke-width',2.5).attr('opacity',0.9);

          var vr=varG.selectAll('circle').data(vn.filter(function(d){{return d.var&&sh.var;}}),function(d){{return 'v'+d.id;}});
          vr.exit().remove();
          vr=vr.enter().append('circle').merge(vr);
          vr.attr('cx',function(d){{return d.x;}}).attr('cy',function(d){{return d.y;}})
            .attr('r',function(d){{return nR(d)+10;}})
            .attr('fill','none').attr('stroke','#0AA5A5').attr('stroke-width',2.5).attr('stroke-dasharray','2 3').attr('opacity',0.9);

          var ci=ndG.selectAll('circle').data(vn,function(d){{return d.id;}});
          ci.exit().remove();
          ci=ci.enter().append('circle')
            .call(d3.drag()
              .on('start',function(ev,d){{if(!ev.active)sim.alphaTarget(0.3).restart();d.fx=d.x;d.fy=d.y;}})
              .on('drag',function(ev,d){{d.fx=ev.x;d.fy=ev.y;}})
              .on('end',function(ev,d){{if(!ev.active)sim.alphaTarget(0);if(!frozen){{d.fx=null;d.fy=null;}}}}))
            .on('click',function(ev,d){{cN(d);if(window.__RIN_ONSELECT)window.__RIN_ONSELECT(d.id);}})
            .merge(ci);
          ci.attr('cx',function(d){{return d.x;}}).attr('cy',function(d){{return d.y;}})
            .attr('r',function(d){{return nR(d);}})
            .attr('fill',function(d){{return nFill(d);}})
            .attr('stroke',function(d){{return nStroke(d);}})
            .attr('stroke-width',function(d){{return sel===d.id?2.5:(d.mut&&sh.mut?2:0.5);}})
            .attr('cursor','grab');

          var sl=document.getElementById('{uid}lb').checked;
          var tx=lbG.selectAll('text').data(sl?vn:[],function(d){{return d.id;}});
          tx.exit().remove();
          tx=tx.enter().append('text').attr('text-anchor','middle').attr('dominant-baseline','central')
            .attr('font-family','sans-serif').attr('pointer-events','none').merge(tx);
          tx.attr('x',function(d){{return d.x;}}).attr('y',function(d){{return d.y;}})
            .attr('font-size',function(d){{return Math.max(7,nR(d)-2)+'px';}})
            .attr('font-weight','500')
            .attr('fill',function(d){{return d.mut&&sh.mut?'#412402':'#fff';}})
            .text(function(d){{return d.aa;}});

          var sp=document.getElementById('{uid}pn').checked;
          var px=psG.selectAll('text').data(sp?vn:[],function(d){{return 'p'+d.id;}});
          px.exit().remove();
          px=px.enter().append('text').attr('text-anchor','middle')
            .attr('font-family','monospace').attr('pointer-events','none')
            .attr('font-size','7px').attr('fill','#999').merge(px);
          px.attr('x',function(d){{return d.x;}}).attr('y',function(d){{return d.y+nR(d)+10;}})
            .text(function(d){{return d.id;}});
        }}

        function cN(d){{
          sel=d.id;var info=document.getElementById('{uid}nfo');
          var catL={{both:'both models',onlyA:'only {lL}',onlyB:'only {lR}'}};
          var h='<b>'+d.rn+(d.mut?' &rarr; '+d.rnB:'')+'</b> pos '+d.id+' ('+catL[d.cat]+', deg '+d.deg+')';
          if(d.mut)h+=' <span style="color:#EF9F27">mutation</span>';
          function pr(type){{
            return allEdges.filter(function(e){{return e.type===type&&(e.source.id===d.id||e.target.id===d.id);}})
              .map(function(e){{var q=e.source.id===d.id?e.target.id:e.source.id;var r=nMap[q];
                var tag=r?({{both:'',onlyA:' [A]',onlyB:' [B]'}})[r.cat]:'';
                return r?r.aa+'<sub>'+q+'</sub>'+tag:q;}});
          }}
          var a=pr('cons'),b=pr('lost'),c=pr('gain');
          if(a.length)h+='<br><span style="color:#1D9E75">Conserved ('+a.length+'):</span> '+a.join(', ');
          if(b.length)h+='<br><span style="color:#E24B4A">Lost ('+b.length+'):</span> '+b.join(', ');
          if(c.length)h+='<br><span style="color:#378ADD">Gained ('+c.length+'):</span> '+c.join(', ');
          info.innerHTML=h;tick();
        }}
        window['{uid}HL']=function(p){{var nd=nMap[p];if(!nd)return;cN(nd);var k=2.5;try{{var t=d3.zoomIdentity.translate(W/2-nd.x*k,H/2-nd.y*k).scale(k);svg.transition().duration(500).call(__zoom.transform,t);}}catch(e){{}}}};
        window.__RIN_HL=window['{uid}HL'];

        window['{uid}T']=function(k,b){{
          sh[k]=sh[k]?0:1;b.style.opacity=sh[k]?1:0.3;
          sim.force('link',d3.forceLink(allEdges.filter(function(e){{return sh[e.type];}})).id(function(d){{return d.id;}}).distance(+document.getElementById('{uid}ld').value).strength(0.35));
          sim.alpha(0.3).restart();
        }};
        window['{uid}R']=function(){{tick();}};
        window['{uid}S']=function(w,v){{
          if(w==='c')sim.force('charge').strength(+v);
          if(w==='d')sim.force('link').distance(+v);
          sim.alpha(0.3).restart();
        }};
        window['{uid}H']=function(){{nodes.forEach(function(d){{d.fx=null;d.fy=null;}});frozen=false;document.getElementById('{uid}fb').textContent='Freeze';sim.alpha(1).restart();}};
        window['{uid}F']=function(){{frozen=!frozen;document.getElementById('{uid}fb').textContent=frozen?'Unfreeze':'Freeze';if(frozen){{sim.stop();nodes.forEach(function(d){{d.fx=d.x;d.fy=d.y;}});}}else{{nodes.forEach(function(d){{d.fx=null;d.fy=null;}});sim.alpha(0.3).restart();}}}};
      }}
      if(typeof d3!=='undefined'){{setTimeout(_run,50);}}
      else if(typeof require!=='undefined'){{try{{require.config({{paths:{{d3:'{asset_url("d3").removesuffix(".js")}'}}}});require(['d3'],function(d){{window.d3=d;_run();}});}}catch(e){{var s=document.createElement('script');s.src='{asset_url("d3")}';s.onload=_run;document.head.appendChild(s);}}}}
      else{{var s=document.createElement('script');s.src='{asset_url("d3")}';s.onload=_run;document.head.appendChild(s);}}
    }})();
    </script>'''


def linked_view_html(diff, GL, GR, lL, lR, pdb_text, pdb_fmt='pdb', chain=None, ptm_pos=None, var_pos=None):
    """Linked view that REUSES the real force_network_html on the left and puts NGL on
    the right, bridged two ways. A residue dropdown highlights in both. RIN 'position'
    == NGL author residue number for the left structure (no SIFTS needed)."""
    pL = diff.get('pos_to_left', {}); pR = diff.get('pos_to_right', {})
    all_pos = sorted(set(pL) | set(pR)); in_left = set(pL)
    opts = ["<option value=''>-- go to residue --</option>"]
    for p in all_pos:
        rn = (GL.nodes[pL[p]].get('resname', '') if p in pL else GR.nodes[pR[p]].get('resname', ''))
        tag = '' if p in in_left else ' (not in 3D)'
        opts.append(f"<option value='{p}'>{rn} {p}{tag}</option>")
    opts_html = ''.join(opts)
    lostBy = {}; gainBy = {}; consBy = {}
    for e in diff['lost']:
        lostBy.setdefault(e[0], []).append(e[1]); lostBy.setdefault(e[1], []).append(e[0])
    for e in diff['gained']:
        gainBy.setdefault(e[0], []).append(e[1]); gainBy.setdefault(e[1], []).append(e[0])
    for e in diff['conserved']:
        consBy.setdefault(e[0], []).append(e[1]); consBy.setdefault(e[1], []).append(e[0])
    ext = 'cif' if str(pdb_fmt).lower().startswith('cif') else 'pdb'
    ch = (chain or '').strip()
    ngl_data = json.dumps({'chain': ch, 'ext': ext, 'lostBy': lostBy, 'gainBy': gainBy, 'consBy': consBy, 'inL': sorted(in_left)})
    pdb_js = json.dumps(pdb_text)
    left_html = force_network_html(diff, GL, GR, lL, lR, ptm_pos=ptm_pos, var_pos=var_pos)

    head = ("<div style='font-family:sans-serif;'>"
            "<div style='display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin:8px 0;'>"
            "<b style='font-size:13px;'>Linked network &harr; structure</b>"
            "<span style='color:#888;font-size:11px;'>click a residue in either panel, or use the dropdown</span>"
            "<select id='LVSEARCH' onchange='LVgo(this.value)' style='font-size:12px;padding:3px 6px;border-radius:5px;'>"
            + opts_html + "</select>"
            "<span id='LVinfo' style='margin-left:auto;color:#333;font-size:12px;'></span></div>"
            "<div style='display:flex;gap:10px;align-items:stretch;flex-wrap:wrap;'>"
            "<div style='flex:1 1 430px;min-width:340px;'>")
    mid = ("</div><div style='flex:1 1 380px;min-width:340px;'>"
           "<div id='LVngl' style='width:100%;height:540px;border:1px solid #ddd;border-radius:6px;background:#fff;'></div>"
           "</div></div>"
           "<div style='font-size:11px;color:#777;background:#f5f5f3;border-radius:6px;padding:5px 8px;margin-top:6px;'>"
           "Left = force network (same as the standalone view). Right = structure (grey). "
           "Selected residue = orange; its lost contacts = dashed red, gained = solid blue on the 3D.</div>")
    script = f'''<script src="{asset_url('ngl')}"></script>
    <script>
    (function(){{
      var DAT={ngl_data}, PDB={pdb_js};
      var CH=DAT.chain, INL=new Set(DAT.inL), lostBy=DAT.lostBy, gainBy=DAT.gainBy, consBy=DAT.consBy;
      var comp=null,baseRep=null,selRep=null,lblRep=null,consLn=null,lostLn=null,gainLn=null,nbLbl=null,stage=null;
      function sfx(p){{return CH?(p+':'+CH):(''+p);}}
      function atom(p){{return sfx(p)+'.CA';}}
      function waitNGL(cb,n){{n=n||0;if(window.NGL){{cb();}}else if(n<200){{setTimeout(function(){{waitNGL(cb,n+1);}},50);}}else{{document.getElementById('LVngl').innerHTML="<p style='color:#a00;padding:20px;'>NGL could not load.</p>";}}}}
      waitNGL(function(){{
        try{{
          stage=new NGL.Stage('LVngl',{{backgroundColor:'white'}});
          window.addEventListener('resize',function(){{try{{stage.handleResize();}}catch(e){{}}}});
          stage.loadFile(new Blob([PDB],{{type:'text/plain'}}),{{ext:DAT.ext}}).then(function(o){{
            comp=o;baseRep=comp.addRepresentation('cartoon',{{color:'#c6cad2'}});comp.autoView();
            stage.signals.clicked.add(function(pp){{if(pp&&pp.atom){{var p=pp.atom.resno;nglHi(p);if(window.__RIN_HL)window.__RIN_HL(p);}}}});
          }});
        }}catch(e){{document.getElementById('LVngl').innerHTML="<p style='color:#a00;padding:20px;'>NGL error: "+e+"</p>";}}
      }});
      function nsel(arr){{arr=(arr||[]).filter(function(q){{return INL.has(q);}});return arr.length?arr.map(function(q){{return sfx(q);}}).join(' or '):'none';}}
      function nglHi(p){{
        p=+p;
        var nc=(consBy[p]||[]).filter(function(q){{return INL.has(q);}}).length;
        var nl=(lostBy[p]||[]).filter(function(q){{return INL.has(q);}}).length;
        var ng=(gainBy[p]||[]).filter(function(q){{return INL.has(q);}}).length;
        document.getElementById('LVinfo').textContent='Residue '+p+(INL.has(p)?(' \u00b7 neighbours: '+nc+' conserved, '+nl+' lost, '+ng+' gained'):' (not in this structure)');
        if(!comp||!INL.has(p))return;
        // recolour the cartoon: selected = orange, first-degree neighbours by edge type
        try{{
          var scheme=NGL.ColormakerRegistry.addSelectionScheme([
            ['#FF7F0E', sfx(p)],
            ['#E24B4A', nsel(lostBy[p])],
            ['#1f77b4', nsel(gainBy[p])],
            ['#1D9E75', nsel(consBy[p])],
            ['#c6cad2', '*']
          ]);
          if(baseRep)comp.removeRepresentation(baseRep);
          baseRep=comp.addRepresentation('cartoon',{{colorScheme:scheme}});
        }}catch(e){{}}
        // orange stick marker + label on the selected residue
        if(selRep)comp.removeRepresentation(selRep);
        selRep=comp.addRepresentation('licorice',{{sele:sfx(p),color:'#FF7F0E'}});
        if(lblRep){{comp.removeRepresentation(lblRep);lblRep=null;}}
        try{{lblRep=comp.addRepresentation('label',{{sele:sfx(p)+' and .CA',labelType:'residue',color:'#111111',scale:1.4,showBackground:true,backgroundColor:'white',backgroundOpacity:0.55}});}}catch(e){{}}
        // dashed Ca-Ca lines to first-degree neighbours, coloured by edge type
        [consLn,lostLn,gainLn,nbLbl].forEach(function(r){{if(r)comp.removeRepresentation(r);}});
        consLn=lostLn=gainLn=nbLbl=null;
        function pairs(arr){{return (arr||[]).filter(function(q){{return INL.has(q);}}).map(function(q){{return [atom(p),atom(q)];}});}}
        var lc=pairs(consBy[p]),ll=pairs(lostBy[p]),lg=pairs(gainBy[p]);
        try{{if(lc.length)consLn=comp.addRepresentation('distance',{{atomPair:lc,color:'#1D9E75',labelVisible:false,useCylinder:false}});}}catch(e){{}}
        try{{if(ll.length)lostLn=comp.addRepresentation('distance',{{atomPair:ll,color:'#E24B4A',labelVisible:false,useCylinder:false}});}}catch(e){{}}
        try{{if(lg.length)gainLn=comp.addRepresentation('distance',{{atomPair:lg,color:'#1f77b4',labelVisible:false,useCylinder:false}});}}catch(e){{}}
        // labels on the neighbour amino acids
        var allNb=[].concat(consBy[p]||[],lostBy[p]||[],gainBy[p]||[]).filter(function(q){{return INL.has(q);}});
        if(allNb.length){{try{{var s='('+allNb.map(function(q){{return sfx(q);}}).join(' or ')+') and .CA';nbLbl=comp.addRepresentation('label',{{sele:s,labelType:'residue',color:'#444444',scale:1.05,showBackground:true,backgroundColor:'white',backgroundOpacity:0.4}});}}catch(e){{}}}}
        try{{var lo=p-12,hi=p+12;var win=CH?(lo+'-'+hi+':'+CH):(lo+'-'+hi);comp.autoView(win,600);}}catch(e){{try{{comp.autoView(sfx(p),600);}}catch(e2){{}}}}
      }}
      window.__RIN_ONSELECT=function(p){{nglHi(p);}};
      window.LVgo=function(v){{if(v===''||v===null)return;var p=+v;nglHi(p);if(window.__RIN_HL)window.__RIN_HL(p);}};
    }})();
    </script></div>'''
    return head + left_html + mid + script
