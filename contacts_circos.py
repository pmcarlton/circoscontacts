#!/usr/bin/env python3
"""
Generate an interactive circos-style contact plot (HTML + SVG export)
from ChimeraX contact files and corresponding mmCIFs.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


try:
    import gemmi  # type: ignore
except Exception as exc:  # pragma: no cover - runtime check
    raise SystemExit(
        "Missing dependency: gemmi. Install with `pip install gemmi`."
    ) from exc


DNA_RESNAMES = {"DA", "DC", "DG", "DT", "DI", "DU"}


def find_atom_site_category(doc: "gemmi.cif.Document") -> Dict[str, List[str]] | None:
    for block in doc:
        cat = block.get_mmcif_category("_atom_site")
        if cat:
            return cat
    return None


def parse_cifs(cif_paths: Iterable[Path]) -> Dict[str, Dict[int, str]]:
    chain_res: Dict[str, Dict[int, str]] = defaultdict(dict)
    conflicts: List[str] = []
    for path in cif_paths:
        doc = gemmi.cif.read(str(path))
        cat = find_atom_site_category(doc)
        if not cat:
            continue
        asym_key = "auth_asym_id" if "auth_asym_id" in cat else "label_asym_id"
        seq_key = "auth_seq_id" if "auth_seq_id" in cat else "label_seq_id"
        comp_key = "label_comp_id"
        asym = cat.get(asym_key, [])
        seq = cat.get(seq_key, [])
        comp = cat.get(comp_key, [])
        for a, s, c in zip(asym, seq, comp):
            if s in (".", "?"):
                continue
            try:
                resnum = int(s)
            except ValueError:
                continue
            existing = chain_res[a].get(resnum)
            if existing and existing != c:
                conflicts.append(
                    f"{path.name}: chain {a} residue {resnum} had {existing} vs {c}"
                )
                continue
            chain_res[a][resnum] = c
    if conflicts:
        sys.stderr.write(
            "Warning: residue name conflicts across CIFs (showing first 5):\n"
        )
        for line in conflicts[:5]:
            sys.stderr.write(f"  {line}\n")
    return chain_res


def detect_dna_chains(chain_res: Dict[str, Dict[int, str]]) -> List[str]:
    dna_chains = []
    for chain_id, residues in chain_res.items():
        if not residues:
            continue
        dna_count = sum(1 for r in residues.values() if r in DNA_RESNAMES)
        if dna_count / max(1, len(residues)) >= 0.5:
            dna_chains.append(chain_id)
    return sorted(dna_chains)


def parse_contacts(contact_paths: Iterable[Path]) -> Counter:
    counts: Counter = Counter()
    for path in contact_paths:
        with path.open() as fh:
            for line in fh:
                if not line.startswith("/"):
                    continue
                parts = line.split()
                if len(parts) < 8:
                    continue
                chain1 = parts[0].lstrip("/")
                chain2 = parts[4].lstrip("/")
                try:
                    res1 = int(parts[2])
                    res2 = int(parts[6])
                except ValueError:
                    continue
                key = canonical_contact(chain1, res1, chain2, res2)
                counts[key] += 1
    return counts


def canonical_contact(chain1: str, res1: int, chain2: str, res2: int) -> Tuple[str, int, str, int]:
    if (chain1, res1) <= (chain2, res2):
        return chain1, res1, chain2, res2
    return chain2, res2, chain1, res1


def build_chain_maps(
    chain_res: Dict[str, Dict[int, str]],
    dna_chains: List[str],
    dna_reverse: str | None,
) -> Tuple[
    Dict[str, Dict[int, int]],
    Dict[str, List[Dict[str, str | int]]],
    Dict[str, str],
    List[str],
    Dict[str, int],
]:
    """
    Returns:
        pos_map: chain -> resnum -> pos (1-based)
        pos_info: display_chain -> list of info dicts (index 0 => pos 1)
        display_chain_of: chain -> display_chain
        display_chains: list of display chain ids
    """
    pos_map: Dict[str, Dict[int, int]] = {}
    pos_info: Dict[str, List[Dict[str, str | int]]] = {}
    display_chain_of: Dict[str, str] = {}
    start_pos_map: Dict[str, int] = {}

    display_chains: List[str] = []
    # Proteins first (alphabetical), DNA will be added as one chain "DNA".
    protein_chains = [c for c in sorted(chain_res) if c not in dna_chains]
    display_chains.extend(protein_chains)

    # Build protein maps
    for chain_id in protein_chains:
        residues = chain_res.get(chain_id, {})
        resnums = sorted(residues)
        pos_map[chain_id] = {resnum: idx + 1 for idx, resnum in enumerate(resnums)}
        info = [
            {
                "chain": chain_id,
                "resnum": resnum,
                "resname": residues[resnum],
                "source_chain": chain_id,
            }
            for resnum in resnums
        ]
        pos_info[chain_id] = info
        display_chain_of[chain_id] = chain_id
        start_pos_map[chain_id] = 1

    if not dna_chains:
        return pos_map, pos_info, display_chain_of, display_chains, start_pos_map

    display_chains.insert(0, "DNA")
    display_chain_of.update({c: "DNA" for c in dna_chains})

    if len(dna_chains) == 1:
        dna_chain = dna_chains[0]
        residues = chain_res.get(dna_chain, {})
        resnums = sorted(residues)
        pos_map[dna_chain] = {resnum: idx + 1 for idx, resnum in enumerate(resnums)}
        pos_info["DNA"] = [
            {
                "chain": "DNA",
                "resnum": resnum,
                "resname": residues[resnum],
                "source_chain": dna_chain,
            }
            for resnum in resnums
        ]
        start_pos_map["DNA"] = 1
        return pos_map, pos_info, display_chain_of, display_chains, start_pos_map

    # Duplex handling: collapse to base-pair positions.
    dna_rev = dna_reverse
    if dna_rev and dna_rev not in dna_chains:
        sys.stderr.write(
            f"Warning: --dna-reverse {dna_rev} not in detected DNA chains; ignoring.\n"
        )
        dna_rev = None

    dna_forward = None
    if dna_rev:
        dna_forward = next(c for c in dna_chains if c != dna_rev)
    else:
        if len(dna_chains) > 2:
            sys.stderr.write(
                "Warning: more than 2 DNA chains detected; using first two.\n"
            )
        dna_rev = dna_chains[0]
        dna_forward = dna_chains[1]

    f_residues = chain_res.get(dna_forward, {})
    r_residues = chain_res.get(dna_rev, {})
    f_resnums = sorted(f_residues)
    r_resnums = sorted(r_residues)

    length = min(len(f_resnums), len(r_resnums))
    if len(f_resnums) != len(r_resnums):
        sys.stderr.write(
            "Warning: DNA strands have different lengths; "
            f"using min length {length}.\n"
        )

    pos_map[dna_forward] = {}
    pos_map[dna_rev] = {}
    dna_info: List[Dict[str, str | int]] = []
    for i in range(length):
        f_resnum = f_resnums[i]
        r_resnum = r_resnums[-(i + 1)]
        pos = i + 1
        pos_map[dna_forward][f_resnum] = pos
        pos_map[dna_rev][r_resnum] = pos
        dna_info.append(
            {
                "chain": "DNA",
                "resnum": pos,
                "resname": "bp",
                "source_chain": "DNA",
                "pair": {
                    "f_chain": dna_forward,
                    "f_resnum": f_resnum,
                    "f_resname": f_residues.get(f_resnum, "?"),
                    "r_chain": dna_rev,
                    "r_resnum": r_resnum,
                    "r_resname": r_residues.get(r_resnum, "?"),
                },
            }
        )
    pos_info["DNA"] = dna_info
    dna_first = sorted(dna_chains)[0]
    start_pos_map["DNA"] = 1 if dna_first == dna_forward else length

    return pos_map, pos_info, display_chain_of, display_chains, start_pos_map


def aggregate_contacts(
    counts: Counter,
    pos_map: Dict[str, Dict[int, int]],
    display_chain_of: Dict[str, str],
) -> Tuple[List[Dict[str, int | str]], int]:
    agg: Counter = Counter()
    skipped = 0
    for chain1, res1, chain2, res2 in counts:
        if chain1 not in pos_map or chain2 not in pos_map:
            skipped += 1
            continue
        if res1 not in pos_map[chain1] or res2 not in pos_map[chain2]:
            skipped += 1
            continue
        a_chain = display_chain_of.get(chain1, chain1)
        b_chain = display_chain_of.get(chain2, chain2)
        a_pos = pos_map[chain1][res1]
        b_pos = pos_map[chain2][res2]
        key = canonical_contact(a_chain, a_pos, b_chain, b_pos)
        agg[key] += counts[(chain1, res1, chain2, res2)]
    if skipped:
        sys.stderr.write(f"Warning: skipped {skipped} contacts without CIF mapping\n")

    max_count = max(agg.values()) if agg else 1
    contacts = [
        {"a": a, "a_pos": pa, "b": b, "b_pos": pb, "count": c}
        for (a, pa, b, pb), c in agg.items()
    ]
    contacts.sort(key=lambda x: x["count"], reverse=True)
    return contacts, max_count


def build_colors(chains: List[str]) -> Dict[str, str]:
    palette = [
        "#1b9e77",
        "#d95f02",
        "#7570b3",
        "#e7298a",
        "#66a61e",
        "#e6ab02",
        "#a6761d",
        "#666666",
        "#1f78b4",
        "#b2df8a",
        "#fb9a99",
        "#fdbf6f",
    ]
    colors = {}
    idx = 0
    for chain in chains:
        if chain == "DNA":
            colors[chain] = "#6d4c41"
            continue
        colors[chain] = palette[idx % len(palette)]
        idx += 1
    return colors


def generate_html(
    title: str,
    chains: List[str],
    chain_lengths: Dict[str, int],
    chain_colors: Dict[str, str],
    pos_info: Dict[str, List[Dict[str, str | int]]],
    contacts: List[Dict[str, int | str]],
    max_count: int,
    default_order: List[str],
    start_pos_map: Dict[str, int],
    output_path: Path,
) -> None:
    data = {
        "title": title,
        "chains": [
            {
                "id": chain,
                "length": chain_lengths[chain],
                "color": chain_colors[chain],
                "resinfo": pos_info.get(chain, []),
                "start_pos": start_pos_map.get(chain, 1),
            }
            for chain in chains
        ],
        "contacts": contacts,
        "max_count": max_count,
        "default_order": default_order,
    }
    data_json = json.dumps(data, separators=(",", ":"))

    template = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@300;400;600&family=IBM+Plex+Serif:wght@400;600&display=swap');
html, body {
  height: 100%;
  margin: 0;
}
body {
  font-family: 'IBM Plex Sans', sans-serif;
  color: #0f172a;
  background: radial-gradient(circle at 20% 20%, #f1f5f9, #e2e8f0 40%, #cbd5f5 100%);
  overflow: hidden;
}
.page {
  display: grid;
  grid-template-columns: minmax(280px, 360px) 1fr;
  gap: 24px;
  padding: 24px;
  box-sizing: border-box;
  height: 100vh;
  align-items: stretch;
}
.panel {
  background: rgba(255,255,255,0.9);
  border-radius: 16px;
  padding: 20px;
  box-shadow: 0 12px 40px rgba(15, 23, 42, 0.12);
  max-height: calc(100vh - 48px);
}
.controls {
  position: sticky;
  top: 24px;
  overflow: auto;
}
h1 {
  margin: 0 0 8px 0;
  font-family: 'IBM Plex Serif', serif;
  font-weight: 600;
  font-size: 22px;
}
.subtitle {
  color: #475569;
  font-size: 14px;
  margin-bottom: 16px;
}
.control {
  margin-bottom: 14px;
}
.control label {
  display: block;
  font-size: 13px;
  margin-bottom: 6px;
  color: #334155;
}
.row {
  display: flex;
  gap: 8px;
  align-items: center;
}
input[type="range"] {
  width: 100%;
}
input[type="number"], input[type="text"] {
  width: 100%;
  padding: 8px 10px;
  border-radius: 8px;
  border: 1px solid #cbd5f5;
  font-size: 13px;
}
button {
  background: #0f172a;
  color: #fff;
  border: none;
  border-radius: 8px;
  padding: 8px 12px;
  cursor: pointer;
  font-size: 13px;
}
button.secondary {
  background: #64748b;
}
.hint {
  font-size: 12px;
  color: #64748b;
}
.hint .hidden-chain {
  color: #dc2626;
  font-weight: 600;
}
.hint .error {
  color: #b91c1c;
  font-weight: 600;
}
.plot-wrap {
  display: flex;
  align-items: center;
  justify-content: center;
  background: rgba(255,255,255,0.85);
  border-radius: 16px;
  padding: 10px;
  box-shadow: inset 0 0 0 1px rgba(148,163,184,0.3);
  height: calc(100vh - 48px);
  overflow: auto;
}
svg {
  display: block;
}
.tooltip {
  position: absolute;
  pointer-events: none;
  background: rgba(15, 23, 42, 0.92);
  color: #fff;
  padding: 8px 10px;
  border-radius: 8px;
  font-size: 12px;
  opacity: 0;
  transition: opacity 0.15s ease;
}
.status {
  font-size: 12px;
  color: #0f172a;
  background: #e2e8f0;
  padding: 6px 8px;
  border-radius: 8px;
  margin-top: 8px;
}
@media (max-width: 960px) {
  .page {
    grid-template-columns: 1fr;
    height: auto;
  }
  body {
    overflow: auto;
  }
  .controls {
    position: static;
    max-height: none;
  }
  .plot-wrap {
    height: auto;
  }
}
</style>
</head>
<body>
<div class="page">
  <div class="panel controls">
    <h1>__TITLE__</h1>
    <div class="subtitle">Interactive circos view of aggregated ChimeraX contacts.</div>
    <div class="control">
      <label for="threshold">Minimum contact visibility threshold</label>
      <div class="row">
        <input id="threshold" type="range" min="0" max="__MAX_COUNT__" step="1" value="2">
        <input id="thresholdInput" type="number" min="0" max="__MAX_COUNT__" step="1" value="2">
      </div>
    </div>
    <div class="control">
      <label for="maxWidth">Maximum line width</label>
      <div class="row">
        <input id="maxWidth" type="range" min="1" max="12" step="0.5" value="6">
        <input id="maxWidthInput" type="number" min="1" max="12" step="0.5" value="6">
      </div>
    </div>
    <div class="control">
      <label for="angle">Angle offset (degrees)</label>
      <div class="row">
        <input id="angle" type="range" min="0" max="360" step="1" value="0">
        <input id="angleInput" type="number" min="0" max="360" step="1" value="0">
      </div>
    </div>
    <div class="control">
      <label for="zoom">Zoom</label>
      <div class="row">
        <input id="zoom" type="range" min="0.6" max="1.6" step="0.05" value="1">
        <input id="zoomInput" type="number" min="0.6" max="1.6" step="0.05" value="1">
      </div>
    </div>
    <div class="control">
      <label for="order">Chain order</label>
      <input id="order" type="text" value="__DEFAULT_ORDER__">
      <div class="row" style="margin-top:8px;">
        <button id="applyOrder">Apply order</button>
        <button id="resetOrder" class="secondary">Reset</button>
      </div>
      <div class="hint" id="orderHint"></div>
    </div>
    <div class="control">
      <button id="downloadSvg">Download SVG</button>
      <div class="status" id="renderStatus">Render status: pending</div>
    </div>
  </div>
  <div class="plot-wrap panel" id="plotWrap" style="position:relative;">
    <svg id="circos" viewBox="0 0 900 900" role="img" aria-label="Circos plot"></svg>
    <div class="tooltip" id="tooltip"></div>
  </div>
</div>
<script>
const DATA = __DATA_JSON__;

const svg = document.getElementById('circos');
const tooltip = document.getElementById('tooltip');
const threshold = document.getElementById('threshold');
const thresholdInput = document.getElementById('thresholdInput');
const maxWidth = document.getElementById('maxWidth');
const maxWidthInput = document.getElementById('maxWidthInput');
const angle = document.getElementById('angle');
const angleInput = document.getElementById('angleInput');
const zoom = document.getElementById('zoom');
const zoomInput = document.getElementById('zoomInput');
const orderInput = document.getElementById('order');
const orderHint = document.getElementById('orderHint');
const renderStatus = document.getElementById('renderStatus');
const plotWrap = document.getElementById('plotWrap');

const chainMap = new Map(DATA.chains.map(c => [c.id, c]));
const defaultOrder = DATA.default_order.slice();

let filterChain = null;
let hoverArcPath = null;
let filterResidue = null;
let lastHover = null;

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

function parseOrder(raw) {
  const tokens = raw.split(/[\\s,]+/).filter(Boolean);
  const lowerMap = new Map(DATA.chains.map(c => [c.id.toLowerCase(), c.id]));
  const seen = new Set();
  const order = [];
  const invalid = [];
  for (const t of tokens) {
    const id = lowerMap.get(t.toLowerCase());
    if (id && !seen.has(id)) {
      order.push(id);
      seen.add(id);
    } else if (!id) {
      invalid.push(t);
    }
  }
  return { order, invalid };
}

function updateOrderHint(order, invalid) {
  const visible = new Set(order);
  const parts = DATA.chains.map(c => {
    const cls = visible.has(c.id) ? '' : 'hidden-chain';
    return `<span class="${cls}">${c.id}</span>`;
  });
  let msg = `Available chains: ${parts.join(', ')}`;
  if (invalid.length) {
    msg += ` <span class="error">Unknown: ${invalid.join(', ')}</span>`;
  }
  orderHint.innerHTML = msg;
}

function polar(cx, cy, r, angle) {
  return [cx + r * Math.cos(angle), cy + r * Math.sin(angle)];
}

function arcPath(cx, cy, rOuter, rInner, start, end) {
  const [x1, y1] = polar(cx, cy, rOuter, start);
  const [x2, y2] = polar(cx, cy, rOuter, end);
  const [x3, y3] = polar(cx, cy, rInner, end);
  const [x4, y4] = polar(cx, cy, rInner, start);
  const largeArc = end - start > Math.PI ? 1 : 0;
  return [
    `M ${x1} ${y1}`,
    `A ${rOuter} ${rOuter} 0 ${largeArc} 1 ${x2} ${y2}`,
    `L ${x3} ${y3}`,
    `A ${rInner} ${rInner} 0 ${largeArc} 0 ${x4} ${y4}`,
    'Z'
  ].join(' ');
}

function styledArcPath(cx, cy, outerR, innerR, start, end, startAtStart, flare, flareAngle) {
  const rCap = (outerR - innerR) / 2;
  const midR = (outerR + innerR) / 2;
  const capAngle = Math.atan(rCap / midR);
  if (startAtStart) {
    const endAdj = Math.max(start, end - capAngle);
    const flareEnd = Math.min(endAdj, start + flareAngle);
    const outerLargeArc = (endAdj - flareEnd) > Math.PI ? 1 : 0;
    const innerLargeArc = (endAdj - start) > Math.PI ? 1 : 0;
    const [sxo, syo] = polar(cx, cy, outerR + flare, start);
    const [fx, fy] = polar(cx, cy, outerR, flareEnd);
    const [exo, eyo] = polar(cx, cy, outerR, endAdj);
    const [exi, eyi] = polar(cx, cy, innerR, endAdj);
    const [sxi, syi] = polar(cx, cy, innerR, start);
    return [
      `M ${sxo} ${syo}`,
      `L ${fx} ${fy}`,
      `A ${outerR} ${outerR} 0 ${outerLargeArc} 1 ${exo} ${eyo}`,
      `A ${rCap} ${rCap} 0 0 1 ${exi} ${eyi}`,
      `A ${innerR} ${innerR} 0 ${innerLargeArc} 0 ${sxi} ${syi}`,
      'Z'
    ].join(' ');
  }
  const startAdj = Math.min(end, start + capAngle);
  const flareStart = Math.max(startAdj, end - flareAngle);
  const outerLargeArc = (flareStart - startAdj) > Math.PI ? 1 : 0;
  const innerLargeArc = (end - startAdj) > Math.PI ? 1 : 0;
  const [sxo, syo] = polar(cx, cy, outerR, startAdj);
  const [fx, fy] = polar(cx, cy, outerR, flareStart);
  const [exo, eyo] = polar(cx, cy, outerR + flare, end);
  const [exi, eyi] = polar(cx, cy, innerR, end);
  const [sxi, syi] = polar(cx, cy, innerR, startAdj);
  return [
    `M ${sxo} ${syo}`,
    `A ${outerR} ${outerR} 0 ${outerLargeArc} 1 ${fx} ${fy}`,
    `L ${exo} ${eyo}`,
    `L ${exi} ${eyi}`,
    `A ${innerR} ${innerR} 0 ${innerLargeArc} 0 ${sxi} ${syi}`,
    `A ${rCap} ${rCap} 0 0 1 ${sxo} ${syo}`,
    'Z'
  ].join(' ');
}

function mixColor(c1, c2) {
  const r1 = parseInt(c1.slice(1,3), 16);
  const g1 = parseInt(c1.slice(3,5), 16);
  const b1 = parseInt(c1.slice(5,7), 16);
  const r2 = parseInt(c2.slice(1,3), 16);
  const g2 = parseInt(c2.slice(3,5), 16);
  const b2 = parseInt(c2.slice(5,7), 16);
  const r = Math.round((r1 + r2) / 2);
  const g = Math.round((g1 + g2) / 2);
  const b = Math.round((b1 + b2) / 2);
  return `rgb(${r}, ${g}, ${b})`;
}

function svgPoint(event) {
  const rect = svg.getBoundingClientRect();
  const vb = svg.viewBox.baseVal;
  const x = (event.clientX - rect.left) * (vb.width / rect.width);
  const y = (event.clientY - rect.top) * (vb.height / rect.height);
  return { x, y };
}

function positionTooltip(event) {
  const rect = plotWrap.getBoundingClientRect();
  const fontSize = parseFloat(getComputedStyle(document.body).fontSize) || 16;
  const offset = fontSize * 10;
  let x = event.clientX - rect.left + offset;
  let y = event.clientY - rect.top + offset * 0.6;
  const tt = tooltip.getBoundingClientRect();
  const maxX = rect.width - tt.width - 8;
  const maxY = rect.height - tt.height - 8;
  if (x > maxX) x = event.clientX - rect.left - offset - tt.width;
  if (y > maxY) y = event.clientY - rect.top - offset - tt.height;
  x = clamp(x, 8, Math.max(8, maxX));
  y = clamp(y, 8, Math.max(8, maxY));
  tooltip.style.left = `${x}px`;
  tooltip.style.top = `${y}px`;
}

function residueAtEvent(chain, chainMeta, event, cx, cy) {
  const { x, y } = svgPoint(event);
  let ang = Math.atan2(y - cy, x - cx);
  if (ang < 0) ang += Math.PI * 2;
  let adj = ang;
  if (adj < chain.start) adj += Math.PI * 2;
  const frac = clamp((adj - chain.start) / (chain.end - chain.start), 0, 1);
  const idx = clamp(Math.floor(frac * chain.length) + 1, 1, chain.length);
  const info = chainMeta.resinfo[idx - 1];
  return { frac, idx, info };
}

function render() {
  const parsed = parseOrder(orderInput.value);
  let order = parsed.order;
  if (!order.length) {
    order = defaultOrder.slice();
    updateOrderHint(order, parsed.invalid.length ? parsed.invalid : ['(none)']);
  } else {
    updateOrderHint(order, parsed.invalid);
  }
  const maxCount = DATA.max_count;
  const gapDeg = 2;
  const gap = gapDeg * Math.PI / 180;
  const angleOffset = Number(angle.value) * Math.PI / 180;
  const thresholdValue = Number(threshold.value);
  const maxStroke = Number(maxWidth.value);
  const zoomValue = Number(zoom.value);
  const baseSize = Math.max(320, Math.floor(Math.min(plotWrap.clientWidth, plotWrap.clientHeight) - 20));
  const svgSize = Math.floor(baseSize * zoomValue);
  const minStroke = 0.2;
  const minAlpha = 0.02;
  const maxAlpha = 0.9;

  svg.setAttribute('width', svgSize.toString());
  svg.setAttribute('height', svgSize.toString());

  const totalLength = order.reduce((sum, id) => sum + chainMap.get(id).length, 0);
  const available = 2 * Math.PI - gap * order.length;
  const scale = available / totalLength;

  const chainAngles = new Map();
  let cursor = angleOffset;
  for (const id of order) {
    const chain = chainMap.get(id);
    const start = cursor;
    const end = start + chain.length * scale;
    chainAngles.set(id, { start, end, length: chain.length, color: chain.color });
    cursor = end + gap;
  }

  while (svg.firstChild) svg.removeChild(svg.firstChild);

  const cx = 450;
  const cy = 450;
  const outerR = 360;
  const innerR = 330;
  const linkR = 320;
  const controlR = 140;

  const arcsGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const linksGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const labelsGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const overlayGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  overlayGroup.setAttribute('pointer-events', 'none');

  for (const id of order) {
    const chain = chainAngles.get(id);
    const chainMeta = chainMap.get(id);
    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    if (chain.length >= 10) {
      const indicatorResidues = Math.max(2, Math.min(6, Math.round(chain.length * 0.03)));
      const segmentAngle = (indicatorResidues / chain.length) * (chain.end - chain.start);
      const startAtStart = (chainMeta.start_pos || 1) === 1;
      const flare = 6;
      path.setAttribute(
        'd',
        styledArcPath(cx, cy, outerR, innerR, chain.start, chain.end, startAtStart, flare, segmentAngle)
      );
    } else {
      path.setAttribute('d', arcPath(cx, cy, outerR, innerR, chain.start, chain.end));
    }
    path.setAttribute('fill', chain.color);
    path.setAttribute('opacity', '0.9');
    path.style.cursor = 'pointer';
    path.addEventListener('pointerdown', (event) => {
      if (event.shiftKey) {
        const res = (lastHover && lastHover.chain === id)
          ? lastHover
          : residueAtEvent(chain, chainMeta, event, cx, cy);
        const info = res.info;
        const label = info.pair
          ? `DNA bp ${info.resnum}`
          : `${info.source_chain}:${info.resname}${info.resnum}`;
        filterResidue = {
          displayChain: id,
          displayPos: res.idx,
          source_chain: info.source_chain,
          resnum: info.resnum,
          isDNA: !!info.pair,
          label
        };
        filterChain = null;
      } else {
        filterChain = id;
        filterResidue = null;
      }
      safeRender();
    });
    path.addEventListener('mousemove', (event) => {
      const { frac, idx, info } = residueAtEvent(chain, chainMeta, event, cx, cy);
      lastHover = { chain: id, idx, info };

      if (hoverArcPath) {
        const highlight = arcPath(cx, cy, innerR - 1, innerR - 7, chain.start, chain.start + frac * (chain.end - chain.start));
        hoverArcPath.setAttribute('d', highlight);
        hoverArcPath.setAttribute('fill', chain.color);
        hoverArcPath.setAttribute('opacity', '0.55');
      }

      const formatInfo = (info) => {
        if (info.pair) {
          return `DNA bp ${info.resnum}: ${info.pair.f_chain}:${info.pair.f_resname}${info.pair.f_resnum} / ${info.pair.r_chain}:${info.pair.r_resname}${info.pair.r_resnum}`;
        }
        return `${info.source_chain}:${info.resname}${info.resnum}`;
      };
      tooltip.innerHTML = `<div><strong>${formatInfo(info)}</strong></div>`;
      tooltip.style.opacity = '1';
      positionTooltip(event);
    });
    path.addEventListener('mouseleave', () => {
      if (hoverArcPath) hoverArcPath.setAttribute('d', '');
      lastHover = null;
      tooltip.style.opacity = '0';
    });
    arcsGroup.appendChild(path);

    const mid = (chain.start + chain.end) / 2;
    const [tx, ty] = polar(cx, cy, outerR + 18, mid);
    const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    label.setAttribute('x', tx);
    label.setAttribute('y', ty);
    label.setAttribute('fill', '#0f172a');
    label.setAttribute('font-size', '12');
    label.setAttribute('font-weight', '600');
    label.setAttribute('text-anchor', 'middle');
    label.textContent = id;
    labelsGroup.appendChild(label);
  }

  hoverArcPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
  overlayGroup.appendChild(hoverArcPath);

  for (const link of DATA.contacts) {
    const a = chainAngles.get(link.a);
    const b = chainAngles.get(link.b);
    if (!a || !b) continue;
    if (filterResidue) {
      const aInfo = chainMap.get(link.a).resinfo[link.a_pos - 1];
      const bInfo = chainMap.get(link.b).resinfo[link.b_pos - 1];
      let match = false;
      if (filterResidue.isDNA) {
        match =
          (link.a === filterResidue.displayChain && link.a_pos === filterResidue.displayPos) ||
          (link.b === filterResidue.displayChain && link.b_pos === filterResidue.displayPos);
      } else {
        match =
          (aInfo.source_chain === filterResidue.source_chain && aInfo.resnum === filterResidue.resnum) ||
          (bInfo.source_chain === filterResidue.source_chain && bInfo.resnum === filterResidue.resnum);
      }
      if (!match) continue;
    } else if (filterChain && link.a !== filterChain && link.b !== filterChain) {
      continue;
    }

    const aPos = link.a_pos - 0.5;
    const bPos = link.b_pos - 0.5;
    const aAngle = a.start + (aPos / a.length) * (a.end - a.start);
    const bAngle = b.start + (bPos / b.length) * (b.end - b.start);
    const [x1, y1] = polar(cx, cy, linkR, aAngle);
    const [x2, y2] = polar(cx, cy, linkR, bAngle);
    const [c1x, c1y] = polar(cx, cy, controlR, aAngle);
    const [c2x, c2y] = polar(cx, cy, controlR, bAngle);

    const count = link.count;
    const t = clamp((count - thresholdValue) / Math.max(1, maxCount - thresholdValue), 0, 1);
    const eased = Math.pow(t, 0.65);
    const alpha = minAlpha + (maxAlpha - minAlpha) * eased;
    const width = minStroke + (maxStroke - minStroke) * eased;
    const color = mixColor(a.color, b.color);

    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    path.setAttribute('d', `M ${x1} ${y1} C ${c1x} ${c1y} ${c2x} ${c2y} ${x2} ${y2}`);
    path.setAttribute('fill', 'none');
    path.setAttribute('stroke', color);
    path.setAttribute('stroke-opacity', alpha.toString());
    path.setAttribute('stroke-width', width.toString());
    path.setAttribute('stroke-linecap', 'round');

    path.addEventListener('mousemove', (event) => {
      const aInfo = chainMap.get(link.a).resinfo[link.a_pos - 1];
      const bInfo = chainMap.get(link.b).resinfo[link.b_pos - 1];
      const formatInfo = (info) => {
        if (info.pair) {
          return `DNA bp ${info.resnum}: ${info.pair.f_chain}:${info.pair.f_resname}${info.pair.f_resnum} / ${info.pair.r_chain}:${info.pair.r_resname}${info.pair.r_resnum}`;
        }
        return `${info.source_chain}:${info.resname}${info.resnum}`;
      };
      tooltip.innerHTML = `
        <div><strong>${formatInfo(aInfo)}</strong> ↔ <strong>${formatInfo(bInfo)}</strong></div>
        <div>Count: ${count}</div>
      `;
      tooltip.style.opacity = '1';
      positionTooltip(event);
    });
    path.addEventListener('mouseleave', () => {
      tooltip.style.opacity = '0';
    });

    linksGroup.appendChild(path);
  }

  svg.appendChild(linksGroup);
  svg.appendChild(arcsGroup);
  svg.appendChild(overlayGroup);
  svg.appendChild(labelsGroup);
  return {
    linkCount: linksGroup.childNodes.length,
    chainCount: order.length
  };
}

function safeRender() {
  try {
    const stats = render();
    let filterNote = '';
    if (filterResidue) {
      filterNote = ` Filtered: ${filterResidue.label}.`;
    } else if (filterChain) {
      filterNote = ` Filtered: ${filterChain}.`;
    }
    renderStatus.textContent = `Rendered ${stats.linkCount} links across ${stats.chainCount} chains.${filterNote}`;
  } catch (err) {
    renderStatus.textContent = `Render error: ${err.message}`;
    console.error(err);
  }
}

function syncRange(range, input) {
  range.addEventListener('input', () => {
    input.value = range.value;
    safeRender();
  });
  input.addEventListener('change', () => {
    range.value = clamp(Number(input.value), Number(range.min), Number(range.max));
    input.value = range.value;
    safeRender();
  });
}

syncRange(threshold, thresholdInput);
syncRange(maxWidth, maxWidthInput);
syncRange(angle, angleInput);
syncRange(zoom, zoomInput);

document.getElementById('applyOrder').addEventListener('click', () => safeRender());
document.getElementById('resetOrder').addEventListener('click', () => {
  orderInput.value = defaultOrder.join(', ');
  safeRender();
});

window.addEventListener('pointerup', () => {
  if (filterChain || filterResidue) {
    filterChain = null;
    filterResidue = null;
    safeRender();
  }
});

document.getElementById('downloadSvg').addEventListener('click', () => {
  const serializer = new XMLSerializer();
  const svgString = serializer.serializeToString(svg);
  const blob = new Blob([svgString], {type: 'image/svg+xml'});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'contacts_circos.svg';
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
});

safeRender();
</script>
</body>
</html>
"""
    html = (
        template.replace("__TITLE__", title)
        .replace("__MAX_COUNT__", str(max_count))
        .replace("__DEFAULT_ORDER__", ", ".join(default_order))
        .replace("__DATA_JSON__", data_json)
    )
    output_path.write_text(html)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create an interactive circos-style contact plot from ChimeraX contacts."
    )
    parser.add_argument(
        "--contacts-dir",
        type=Path,
        default=Path("data/chimerax/SPORM/spo-dsb-full-cifs/correct"),
        help="Directory containing .contacts and .cif files.",
    )
    parser.add_argument("--contacts-glob", default="*.contacts")
    parser.add_argument("--cif-glob", default="*.cif")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output HTML path (default: contacts_circos.html in contacts dir).",
    )
    parser.add_argument("--title", default="ChimeraX Contact Circos")
    parser.add_argument(
        "--dna-chains",
        nargs="*",
        help="Optional explicit DNA chain IDs (overrides auto-detect).",
    )
    parser.add_argument(
        "--dna-reverse",
        help=(
            "DNA chain ID to place in reverse orientation in the combined DNA chromosome "
            "(duplex mode)."
        ),
    )
    args = parser.parse_args()

    contacts_dir: Path = args.contacts_dir
    contact_paths = sorted(contacts_dir.glob(args.contacts_glob))
    cif_paths = sorted(contacts_dir.glob(args.cif_glob))
    if not contact_paths:
        raise SystemExit(f"No contact files found in {contacts_dir}")
    if not cif_paths:
        raise SystemExit(f"No CIF files found in {contacts_dir}")

    chain_res = parse_cifs(cif_paths)
    dna_chains = args.dna_chains or detect_dna_chains(chain_res)
    if dna_chains:
        sys.stderr.write(f"DNA chains: {', '.join(dna_chains)}\n")
    else:
        sys.stderr.write("No DNA chains detected; plotting all chains independently.\n")

    pos_map, pos_info, display_chain_of, display_chains, start_pos_map = build_chain_maps(
        chain_res, dna_chains, args.dna_reverse
    )

    raw_counts = parse_contacts(contact_paths)
    contacts, max_count = aggregate_contacts(raw_counts, pos_map, display_chain_of)

    chain_lengths = {
        chain: len(pos_info.get(chain, [])) for chain in display_chains
    }
    chain_colors = build_colors(display_chains)
    default_order = display_chains[:]

    output_path = args.output or contacts_dir / "contacts_circos.html"
    generate_html(
        args.title,
        display_chains,
        chain_lengths,
        chain_colors,
        pos_info,
        contacts,
        max_count,
        default_order,
        start_pos_map,
        output_path,
    )
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
