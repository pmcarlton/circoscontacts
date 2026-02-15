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
DNA_BASE_MAP = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "DI": "I", "DU": "U"}


def find_atom_site_category(doc: "gemmi.cif.Document") -> Dict[str, List[str]] | None:
    for block in doc:
        cat = block.get_mmcif_category("_atom_site")
        if cat:
            return cat
    return None


def parse_cifs(cif_paths: Iterable[Path]) -> Tuple[Dict[str, Dict[int, str]], Dict[str, List[int]]]:
    chain_res: Dict[str, Dict[int, str]] = defaultdict(dict)
    chain_resnums: Dict[str, List[int]] = defaultdict(list)
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
    for chain_id, residues in chain_res.items():
        chain_resnums[chain_id] = sorted(residues)
    if conflicts:
        sys.stderr.write(
            "Warning: residue name conflicts across CIFs (showing first 5):\n"
        )
        for line in conflicts[:5]:
            sys.stderr.write(f"  {line}\n")
    return chain_res, chain_resnums


def detect_dna_chains(chain_res: Dict[str, Dict[int, str]]) -> List[str]:
    dna_chains = []
    for chain_id, residues in chain_res.items():
        if not residues:
            continue
        dna_count = sum(1 for r in residues.values() if r in DNA_RESNAMES)
        if dna_count / max(1, len(residues)) >= 0.5:
            dna_chains.append(chain_id)
    return sorted(dna_chains)


def dna_sequence(residues: Dict[int, str], resnums: List[int]) -> str:
    seq = []
    for r in resnums:
        seq.append(DNA_BASE_MAP.get(residues.get(r, ""), "N"))
    return "".join(seq)


def reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "I": "I", "U": "A", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def count_mismatches(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def detect_split_dna(
    dna_chains: List[str],
    chain_res: Dict[str, Dict[int, str]],
    chain_resnums: Dict[str, List[int]],
    mismatches: int,
) -> Dict[str, object] | None:
    if len(dna_chains) != 3:
        return None
    lengths = {c: len(chain_resnums.get(c, [])) for c in dna_chains}
    long_chain = max(lengths, key=lengths.get)
    if list(lengths.values()).count(lengths[long_chain]) > 1:
        return None
    short_chains = [c for c in dna_chains if c != long_chain]
    if lengths[short_chains[0]] + lengths[short_chains[1]] != lengths[long_chain]:
        return None

    long_seq = dna_sequence(chain_res[long_chain], chain_resnums[long_chain])
    short_seqs = {
        c: dna_sequence(chain_res[c], chain_resnums[c]) for c in short_chains
    }

    order1 = short_chains
    order2 = short_chains[::-1]
    seq1 = short_seqs[order1[0]] + short_seqs[order1[1]]
    seq2 = short_seqs[order2[0]] + short_seqs[order2[1]]

    rc_long = reverse_complement(long_seq)
    mism1 = count_mismatches(seq1, rc_long)
    mism2 = count_mismatches(seq2, rc_long)

    best = None
    if mism1 <= mismatches:
        best = (order1, mism1, seq1)
    if mism2 <= mismatches and (best is None or mism2 < best[1]):
        best = (order2, mism2, seq2)
    if best is None:
        return None

    order, mism, comp_seq = best
    return {
        "long_chain": long_chain,
        "short_order": order,
        "mismatches": mism,
        "long_seq": long_seq,
        "comp_seq": comp_seq,
    }


def parse_contacts(contact_paths: Iterable[Path]) -> Tuple[Counter, Counter]:
    counts_atom: Counter = Counter()
    counts_residue: Counter = Counter()
    for path in contact_paths:
        per_file_pairs = set()
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
                counts_atom[key] += 1
                per_file_pairs.add(key)
        for key in per_file_pairs:
            counts_residue[key] += 1
    return counts_atom, counts_residue


def canonical_contact(chain1: str, res1: int, chain2: str, res2: int) -> Tuple[str, int, str, int]:
    if (chain1, res1) <= (chain2, res2):
        return chain1, res1, chain2, res2
    return chain2, res2, chain1, res1


def build_chain_maps(
    chain_res: Dict[str, Dict[int, str]],
    chain_resnums: Dict[str, List[int]],
    dna_chains: List[str],
    dna_reverse: str | None,
    dna_mismatches: int,
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
        resnums = chain_resnums.get(chain_id, sorted(residues))
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
        resnums = chain_resnums.get(dna_chain, sorted(residues))
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

    split = detect_split_dna(dna_chains, chain_res, chain_resnums, dna_mismatches)
    if split:
        long_chain = split["long_chain"]
        short_order = split["short_order"]
        long_resnums = chain_resnums.get(long_chain, sorted(chain_res[long_chain]))
        comp_entries = []
        for chain_id in short_order:
            resnums = chain_resnums.get(chain_id, sorted(chain_res[chain_id]))
            for resnum in resnums:
                comp_entries.append(
                    {
                        "chain": chain_id,
                        "resnum": resnum,
                        "resname": chain_res[chain_id][resnum],
                    }
                )
        length = min(len(long_resnums), len(comp_entries))
        pos_map[long_chain] = {}
        for idx, resnum in enumerate(long_resnums[:length]):
            pos_map[long_chain][resnum] = idx + 1
        for entry_idx, entry in enumerate(comp_entries[:length]):
            pos = length - entry_idx
            pos_map.setdefault(entry["chain"], {})[entry["resnum"]] = pos

        dna_info: List[Dict[str, str | int]] = []
        for i in range(length):
            long_resnum = long_resnums[i]
            comp_entry = comp_entries[length - i - 1]
            dna_info.append(
                {
                    "chain": "DNA",
                    "resnum": i + 1,
                    "resname": "bp",
                    "source_chain": "DNA",
                    "pair": {
                        "f_chain": long_chain,
                        "f_resnum": long_resnum,
                        "f_resname": chain_res[long_chain].get(long_resnum, "?"),
                        "r_chain": comp_entry["chain"],
                        "r_resnum": comp_entry["resnum"],
                        "r_resname": comp_entry["resname"],
                    },
                }
            )
        pos_info["DNA"] = dna_info
        dna_first = sorted(dna_chains)[0]
        first_resnum = 1
        if first_resnum not in chain_res.get(dna_first, {}):
            first_resnum = min(chain_res.get(dna_first, {1: ""}).keys())
        start_pos_map["DNA"] = pos_map[dna_first].get(first_resnum, 1)
        if dna_info:
            first = dna_info[0]["pair"]
            last = dna_info[-1]["pair"]
            sys.stderr.write(
                "Detected split DNA: long="
                f"{long_chain} shorts={','.join(short_order)} mismatches={split['mismatches']} "
                f"bp1={first['f_chain']}:{first['f_resnum']}/{first['r_chain']}:{first['r_resnum']} "
                f"bpN={last['f_chain']}:{last['f_resnum']}/{last['r_chain']}:{last['r_resnum']}\n"
            )
        else:
            sys.stderr.write(
                "Detected split DNA: long="
                f"{long_chain} shorts={','.join(short_order)} mismatches={split['mismatches']}\n"
            )
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
    f_resnums = chain_resnums.get(dna_forward, sorted(f_residues))
    r_resnums = chain_resnums.get(dna_rev, sorted(r_residues))

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
    first_resnum = 1
    if first_resnum not in chain_res.get(dna_first, {}):
        first_resnum = min(chain_res.get(dna_first, {1: ""}).keys())
    start_pos_map["DNA"] = pos_map[dna_first].get(first_resnum, 1)

    return pos_map, pos_info, display_chain_of, display_chains, start_pos_map


def aggregate_contacts(
    counts_atom: Counter,
    counts_residue: Counter,
    pos_map: Dict[str, Dict[int, int]],
    display_chain_of: Dict[str, str],
) -> Tuple[List[Dict[str, int | str | list]], int, int]:
    agg_atom: Counter = Counter()
    agg_residue: Counter = Counter()
    sources: Dict[Tuple[str, int, str, int], Dict[str, set]] = defaultdict(
        lambda: {"a": set(), "b": set()}
    )
    skipped = 0
    for chain1, res1, chain2, res2 in counts_atom:
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
        swapped = key != (a_chain, a_pos, b_chain, b_pos)
        if swapped:
            sources[key]["a"].add((chain2, res2))
            sources[key]["b"].add((chain1, res1))
        else:
            sources[key]["a"].add((chain1, res1))
            sources[key]["b"].add((chain2, res2))
        agg_atom[key] += counts_atom[(chain1, res1, chain2, res2)]

    for chain1, res1, chain2, res2 in counts_residue:
        if chain1 not in pos_map or chain2 not in pos_map:
            continue
        if res1 not in pos_map[chain1] or res2 not in pos_map[chain2]:
            continue
        a_chain = display_chain_of.get(chain1, chain1)
        b_chain = display_chain_of.get(chain2, chain2)
        a_pos = pos_map[chain1][res1]
        b_pos = pos_map[chain2][res2]
        key = canonical_contact(a_chain, a_pos, b_chain, b_pos)
        agg_residue[key] += counts_residue[(chain1, res1, chain2, res2)]
    if skipped:
        sys.stderr.write(f"Warning: skipped {skipped} contacts without CIF mapping\n")

    max_count_atom = max(agg_atom.values()) if agg_atom else 1
    max_count_residue = max(agg_residue.values()) if agg_residue else 1
    contacts = []
    for (a, pa, b, pb), c in agg_atom.items():
        src = sources.get((a, pa, b, pb), {"a": set(), "b": set()})
        contacts.append(
            {
                "a": a,
                "a_pos": pa,
                "b": b,
                "b_pos": pb,
                "count_atom": c,
                "count_residue": agg_residue.get((a, pa, b, pb), 0),
                "sources_a": [
                    {"chain": chain, "resnum": resnum}
                    for chain, resnum in sorted(src["a"])
                ],
                "sources_b": [
                    {"chain": chain, "resnum": resnum}
                    for chain, resnum in sorted(src["b"])
                ],
            }
        )
    contacts.sort(key=lambda x: x["count_atom"], reverse=True)
    return contacts, max_count_atom, max_count_residue


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
    max_count_atom: int,
    max_count_residue: int,
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
        "max_count_atom": max_count_atom,
        "max_count_residue": max_count_residue,
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
.hint .chain-toggle {
  display: flex;
  align-items: center;
  gap: 8px;
  margin: 4px 0;
}
.hint .flip-toggle {
  display: inline-flex;
  align-items: center;
  margin-left: 4px;
  font-size: 12px;
  cursor: pointer;
  user-select: none;
}
.hint .flip-toggle.active {
  font-weight: 600;
}
.hint .orig-label {
  min-width: 22px;
  font-weight: 600;
}
.hint .name-input {
  flex: 1;
  min-width: 90px;
  padding: 4px 6px;
  border: 1px solid #cbd5f5;
  border-radius: 6px;
  font-size: 12px;
}
.hint input[type="checkbox"] {
  transform: translateY(1px);
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
.tooltip.clickable {
  background: dodgerblue;
  box-shadow: 0 10px 24px rgba(30, 144, 255, 0.45);
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
      <label>Contact mode</label>
      <div class="row">
        <label class="chain-toggle"><input type="radio" name="countMode" value="atom" checked>Atom</label>
        <label class="chain-toggle"><input type="radio" name="countMode" value="residue">Residue</label>
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
      <label for="order">Chain Presence/Order</label>
      <input id="order" type="text" value="__DEFAULT_ORDER__">
      <div class="row" style="margin-top:8px;">
        <button id="applyOrder">Apply order</button>
        <button id="resetOrder" class="secondary">Reset</button>
      </div>
      <div class="hint" id="orderHint"></div>
      <label for="labelSize" style="margin-top:8px;">Label font size</label>
      <div class="row">
        <input id="labelSize" type="range" min="8" max="24" step="1" value="12">
        <input id="labelSizeInput" type="number" min="8" max="24" step="1" value="12">
      </div>
    </div>
    <div class="control">
      <button id="downloadSvg">Download SVG</button>
      <button id="exportCxc" class="secondary" style="margin-left:8px;">ChimeraX Colors</button>
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
const labelSize = document.getElementById('labelSize');
const labelSizeInput = document.getElementById('labelSizeInput');
const orderInput = document.getElementById('order');
const orderHint = document.getElementById('orderHint');
const renderStatus = document.getElementById('renderStatus');
const plotWrap = document.getElementById('plotWrap');
const exportCxc = document.getElementById('exportCxc');

const chainMap = new Map(DATA.chains.map(c => [c.id, c]));
const defaultOrder = DATA.default_order.slice();
let currentOrder = defaultOrder.slice();
const contactVisibility = new Map(DATA.chains.map(c => [c.id, true]));
const chainFlip = new Map(DATA.chains.map(c => [c.id, false]));
const chainDisplayName = new Map(DATA.chains.map(c => [c.id, c.id]));
let countMode = 'atom';

let filterChain = null;
let hoverArcPath = null;
let filterResidue = null;
let lastHover = null;
let clickableResidues = new Map();

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

function escapeHtml(value) {
  return value
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/\"/g, '&quot;');
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
    const inOrder = visible.has(c.id);
    const cls = inOrder ? '' : 'hidden-chain';
    const checked = inOrder && contactVisibility.get(c.id) !== false ? 'checked' : '';
    const disabled = inOrder ? '' : 'disabled';
    const flipped = chainFlip.get(c.id) ? 'active' : '';
    const symbol = chainFlip.get(c.id) ? '◀' : '▶';
    const flipDisabled = inOrder ? '' : 'aria-disabled="true"';
    const displayName = escapeHtml(chainDisplayName.get(c.id) || c.id);
    const flipColor = c.color || '#475569';
    return `<div class="chain-toggle ${cls}">
      <input type="checkbox" data-role="contacts" data-chain="${c.id}" ${checked} ${disabled}>
      <span class="orig-label">${c.id}</span>
      <span class="flip-toggle ${flipped}" data-role="flip" data-chain="${c.id}" style="color:${flipColor}" ${flipDisabled}>${symbol}</span>
      <input class="name-input" type="text" data-role="name" data-chain="${c.id}" value="${displayName}" ${disabled}>
    </div>`;
  });
  let msg = `Chains:${parts.join('')}`;
  if (invalid.length) {
    msg += ` <span class="error">Unknown: ${invalid.join(', ')}</span>`;
  }
  orderHint.innerHTML = msg;
  orderHint.querySelectorAll('input[type="checkbox"][data-role="contacts"]').forEach((box) => {
    box.addEventListener('change', (event) => {
      const target = event.currentTarget;
      const chainId = target.getAttribute('data-chain');
      const checked = target.checked;
      contactVisibility.set(chainId, checked);
      safeRender();
    });
  });
  orderHint.querySelectorAll('.flip-toggle[data-role="flip"]').forEach((node) => {
    node.addEventListener('click', (event) => {
      const target = event.currentTarget;
      if (target.getAttribute('aria-disabled') === 'true') return;
      const chainId = target.getAttribute('data-chain');
      const next = !(chainFlip.get(chainId) === true);
      chainFlip.set(chainId, next);
      safeRender();
    });
  });
  orderHint.querySelectorAll('input[data-role="name"]').forEach((node) => {
    const commitName = (target) => {
      const chainId = target.getAttribute('data-chain');
      const lastValid = chainDisplayName.get(chainId) || chainId;
      const next = target.value.trim();
      if (!next) {
        target.value = lastValid;
        return false;
      }
      if (next === lastValid) {
        return false;
      }
      chainDisplayName.set(chainId, next);
      return true;
    };

    node.addEventListener('keydown', (event) => {
      const target = event.currentTarget;
      if (event.key === 'Enter') {
        event.preventDefault();
        if (commitName(target)) {
          safeRender();
        } else {
          target.value = chainDisplayName.get(target.getAttribute('data-chain')) || target.getAttribute('data-chain');
        }
      } else if (event.key === 'Tab') {
        if (commitName(target)) {
          safeRender();
        }
      }
    });

    node.addEventListener('blur', (event) => {
      const target = event.currentTarget;
      const chainId = target.getAttribute('data-chain');
      target.value = chainDisplayName.get(chainId) || chainId;
    });
  });
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

function mixHexColor(c1, c2) {
  const r1 = parseInt(c1.slice(1,3), 16);
  const g1 = parseInt(c1.slice(3,5), 16);
  const b1 = parseInt(c1.slice(5,7), 16);
  const r2 = parseInt(c2.slice(1,3), 16);
  const g2 = parseInt(c2.slice(3,5), 16);
  const b2 = parseInt(c2.slice(5,7), 16);
  const r = Math.round((r1 + r2) / 2);
  const g = Math.round((g1 + g2) / 2);
  const b = Math.round((b1 + b2) / 2);
  return `#${r.toString(16).padStart(2, '0')}${g.toString(16).padStart(2, '0')}${b.toString(16).padStart(2, '0')}`;
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
  const rawFrac = clamp((adj - chain.start) / (chain.end - chain.start), 0, 1);
  let idx = clamp(Math.floor(rawFrac * chain.length) + 1, 1, chain.length);
  let frac = rawFrac;
  if (chain.flip) {
    idx = chain.length - idx + 1;
    frac = 1 - rawFrac;
  }
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
  currentOrder = order.slice();
  for (const chainId of order) {
    if (!contactVisibility.has(chainId)) {
      contactVisibility.set(chainId, true);
    }
  }
  const maxCount = countMode === 'atom' ? DATA.max_count_atom : DATA.max_count_residue;
  threshold.max = maxCount;
  thresholdInput.max = maxCount;
  if (Number(threshold.value) > maxCount) {
    threshold.value = maxCount;
    thresholdInput.value = maxCount;
  }
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
  const visibleChains = new Set(order);

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
    chainAngles.set(id, {
      start,
      end,
      length: chain.length,
      color: chain.color,
      flip: chainFlip.get(id) === true
    });
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

  const displayPos = (chainId, pos) => {
    const chain = chainAngles.get(chainId);
    if (!chain) return pos;
    return chain.flip ? (chain.length - pos + 1) : pos;
  };

  const residueFilterMatch = (link) => {
    if (!filterResidue) return true;
    const aInfo = chainMap.get(link.a).resinfo[link.a_pos - 1];
    const bInfo = chainMap.get(link.b).resinfo[link.b_pos - 1];
    if (filterResidue.isDNA) {
      return (
        (link.a === filterResidue.displayChain && displayPos(link.a, link.a_pos) === filterResidue.displayPos) ||
        (link.b === filterResidue.displayChain && displayPos(link.b, link.b_pos) === filterResidue.displayPos)
      );
    }
    return (
      (aInfo.source_chain === filterResidue.source_chain && aInfo.resnum === filterResidue.resnum) ||
      (bInfo.source_chain === filterResidue.source_chain && bInfo.resnum === filterResidue.resnum)
    );
  };

  const chainFilterMatch = (link) => {
    if (!filterChain) return true;
    return link.a === filterChain || link.b === filterChain;
  };

  const countFor = (link) => countMode === 'atom' ? link.count_atom : link.count_residue;

  const linkVisible = (link) => {
    if (countFor(link) < thresholdValue) return false;
    if (!visibleChains.has(link.a) || !visibleChains.has(link.b)) return false;
    if (contactVisibility.get(link.a) === false || contactVisibility.get(link.b) === false) return false;
    if (!residueFilterMatch(link)) return false;
    if (!chainFilterMatch(link)) return false;
    return true;
  };

  clickableResidues = new Map(order.map(id => [id, new Set()]));
  for (const link of DATA.contacts) {
    if (!linkVisible(link)) continue;
    if (clickableResidues.has(link.a)) clickableResidues.get(link.a).add(displayPos(link.a, link.a_pos));
    if (clickableResidues.has(link.b)) clickableResidues.get(link.b).add(displayPos(link.b, link.b_pos));
  }

  for (const id of order) {
    const chain = chainAngles.get(id);
    const chainMeta = chainMap.get(id);
    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    if (chain.length >= 10) {
      const indicatorResidues = Math.max(2, Math.min(6, Math.round(chain.length * 0.03)));
      const segmentAngle = (indicatorResidues / chain.length) * (chain.end - chain.start);
      const startAtStart = ((chainMeta.start_pos || 1) === 1) !== (chain.flip === true);
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
      const clickable = clickableResidues.get(id) && clickableResidues.get(id).has(idx);
      tooltip.classList.toggle('clickable', Boolean(clickable));
      positionTooltip(event);
    });
    path.addEventListener('mouseleave', () => {
      if (hoverArcPath) hoverArcPath.setAttribute('d', '');
      lastHover = null;
      tooltip.classList.remove('clickable');
      tooltip.style.opacity = '0';
    });
    arcsGroup.appendChild(path);

    const mid = (chain.start + chain.end) / 2;
    const [tx, ty] = polar(cx, cy, outerR + 18, mid);
    const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    label.setAttribute('x', tx);
    label.setAttribute('y', ty);
    label.setAttribute('fill', '#0f172a');
    label.setAttribute('font-size', labelSize.value);
    label.setAttribute('font-weight', '600');
    label.setAttribute('text-anchor', 'middle');
    label.textContent = chainDisplayName.get(id) || id;
    labelsGroup.appendChild(label);
  }

  hoverArcPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
  overlayGroup.appendChild(hoverArcPath);

  for (const link of DATA.contacts) {
    const a = chainAngles.get(link.a);
    const b = chainAngles.get(link.b);
    if (!a || !b) continue;
    if (!linkVisible(link)) continue;

    const aPos = displayPos(link.a, link.a_pos) - 0.5;
    const bPos = displayPos(link.b, link.b_pos) - 0.5;
    const aAngle = a.start + (aPos / a.length) * (a.end - a.start);
    const bAngle = b.start + (bPos / b.length) * (b.end - b.start);
    const [x1, y1] = polar(cx, cy, linkR, aAngle);
    const [x2, y2] = polar(cx, cy, linkR, bAngle);
    const [c1x, c1y] = polar(cx, cy, controlR, aAngle);
    const [c2x, c2y] = polar(cx, cy, controlR, bAngle);

    const count = countFor(link);
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
        <div>Count (${countMode}): ${count}</div>
      `;
      tooltip.style.opacity = '1';
      tooltip.classList.remove('clickable');
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
syncRange(labelSize, labelSizeInput);

document.getElementById('applyOrder').addEventListener('click', () => safeRender());
document.getElementById('resetOrder').addEventListener('click', () => {
  orderInput.value = defaultOrder.join(', ');
  safeRender();
});

document.querySelectorAll('input[name="countMode"]').forEach((radio) => {
  radio.addEventListener('change', (event) => {
    countMode = event.currentTarget.value;
    safeRender();
  });
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

exportCxc.addEventListener('click', () => {
  const thresholdValue = Number(threshold.value);
  const visible = new Set(currentOrder);
  const lines = [];
  lines.push(`# ChimeraX colors from contacts_circos`);
  lines.push(`# threshold: ${thresholdValue}`);
  lines.push(`# chains: ${currentOrder.join(', ')}`);
  lines.push('color light gray');
  lines.push('color light gray target s');
  lines.push('select clear');

  for (const link of DATA.contacts) {
    const count = countMode === 'atom' ? link.count_atom : link.count_residue;
    if (count < thresholdValue) continue;
    if (!visible.has(link.a) || !visible.has(link.b)) continue;
    if (contactVisibility.get(link.a) === false || contactVisibility.get(link.b) === false) continue;
    const selectors = new Set();

    if (link.sources_a && link.sources_b) {
      for (const src of link.sources_a) {
        selectors.add(`/${src.chain}:${src.resnum}`);
      }
      for (const src of link.sources_b) {
        selectors.add(`/${src.chain}:${src.resnum}`);
      }
    } else {
      const addSelectors = (chainId, pos) => {
        const info = chainMap.get(chainId).resinfo[pos - 1];
        if (info && info.pair) {
          selectors.add(`/${info.pair.f_chain}:${info.pair.f_resnum}`);
          selectors.add(`/${info.pair.r_chain}:${info.pair.r_resnum}`);
        } else if (info) {
          selectors.add(`/${info.source_chain}:${info.resnum}`);
        }
      };

      addSelectors(link.a, link.a_pos);
      addSelectors(link.b, link.b_pos);
    }

    if (selectors.size === 0) continue;
    const color = mixHexColor(chainMap.get(link.a).color, chainMap.get(link.b).color);
    lines.push(`select ${Array.from(selectors).join(' ')}`);
    lines.push(`color sel ${color} target rs`);
    lines.push('select clear');
  }

  const blob = new Blob([lines.join('\\n') + '\\n'], {type: 'text/plain'});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'contacts_circos_colors.cxc';
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
        .replace("__MAX_COUNT__", str(max_count_atom))
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
        default=Path("."),
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
    parser.add_argument(
        "--dna-mismatches",
        type=int,
        default=0,
        help="Allowed mismatches when detecting split DNA (default: 0).",
    )
    args = parser.parse_args()

    contacts_dir: Path = args.contacts_dir
    contact_paths = sorted(contacts_dir.glob(args.contacts_glob))
    cif_paths = sorted(contacts_dir.glob(args.cif_glob))
    if not contact_paths:
        raise SystemExit(f"No contact files found in {contacts_dir}")
    if not cif_paths:
        raise SystemExit(f"No CIF files found in {contacts_dir}")

    chain_res, chain_resnums = parse_cifs(cif_paths)
    dna_chains = args.dna_chains or detect_dna_chains(chain_res)
    if dna_chains:
        sys.stderr.write(f"DNA chains: {', '.join(dna_chains)}\n")
    else:
        sys.stderr.write("No DNA chains detected; plotting all chains independently.\n")

    pos_map, pos_info, display_chain_of, display_chains, start_pos_map = build_chain_maps(
        chain_res, chain_resnums, dna_chains, args.dna_reverse, args.dna_mismatches
    )

    raw_counts_atom, raw_counts_residue = parse_contacts(contact_paths)
    contacts, max_count_atom, max_count_residue = aggregate_contacts(
        raw_counts_atom, raw_counts_residue, pos_map, display_chain_of
    )

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
        max_count_atom,
        max_count_residue,
        default_order,
        start_pos_map,
        output_path,
    )
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
