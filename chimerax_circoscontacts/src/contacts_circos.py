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
except Exception:  # pragma: no cover - optional for ChimeraX plugin usage
    gemmi = None


DNA_RESNAMES = {"DA", "DC", "DG", "DT", "DI", "DU"}
DNA_BASE_MAP = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "DI": "I", "DU": "U"}


def find_atom_site_category(doc: "gemmi.cif.Document") -> Dict[str, List[str]] | None:
    for block in doc:
        cat = block.get_mmcif_category("_atom_site")
        if cat:
            return cat
    return None


def parse_cifs(cif_paths: Iterable[Path]) -> Tuple[Dict[str, Dict[int, str]], Dict[str, List[int]]]:
    if gemmi is None:
        raise SystemExit("Missing dependency: gemmi. Install with `pip install gemmi`.")
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
    def extract_chain(token: str) -> str:
        if token.startswith("/"):
            return token[1:]
        if "/" in token:
            return token.rsplit("/", 1)[-1]
        return token

    counts_atom: Counter = Counter()
    counts_residue: Counter = Counter()
    for path in contact_paths:
        per_file_pairs = set()
        with path.open() as fh:
            for line in fh:
                parts = line.split()
                chain1 = chain2 = None
                res1 = res2 = None

                # Legacy ChimeraX contacts format:
                # /E ALA 204 O     /J DT 28 OP1      0.219    2.761
                if len(parts) >= 8 and parts[0].startswith("/") and parts[4].startswith("/"):
                    chain1 = extract_chain(parts[0])
                    chain2 = extract_chain(parts[4])
                    try:
                        res1 = int(parts[2])
                        res2 = int(parts[6])
                    except ValueError:
                        continue
                # Newer ChimeraX contacts format:
                # <file> #1/T GLN 305 OE1 <file> #1/S TYR 119 O 0.313 2.647
                elif len(parts) >= 10 and "/" in parts[1] and "/" in parts[6]:
                    chain1 = extract_chain(parts[1])
                    chain2 = extract_chain(parts[6])
                    try:
                        res1 = int(parts[3])
                        res2 = int(parts[8])
                    except ValueError:
                        continue
                else:
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
    dna_mode: str = "auto",
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
    if dna_mode == "split":
        all_chains = sorted(chain_res)
        for chain_id in all_chains:
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
            display_chains.append(chain_id)
        return pos_map, pos_info, display_chain_of, display_chains, start_pos_map

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
    active_ranges: Dict[str, List[Tuple[int, int]]],
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
                "active_ranges": active_ranges.get(chain, []),
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
  cursor: pointer;
  user-select: none;
}
.hint .orig-label.locked {
  color: #0b4ea2;
  text-decoration: underline;
  text-underline-offset: 2px;
}
.hint .orig-label[aria-disabled="true"] {
  cursor: default;
  opacity: 0.6;
  text-decoration: none;
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
.selection-menu {
  position: absolute;
  background: rgba(255, 255, 255, 0.97);
  border: 1px solid #cbd5e1;
  border-radius: 10px;
  box-shadow: 0 12px 30px rgba(15, 23, 42, 0.18);
  padding: 8px;
  min-width: 230px;
  z-index: 10;
  display: none;
}
.selection-menu .menu-title {
  font-size: 12px;
  color: #334155;
  margin: 0 0 6px 0;
  font-weight: 600;
}
.selection-menu .menu-block {
  border-top: 1px solid #e2e8f0;
  padding-top: 6px;
  margin-top: 6px;
}
.selection-menu .menu-row {
  display: flex;
  gap: 6px;
  align-items: center;
  margin-top: 4px;
}
.selection-menu button {
  font-size: 12px;
  padding: 5px 8px;
}
.selection-menu button:disabled {
  background: #cbd5e1;
  color: #64748b;
  cursor: not-allowed;
}
.selection-menu input[type="color"] {
  width: 34px;
  height: 24px;
  border: none;
  background: transparent;
  padding: 0;
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
      <label for="plotTitle">Plot title</label>
      <input id="plotTitle" type="text" value="__TITLE__">
    </div>
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
      <button id="saveHtml" class="secondary" style="margin-left:8px;">Save HTML</button>
      <button id="exportCxc" class="secondary" style="margin-left:8px;">ChimeraX Colors</button>
      <div class="row" style="margin-top:8px;">
        <button id="saveSession" class="secondary">Save Session</button>
        <button id="loadSession" class="secondary">Load Session</button>
        <input id="loadSessionFile" type="file" accept=".json,application/json" style="display:none;">
      </div>
      <div class="status" id="renderStatus">Render status: pending</div>
    </div>
  </div>
  <div class="plot-wrap panel" id="plotWrap" style="position:relative;">
    <svg id="circos" viewBox="0 0 1400 900" role="img" aria-label="Circos plot"></svg>
    <div class="tooltip" id="tooltip"></div>
    <div class="selection-menu" id="selectionMenu"></div>
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
const labelSize = document.getElementById('labelSize');
const labelSizeInput = document.getElementById('labelSizeInput');
const orderInput = document.getElementById('order');
const orderHint = document.getElementById('orderHint');
const renderStatus = document.getElementById('renderStatus');
const plotWrap = document.getElementById('plotWrap');
const saveHtml = document.getElementById('saveHtml');
const exportCxc = document.getElementById('exportCxc');
const selectionMenu = document.getElementById('selectionMenu');
const saveSessionBtn = document.getElementById('saveSession');
const loadSessionBtn = document.getElementById('loadSession');
const loadSessionFile = document.getElementById('loadSessionFile');
const plotTitleInput = document.getElementById('plotTitle');

const chainMap = new Map(DATA.chains.map(c => [c.id, c]));
const defaultOrder = DATA.default_order.slice();
let currentOrder = defaultOrder.slice();
const contactVisibility = new Map(DATA.chains.map(c => [c.id, true]));
const chainFlip = new Map(DATA.chains.map(c => [c.id, false]));
const chainDisplayName = new Map(DATA.chains.map(c => [c.id, c.id]));
let lockedBottomChain = null;
let countMode = 'atom';

let filterChain = null;
let hoverArcPath = null;
let filterResidue = null;
let lastHover = null;
let clickableResidues = new Map();
let arcSelections = [];
let nextSelectionId = 1;
let dragSelection = null;
let suppressNextClick = false;
let activeCalloutDrag = null;
const dnaSeqSideCache = new Map();
const CANVAS_W = 1400;
const CANVAS_H = 900;

const AA1 = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G',
  HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
  THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V', SEC: 'U', PYL: 'O'
};

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
    const labelDisabled = inOrder ? '' : 'aria-disabled="true"';
    const locked = lockedBottomChain === c.id ? 'locked' : '';
    const flipped = chainFlip.get(c.id) ? 'active' : '';
    const symbol = chainFlip.get(c.id) ? '◀' : '▶';
    const flipDisabled = inOrder ? '' : 'aria-disabled="true"';
    const displayName = escapeHtml(chainDisplayName.get(c.id) || c.id);
    const flipColor = c.color || '#475569';
    return `<div class="chain-toggle ${cls}">
      <input type="checkbox" data-role="contacts" data-chain="${c.id}" ${checked} ${disabled}>
      <span class="orig-label ${locked}" data-role="lock" data-chain="${c.id}" ${labelDisabled} title="Lock this chain midpoint at the bottom">${c.id}</span>
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
  orderHint.querySelectorAll('.orig-label[data-role="lock"]').forEach((node) => {
    node.addEventListener('click', (event) => {
      const target = event.currentTarget;
      if (target.getAttribute('aria-disabled') === 'true') return;
      const chainId = target.getAttribute('data-chain');
      if (!chainId || !chainMap.has(chainId)) return;
      lockedBottomChain = chainId;
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

function arcLinePath(cx, cy, r, start, end, reverse = false) {
  const s = reverse ? end : start;
  const e = reverse ? start : end;
  const [x1, y1] = polar(cx, cy, r, s);
  const [x2, y2] = polar(cx, cy, r, e);
  const largeArc = Math.abs(end - start) > Math.PI ? 1 : 0;
  const sweep = reverse ? 0 : 1;
  return `M ${x1} ${y1} A ${r} ${r} 0 ${largeArc} ${sweep} ${x2} ${y2}`;
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

function lightenHexColor(hex, amount = 0.45) {
  const c = String(hex || '#888888');
  const r = parseInt(c.slice(1, 3), 16);
  const g = parseInt(c.slice(3, 5), 16);
  const b = parseInt(c.slice(5, 7), 16);
  const nr = Math.round(r + (255 - r) * amount);
  const ng = Math.round(g + (255 - g) * amount);
  const nb = Math.round(b + (255 - b) * amount);
  return `#${nr.toString(16).padStart(2, '0')}${ng.toString(16).padStart(2, '0')}${nb.toString(16).padStart(2, '0')}`;
}

function svgPoint(event) {
  if (typeof svg.createSVGPoint === 'function' && svg.getScreenCTM) {
    const pt = svg.createSVGPoint();
    pt.x = event.clientX;
    pt.y = event.clientY;
    const ctm = svg.getScreenCTM();
    if (ctm) {
      const p = pt.matrixTransform(ctm.inverse());
      return { x: p.x, y: p.y };
    }
  }
  const rect = svg.getBoundingClientRect();
  const vb = svg.viewBox.baseVal;
  return {
    x: (event.clientX - rect.left) * (vb.width / rect.width),
    y: (event.clientY - rect.top) * (vb.height / rect.height)
  };
}

function positionTooltip(event) {
  const rect = plotWrap.getBoundingClientRect();
  const tt = tooltip.getBoundingClientRect();
  let x = (rect.width - tt.width) / 2;
  let y = (rect.height - tt.height) / 2;
  const maxX = rect.width - tt.width - 8;
  const maxY = rect.height - tt.height - 8;
  x = clamp(x, 8, Math.max(8, maxX));
  y = clamp(y, 8, Math.max(8, maxY));
  tooltip.style.left = `${x}px`;
  tooltip.style.top = `${y}px`;
}

function unwrapAngleForSpan(angle, start, end) {
  const twoPi = Math.PI * 2;
  const candidates = [angle - twoPi, angle, angle + twoPi, angle + twoPi * 2];
  for (const candidate of candidates) {
    if (candidate >= start && candidate <= end) return candidate;
  }
  const center = (start + end) / 2;
  let best = candidates[0];
  let bestDist = Math.abs(candidates[0] - center);
  for (const candidate of candidates.slice(1)) {
    const dist = Math.abs(candidate - center);
    if (dist < bestDist) {
      best = candidate;
      bestDist = dist;
    }
  }
  return best;
}

function residueAtEvent(chain, chainMeta, event, cx, cy) {
  const { x, y } = svgPoint(event);
  let ang = Math.atan2(y - cy, x - cx);
  if (ang < 0) ang += Math.PI * 2;
  const adj = unwrapAngleForSpan(ang, chain.start, chain.end);
  const rawFrac = clamp((adj - chain.start) / (chain.end - chain.start), 0, 1);
  let idx = clamp(Math.floor(rawFrac * chain.length) + 1, 1, chain.length);
  let frac = rawFrac;
  if (chain.flip) {
    idx = chain.length - idx + 1;
    frac = 1 - rawFrac;
  }
  const info = chainMeta.resinfo[idx - 1];
  return { rawFrac, frac, idx, info };
}

function normalizeRange(startPos, endPos) {
  return [Math.min(startPos, endPos), Math.max(startPos, endPos)];
}

function linkMatchesSelection(link, sel) {
  const [s, e] = normalizeRange(sel.startPos, sel.endPos);
  return (
    (link.a === sel.chainId && link.a_pos >= s && link.a_pos <= e) ||
    (link.b === sel.chainId && link.b_pos >= s && link.b_pos <= e)
  );
}

function selectionLength(sel) {
  const [s, e] = normalizeRange(sel.startPos, sel.endPos);
  return e - s + 1;
}

function baseCode(resname) {
  if (!resname) return 'N';
  const up = String(resname).toUpperCase();
  if (up === 'DT') return 'T';
  if (up === 'DA') return 'A';
  if (up === 'DC') return 'C';
  if (up === 'DG') return 'G';
  if (up === 'DU') return 'U';
  return up.length ? up[up.length - 1] : 'N';
}

function dnaPreferredSide(chainId, chainMeta) {
  if (dnaSeqSideCache.has(chainId)) return dnaSeqSideCache.get(chainId);
  const pairs = (chainMeta.resinfo || []).filter((x) => x && x.pair);
  if (!pairs.length) return null;
  const fChains = new Set(pairs.map((x) => x.pair.f_chain));
  const rChains = new Set(pairs.map((x) => x.pair.r_chain));
  let side = 'f';
  if (!(fChains.size === 1 && rChains.size > 1)) {
    const sample = pairs[0].pair;
    side = String(sample.f_chain).localeCompare(String(sample.r_chain)) <= 0 ? 'f' : 'r';
  }
  dnaSeqSideCache.set(chainId, side);
  return side;
}

function selectionSequenceText(sel) {
  const chainMeta = chainMap.get(sel.chainId);
  if (!chainMeta) return '';
  const dnaSide = dnaPreferredSide(sel.chainId, chainMeta);
  const [s, e] = normalizeRange(sel.startPos, sel.endPos);
  const tokens = [];
  for (let pos = s; pos <= e; pos += 1) {
    const info = chainMeta.resinfo[pos - 1];
    if (!info) continue;
    if (info.pair) {
      const one = dnaSide === 'r'
        ? baseCode(info.pair.r_resname)
        : baseCode(info.pair.f_resname);
      tokens.push(one);
    } else {
      const aa = AA1[String(info.resname || '').toUpperCase()] || 'X';
      tokens.push(aa);
    }
  }
  return tokens.join('');
}

function staticLinkVisible(link, thresholdValue, visibleChains) {
  const count = countMode === 'atom' ? link.count_atom : link.count_residue;
  if (count < thresholdValue) return false;
  if (!visibleChains.has(link.a) || !visibleChains.has(link.b)) return false;
  if (contactVisibility.get(link.a) === false || contactVisibility.get(link.b) === false) return false;
  return true;
}

function trimSelectionToContacts(sel, thresholdValue, visibleChains) {
  const [s, e] = normalizeRange(sel.startPos, sel.endPos);
  const touching = [];
  for (const link of DATA.contacts) {
    if (!staticLinkVisible(link, thresholdValue, visibleChains)) continue;
    if (link.a === sel.chainId && link.a_pos >= s && link.a_pos <= e) touching.push(link.a_pos);
    if (link.b === sel.chainId && link.b_pos >= s && link.b_pos <= e) touching.push(link.b_pos);
  }
  if (!touching.length) return null;
  return [Math.min(...touching), Math.max(...touching)];
}

function selectionDisplayBounds(sel, chain) {
  const a = chain.flip ? (chain.length - sel.startPos + 1) : sel.startPos;
  const b = chain.flip ? (chain.length - sel.endPos + 1) : sel.endPos;
  return normalizeRange(a, b);
}

function closeSelectionMenu() {
  selectionMenu.style.display = 'none';
  selectionMenu.innerHTML = '';
}

function wrapLines(text, maxWidthPx) {
  const charPx = 7.2;
  const maxChars = Math.max(1, Math.floor(maxWidthPx / charPx));
  const rawLines = String(text || '').split('\\n');
  const out = [];
  for (const raw of rawLines) {
    const words = raw.split(/\\s+/).filter(Boolean);
    if (!words.length) {
      out.push('');
      continue;
    }
    let line = '';
    for (const rawWord of words) {
      let word = rawWord;
      while (word.length > maxChars) {
        if (line.length) {
          out.push(line);
          line = '';
        }
        out.push(word.slice(0, maxChars));
        word = word.slice(maxChars);
      }
      if (!line.length) {
        line = word;
        continue;
      }
      const next = `${line} ${word}`;
      if (next.length <= maxChars) line = next;
      else {
        out.push(line);
        line = word;
      }
    }
    if (line.length || !words.length) out.push(line);
  }
  return out.length ? out : [''];
}

function nearestPointOnRect(x, y, rx, ry, rw, rh) {
  const nx = Math.max(rx, Math.min(x, rx + rw));
  const ny = Math.max(ry, Math.min(y, ry + rh));
  return [nx, ny];
}

function sanitizeEditableHtml(html) {
  const holder = document.createElement('div');
  holder.innerHTML = String(html || '');
  holder.querySelectorAll('script,style,iframe,object,embed').forEach((n) => n.remove());
  return holder.innerHTML;
}

function calloutRangeTitle(sel) {
  const [a, b] = normalizeRange(sel.startPos, sel.endPos);
  return `${a}\u2013${b}`;
}

function normalizeRanges(ranges, maxLen) {
  const cleaned = [];
  for (const range of ranges || []) {
    if (!Array.isArray(range) || range.length < 2) continue;
    let s = Math.max(1, Math.min(maxLen, Number(range[0])));
    let e = Math.max(1, Math.min(maxLen, Number(range[1])));
    if (Number.isNaN(s) || Number.isNaN(e)) continue;
    if (e < s) [s, e] = [e, s];
    cleaned.push([s, e]);
  }
  cleaned.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  if (!cleaned.length) return cleaned;
  const merged = [cleaned[0]];
  for (let i = 1; i < cleaned.length; i += 1) {
    const [s, e] = cleaned[i];
    const last = merged[merged.length - 1];
    if (s <= last[1] + 1) {
      last[1] = Math.max(last[1], e);
    } else {
      merged.push([s, e]);
    }
  }
  return merged;
}

function displayRangesForChain(chain, chainMeta) {
  const ranges = normalizeRanges(chainMeta.active_ranges || [], chain.length);
  if (!chain.flip) return ranges;
  const flipped = ranges.map(([s, e]) => [chain.length - e + 1, chain.length - s + 1]);
  return normalizeRanges(flipped, chain.length);
}

function render() {
  const parsed = parseOrder(orderInput.value);
  let order = parsed.order;
  const usingDefaultOrder = !order.length;
  if (usingDefaultOrder) {
    order = defaultOrder.slice();
  }
  if (lockedBottomChain && !order.includes(lockedBottomChain)) lockedBottomChain = null;
  if (usingDefaultOrder) updateOrderHint(order, parsed.invalid.length ? parsed.invalid : ['(none)']);
  else updateOrderHint(order, parsed.invalid);
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
  const thresholdValue = Number(threshold.value);
  const maxStroke = Number(maxWidth.value);
  const wrapRect = plotWrap.getBoundingClientRect();
  const availW = Math.max(320, Math.floor(wrapRect.width - 20));
  const availH = Math.max(320, Math.floor(wrapRect.height - 20));
  const fitScale = Math.min(availW / CANVAS_W, availH / CANVAS_H);
  // Keep callout text legible (~8pt+) while still auto-fitting the plot.
  const renderScale = Math.max(0.9, fitScale * 1.2);
  const svgW = Math.floor(CANVAS_W * renderScale);
  const svgH = Math.floor(CANVAS_H * renderScale);
  const minStroke = 0.2;
  const minAlpha = 0.02;
  const maxAlpha = 0.9;
  const visibleChains = new Set(order);

  svg.setAttribute('width', svgW.toString());
  svg.setAttribute('height', svgH.toString());

  const totalLength = order.reduce((sum, id) => sum + chainMap.get(id).length, 0);
  const available = 2 * Math.PI - gap * order.length;
  const scale = available / totalLength;
  let angleOffset = Number(angle.value) * Math.PI / 180;
  if (lockedBottomChain && order.includes(lockedBottomChain)) {
    let prefix = 0;
    for (const id of order) {
      if (id === lockedBottomChain) break;
      prefix += chainMap.get(id).length * scale + gap;
    }
    const halfSpan = chainMap.get(lockedBottomChain).length * scale * 0.5;
    angleOffset = (Math.PI / 2) - (prefix + halfSpan);
  }

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

  const cx = CANVAS_W / 2;
  const cy = CANVAS_H / 2;
  const outerR = 360;
  const innerR = 330;
  const linkR = 320;
  const controlR = 140;

  const arcsGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const selectionFillGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const selectionStrokeGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const calloutGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const linksGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const labelsGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const labelDefs = document.createElementNS('http://www.w3.org/2000/svg', 'defs');
  const overlayGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const targetGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  labelsGroup.setAttribute('pointer-events', 'none');
  overlayGroup.setAttribute('pointer-events', 'none');

  const plotTitle = (plotTitleInput.value || '').trim();
  if (plotTitle) {
    const titleNode = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    titleNode.setAttribute('x', `${cx}`);
    titleNode.setAttribute('y', '38');
    titleNode.setAttribute('text-anchor', 'middle');
    titleNode.setAttribute('fill', '#0f172a');
    titleNode.setAttribute('font-size', '24');
    titleNode.setAttribute('font-weight', '700');
    titleNode.setAttribute('font-family', 'IBM Plex Sans, sans-serif');
    titleNode.textContent = plotTitle;
    svg.appendChild(titleNode);
  }

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
        (link.a === filterResidue.displayChain && link.a_pos === filterResidue.dataPos) ||
        (link.b === filterResidue.displayChain && link.b_pos === filterResidue.dataPos)
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
    if (clickableResidues.has(link.a)) clickableResidues.get(link.a).add(link.a_pos);
    if (clickableResidues.has(link.b)) clickableResidues.get(link.b).add(link.b_pos);
  }

  const anyActiveMask = order.some((chainId) => {
    const meta = chainMap.get(chainId);
    const chain = chainAngles.get(chainId);
    return displayRangesForChain(chain, meta).length > 0;
  });

  for (const id of order) {
    const chain = chainAngles.get(id);
    const chainMeta = chainMap.get(id);
    const activeRanges = displayRangesForChain(chain, chainMeta);
    const hasActiveMask = activeRanges.length > 0;
    const fullyActive =
      hasActiveMask &&
      activeRanges.length === 1 &&
      activeRanges[0][0] === 1 &&
      activeRanges[0][1] === chain.length;

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
    const baseOpacity = (hasActiveMask && !fullyActive) || (anyActiveMask && !hasActiveMask) ? '0.45' : '0.9';
    path.setAttribute('opacity', baseOpacity);
    path.style.cursor = 'pointer';
    path.addEventListener('pointerdown', (event) => {
      closeSelectionMenu();
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
          dataPos: res.idx,
          source_chain: info.source_chain,
          resnum: info.resnum,
          isDNA: !!info.pair,
          label
        };
        filterChain = null;
      } else {
        const res = residueAtEvent(chain, chainMeta, event, cx, cy);
        dragSelection = {
          chainId: id,
          startPos: res.idx,
          currentPos: res.idx,
          moved: false
        };
        filterChain = id;
        filterResidue = null;
      }
      safeRender();
    });
    path.addEventListener('mousemove', (event) => {
      const { rawFrac, idx, info } = residueAtEvent(chain, chainMeta, event, cx, cy);
      lastHover = { chain: id, idx, info };
      if (dragSelection && dragSelection.chainId === id) {
        if (dragSelection.currentPos !== idx) dragSelection.moved = true;
        dragSelection.currentPos = idx;
        safeRender();
        return;
      }

      if (hoverArcPath) {
        const span = chain.end - chain.start;
        const rawAngle = chain.start + rawFrac * span;
        const hStart = chain.flip ? rawAngle : chain.start;
        const hEnd = chain.flip ? chain.end : rawAngle;
        const highlight = hEnd > hStart
          ? arcPath(cx, cy, innerR - 1, innerR - 7, hStart, hEnd)
          : '';
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
      if (dragSelection && dragSelection.chainId === id) return;
      tooltip.classList.remove('clickable');
      tooltip.style.opacity = '0';
    });
    arcsGroup.appendChild(path);
    if (hasActiveMask && !fullyActive) {
      const span = chain.end - chain.start;
      for (const [startRes, endRes] of activeRanges) {
        const segStart = chain.start + ((startRes - 1) / chain.length) * span;
        const segEnd = chain.start + (endRes / chain.length) * span;
        if (segEnd <= segStart) continue;
        const activePath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        activePath.setAttribute('d', arcPath(cx, cy, outerR, innerR, segStart, segEnd));
        activePath.setAttribute('fill', chain.color);
        activePath.setAttribute('opacity', '0.9');
        activePath.setAttribute('pointer-events', 'none');
        arcsGroup.appendChild(activePath);
      }
    }

    const mid = (chain.start + chain.end) / 2;
    const reverse = Math.sin(mid) > 0;
    const labelPathOffset = reverse ? 18 : 10;
    const labelExportOffset = reverse ? 26 : 18;
    const labelPathId = `label-path-${id.replace(/[^A-Za-z0-9_-]/g, '_')}-${Math.round(chain.start * 1000)}`;
    const labelPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    labelPath.setAttribute('id', labelPathId);
    labelPath.setAttribute('d', arcLinePath(cx, cy, outerR + labelPathOffset, chain.start, chain.end, reverse));
    labelDefs.appendChild(labelPath);
    const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    const [labelX, labelY] = polar(cx, cy, outerR + labelExportOffset, mid);
    label.setAttribute('fill', '#0f172a');
    label.setAttribute('font-size', labelSize.value);
    label.setAttribute('font-weight', '600');
    label.setAttribute('text-anchor', 'middle');
    label.setAttribute('dominant-baseline', 'middle');
    label.setAttribute('data-export-label', 'arc');
    label.setAttribute('data-export-x', `${labelX}`);
    label.setAttribute('data-export-y', `${labelY}`);
    label.setAttribute('data-export-rotation', `${(mid * 180 / Math.PI) + 90 + (reverse ? 180 : 0)}`);
    const textPath = document.createElementNS('http://www.w3.org/2000/svg', 'textPath');
    textPath.setAttributeNS('http://www.w3.org/1999/xlink', 'xlink:href', `#${labelPathId}`);
    textPath.setAttribute('href', `#${labelPathId}`);
    textPath.setAttribute('startOffset', '50%');
    textPath.textContent = chainDisplayName.get(id) || id;
    label.appendChild(textPath);
    labelsGroup.appendChild(label);
  }

  hoverArcPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
  overlayGroup.appendChild(hoverArcPath);

  const openSelectionMenuAt = (event, chainId, pos) => {
    event.stopPropagation();
    const clicked = arcSelections.filter((x) => {
      if (x.chainId !== chainId) return false;
      const [a, b] = normalizeRange(x.startPos, x.endPos);
      return pos >= a && pos <= b;
    }).sort((a, b) => b.id - a.id);
    if (!clicked.length) return;
    const thresholdValue = Number(threshold.value);
    const visibleChains = new Set(currentOrder);
    const blocks = clicked.map((item) => {
      const [a, b] = normalizeRange(item.startPos, item.endPos);
      const seqText = selectionSequenceText(item);
      const seqEnabled = seqText.length > 0;
      const canTrim = trimSelectionToContacts(item, thresholdValue, visibleChains) !== null;
      const maxLen = chainMap.get(item.chainId).length;
      return `
        <div class="menu-block">
          <div class="menu-title">${item.chainId}:${a}-${b}</div>
          <div class="menu-row">
            <button data-action="callout" data-id="${item.id}" ${seqEnabled ? '' : 'disabled'}>Sequence callout</button>
            <button data-action="comment" data-id="${item.id}">Comment</button>
            <button data-action="autotrim" data-id="${item.id}" ${canTrim ? '' : 'disabled'}>Autotrim</button>
            <button data-action="clear" data-id="${item.id}">Clear</button>
          </div>
          <div class="menu-row">
            <span style="font-size:12px;color:#334155;">Edit</span>
            <input type="number" data-action="start" data-id="${item.id}" value="${a}" min="1" max="${maxLen}" style="width:64px;">
            <input type="range" data-action="startRange" data-id="${item.id}" value="${a}" min="1" max="${maxLen}" style="width:90px;">
          </div>
          <div class="menu-row">
            <span style="font-size:12px;color:#334155;">&nbsp;</span>
            <input type="number" data-action="end" data-id="${item.id}" value="${b}" min="1" max="${maxLen}" style="width:64px;">
            <input type="range" data-action="endRange" data-id="${item.id}" value="${b}" min="1" max="${maxLen}" style="width:90px;">
            <button data-action="applyEdit" data-id="${item.id}">Apply</button>
          </div>
          <div class="menu-row">
            <span style="font-size:12px;color:#334155;">Color</span>
            <input type="color" data-action="color" data-id="${item.id}" value="${item.color || chainMap.get(item.chainId).color}">
          </div>
          <div class="menu-row">
            <span style="font-size:12px;color:#334155;">Link color</span>
            <input type="color" data-action="linkColor" data-id="${item.id}" value="${item.linkColor || '#0f172a'}">
            <button data-action="clearLinkColor" data-id="${item.id}">Default</button>
          </div>
        </div>
      `;
    }).join('');
    selectionMenu.innerHTML = `<div class="menu-title">Selections (${clicked.length})</div>${blocks}`;
    const rect = plotWrap.getBoundingClientRect();
    selectionMenu.style.left = `${clamp(event.clientX - rect.left + 10, 6, rect.width - 250)}px`;
    selectionMenu.style.top = `${clamp(event.clientY - rect.top + 10, 6, rect.height - 220)}px`;
    selectionMenu.style.display = 'block';

    selectionMenu.querySelectorAll('[data-action]').forEach((node) => {
      if (node.tagName === 'BUTTON') {
        node.addEventListener('click', (evt) => {
        evt.stopPropagation();
        const idVal = Number(evt.currentTarget.getAttribute('data-id'));
        const act = evt.currentTarget.getAttribute('data-action');
        const selObj = arcSelections.find((x) => x.id === idVal);
        if (!selObj) return;
        if (act === 'callout') {
          const text = selectionSequenceText(selObj);
          if (!text) return;
          selObj.callout = selObj.callout || { kind: 'sequence', text, html: '', x: null, y: null, w: 220, h: 80, title: '' };
          selObj.callout.kind = 'sequence';
          selObj.callout.text = text;
          selObj.callout.html = '';
          selObj.callout.title = calloutRangeTitle(selObj);
        } else if (act === 'comment') {
          selObj.callout = selObj.callout || { kind: 'comment', text: '', html: '', x: null, y: null, w: 220, h: 80, title: '' };
          selObj.callout.kind = 'comment';
          selObj.callout.text = selObj.comment || selObj.callout.text || '';
          selObj.callout.html = selObj.callout.html || '';
          selObj.callout.title = '';
        } else if (act === 'autotrim') {
          const trimmed = trimSelectionToContacts(selObj, thresholdValue, visibleChains);
          if (trimmed) {
            selObj.startPos = trimmed[0];
            selObj.endPos = trimmed[1];
            if (selObj.callout && selObj.callout.kind === 'sequence') {
              selObj.callout.title = calloutRangeTitle(selObj);
            }
          }
        } else if (act === 'applyEdit') {
          const maxLen = chainMap.get(selObj.chainId).length;
          const s = clamp(Number(selObj._editStart ?? selObj.startPos), 1, maxLen);
          const e = clamp(Number(selObj._editEnd ?? selObj.endPos), 1, maxLen);
          if (s <= e) {
            selObj.startPos = s;
            selObj.endPos = e;
            if (selObj.callout && selObj.callout.kind === 'sequence') {
              selObj.callout.title = calloutRangeTitle(selObj);
            }
          }
        } else if (act === 'clear') {
          arcSelections = arcSelections.filter((x) => x.id !== idVal);
        } else if (act === 'clearLinkColor') {
          selObj.linkColor = null;
        }
        if (act !== 'autotrim' && act !== 'applyEdit') {
          closeSelectionMenu();
        }
        safeRender();
        });
      }
      if (node.tagName === 'INPUT') {
        node.addEventListener('input', (evt) => {
          evt.stopPropagation();
          const idVal = Number(evt.currentTarget.getAttribute('data-id'));
          const act = evt.currentTarget.getAttribute('data-action');
          const selObj = arcSelections.find((x) => x.id === idVal);
          if (!selObj) return;
          if (act === 'color') {
            selObj.color = evt.currentTarget.value;
          } else if (act === 'linkColor') {
            selObj.linkColor = evt.currentTarget.value;
          } else if (act === 'start' || act === 'startRange') {
            const v = clamp(Number(evt.currentTarget.value), 1, chainMap.get(selObj.chainId).length);
            selObj._editStart = v;
            selectionMenu.querySelector(`input[data-action="start"][data-id="${idVal}"]`).value = v;
            selectionMenu.querySelector(`input[data-action="startRange"][data-id="${idVal}"]`).value = v;
          } else if (act === 'end' || act === 'endRange') {
            const v = clamp(Number(evt.currentTarget.value), 1, chainMap.get(selObj.chainId).length);
            selObj._editEnd = v;
            selectionMenu.querySelector(`input[data-action="end"][data-id="${idVal}"]`).value = v;
            selectionMenu.querySelector(`input[data-action="endRange"][data-id="${idVal}"]`).value = v;
          }
          safeRender();
        });
      }
    });
  };

  const renderSel = (sel, idxInStack = 0, preview = false) => {
    const chain = chainAngles.get(sel.chainId);
    if (!chain) return;
    const [ds, de] = selectionDisplayBounds(sel, chain);
    const span = chain.end - chain.start;
    const segStart = chain.start + ((ds - 1) / chain.length) * span;
    const segEnd = chain.start + (de / chain.length) * span;
    if (segEnd <= segStart) return;
    const color = sel.color || chain.color;

    const fillPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    fillPath.setAttribute('d', arcPath(cx, cy, outerR, innerR, segStart, segEnd));
    fillPath.setAttribute('fill', color);
    fillPath.setAttribute('opacity', preview ? '0.20' : '0.28');
    fillPath.setAttribute('pointer-events', 'none');
    selectionFillGroup.appendChild(fillPath);

    const strokePath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    strokePath.setAttribute('d', arcPath(cx, cy, outerR + 0.5, innerR - 0.5, segStart, segEnd));
    strokePath.setAttribute('fill', 'none');
    strokePath.setAttribute('stroke', lightenHexColor(color, 0.5));
    strokePath.setAttribute('stroke-width', preview ? '2.0' : '3.6');
    strokePath.setAttribute('stroke-dasharray', 'none');
    strokePath.setAttribute('stroke-linecap', 'round');
    strokePath.setAttribute('pointer-events', 'none');
    selectionStrokeGroup.appendChild(strokePath);

    if (sel.callout && !preview) {
      const mid = (ds + de) / 2;
      const midAngle = chain.start + ((mid - 0.5) / chain.length) * (chain.end - chain.start);
      const [axRaw, ayRaw] = polar(cx, cy, outerR + 92 + idxInStack * 28, midAngle);
      const callout = sel.callout;
      if (callout.x == null || callout.y == null) {
        callout.w = callout.w || 220;
        callout.h = callout.h || 80;
        callout.x = clamp(axRaw - callout.w / 2, 0, CANVAS_W - callout.w);
        callout.y = clamp(ayRaw - callout.h / 2, 0, CANVAS_H - callout.h);
      }
      const width = clamp(callout.w || 220, 80, 760);
      const seqDefault = selectionSequenceText(sel) || `${sel.chainId}:${sel.startPos}-${sel.endPos}`;
      if (callout.kind === 'sequence' && !callout.text) {
        callout.text = seqDefault;
      }
      if (callout.kind === 'sequence' && !callout.title) {
        callout.title = calloutRangeTitle(sel);
      }
      const textVal = callout.kind === 'comment'
        ? (callout.text || '')
        : (callout.text || seqDefault);
      const htmlVal = callout.html && callout.html.length
        ? callout.html
        : escapeHtml(textVal).replace(/\\n/g, '<br>');
      const titleText = callout.kind === 'sequence'
        ? (callout.title || calloutRangeTitle(sel))
        : '';
      const titleBarH = titleText ? 16 : 4;
      const minH = 40;
      const height = clamp(Math.max(callout.h || minH, minH), minH, 620);
      const boxX = clamp(callout.x, 0, CANVAS_W - width);
      const boxY = clamp(callout.y, 0, CANVAS_H - height);
      callout.x = boxX;
      callout.y = boxY;
      callout.w = width;
      callout.h = height;

      const anchorX = cx + (outerR + 8) * Math.cos(midAngle);
      const anchorY = cy + (outerR + 8) * Math.sin(midAngle);
      const [toX, toY] = nearestPointOnRect(anchorX, anchorY, boxX, boxY, width, height);
      const leader = document.createElementNS('http://www.w3.org/2000/svg', 'line');
      leader.setAttribute('x1', anchorX);
      leader.setAttribute('y1', anchorY);
      leader.setAttribute('x2', toX);
      leader.setAttribute('y2', toY);
      leader.setAttribute('stroke', color);
      leader.setAttribute('stroke-width', '1.4');
      calloutGroup.appendChild(leader);

      const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      rect.setAttribute('x', boxX);
      rect.setAttribute('y', boxY);
      rect.setAttribute('rx', '6');
      rect.setAttribute('ry', '6');
      rect.setAttribute('width', `${width}`);
      rect.setAttribute('height', `${height}`);
      rect.setAttribute('fill', 'rgba(255,255,255,0.92)');
      rect.setAttribute('stroke', color);
      rect.setAttribute('stroke-width', '1');
      rect.setAttribute('data-callout-box', `${sel.id}`);
      calloutGroup.appendChild(rect);

      const dragHandle = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      dragHandle.setAttribute('x', boxX + 2);
      dragHandle.setAttribute('y', boxY + 2);
      dragHandle.setAttribute('width', `${Math.max(12, width - 4)}`);
      dragHandle.setAttribute('height', '12');
      dragHandle.setAttribute('fill', 'rgba(148,163,184,0.18)');
      dragHandle.setAttribute('data-export-ignore', 'true');
      dragHandle.style.cursor = 'move';
      dragHandle.addEventListener('pointerdown', (ev) => {
        ev.stopPropagation();
        const p = svgPoint(ev);
        activeCalloutDrag = { id: sel.id, mode: 'move', sx: p.x, sy: p.y, ox: boxX, oy: boxY };
      });
      calloutGroup.appendChild(dragHandle);

      if (titleText) {
        const title = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        title.setAttribute('x', boxX + 6);
        title.setAttribute('y', boxY + 11);
        title.setAttribute('fill', '#334155');
        title.setAttribute('font-size', '11');
        title.setAttribute('font-weight', '600');
        title.setAttribute('font-family', 'IBM Plex Sans, sans-serif');
        title.textContent = titleText;
        title.setAttribute('pointer-events', 'none');
        title.setAttribute('data-callout-title', `${sel.id}`);
        calloutGroup.appendChild(title);
      }

      const fo = document.createElementNS('http://www.w3.org/2000/svg', 'foreignObject');
      fo.setAttribute('x', boxX + 2);
      fo.setAttribute('y', boxY + titleBarH);
      fo.setAttribute('width', `${Math.max(4, width - 4)}`);
      fo.setAttribute('height', `${Math.max(4, height - titleBarH - 2)}`);
      fo.setAttribute('data-export-text', encodeURIComponent(sel.callout.text || ''));
      fo.setAttribute('data-callout-fo', `${sel.id}`);
      const div = document.createElement('div');
      div.setAttribute('xmlns', 'http://www.w3.org/1999/xhtml');
      div.setAttribute(
        'style',
        'width:100%;height:100%;overflow:auto;white-space:pre-wrap;word-break:break-word;'
        + 'font:12px \"IBM Plex Sans\",sans-serif;color:#0f172a;outline:none;padding:0;margin:0;'
      );
      div.contentEditable = 'true';
      div.innerHTML = htmlVal;
      div.addEventListener('pointerdown', (ev) => ev.stopPropagation());
      div.addEventListener('input', () => {
        sel.callout.html = sanitizeEditableHtml(div.innerHTML);
        sel.callout.text = div.innerText;
        fo.setAttribute('data-export-text', encodeURIComponent(sel.callout.text || ''));
        if (sel.callout.kind !== 'comment' && sel.callout.text !== selectionSequenceText(sel)) {
          sel.callout.kind = 'comment';
          sel.callout.title = '';
          sel.comment = sel.callout.text;
        }
        if (sel.callout.kind === 'comment') {
          sel.comment = sel.callout.text;
        }
      });
      div.addEventListener('blur', () => safeRender());
      fo.appendChild(div);
      calloutGroup.appendChild(fo);

      const corners = [
        ['nw', boxX, boxY],
        ['ne', boxX + width, boxY],
        ['sw', boxX, boxY + height],
        ['se', boxX + width, boxY + height],
      ];
      for (const [mode, hx, hy] of corners) {
        const h = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        h.setAttribute('x', hx - 6);
        h.setAttribute('y', hy - 6);
        h.setAttribute('width', '12');
        h.setAttribute('height', '12');
        h.setAttribute('fill', '#ffffff');
        h.setAttribute('fill-opacity', '0');
        h.setAttribute('stroke', color);
        h.setAttribute('stroke-width', '1');
        h.setAttribute('stroke-opacity', '0');
        h.setAttribute('data-export-ignore', 'true');
        h.style.cursor = `${mode}-resize`;
        h.addEventListener('pointerdown', (ev) => {
          ev.stopPropagation();
          const p = svgPoint(ev);
          activeCalloutDrag = { id: sel.id, mode, sx: p.x, sy: p.y, ox: boxX, oy: boxY, ow: width, oh: height };
        });
        calloutGroup.appendChild(h);
      }
    }
  };

  arcSelections.forEach((sel, idx) => renderSel(sel, idx, false));
  if (dragSelection && dragSelection.moved) {
    renderSel(
      {
        chainId: dragSelection.chainId,
        startPos: dragSelection.startPos,
        endPos: dragSelection.currentPos,
      },
      0,
      true
    );
  }

  let regionTargetInfo = null;
  const regionTarget = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
  regionTarget.setAttribute('r', '8');
  regionTarget.setAttribute('fill', '#1d4ed8');
  regionTarget.setAttribute('opacity', '0.9');
  regionTarget.style.display = 'none';
  regionTarget.style.cursor = 'pointer';
  regionTarget.addEventListener('click', (event) => {
    if (!regionTargetInfo) return;
    openSelectionMenuAt(event, regionTargetInfo.chainId, regionTargetInfo.pos);
  });
  targetGroup.appendChild(regionTarget);

  const pickRegionTarget = (event) => {
    const { x, y } = svgPoint(event);
    const dx = x - cx;
    const dy = y - cy;
    const r = Math.hypot(dx, dy);
    if (Math.abs(r - outerR) > 18) return null;
    let ang = Math.atan2(dy, dx);
    if (ang < 0) ang += Math.PI * 2;
    const sorted = arcSelections.slice().sort((a, b) => b.id - a.id);
    for (const sel of sorted) {
      const chain = chainAngles.get(sel.chainId);
      if (!chain) continue;
      const adj = unwrapAngleForSpan(ang, chain.start, chain.end);
      if (adj < chain.start || adj > chain.end) continue;
      const rawFrac = clamp((adj - chain.start) / (chain.end - chain.start), 0, 1);
      let idx = clamp(Math.floor(rawFrac * chain.length) + 1, 1, chain.length);
      if (chain.flip) idx = chain.length - idx + 1;
      const [s, e] = normalizeRange(sel.startPos, sel.endPos);
      if (idx < s || idx > e) continue;
      const tx = cx + (outerR + 10) * Math.cos(ang);
      const ty = cy + (outerR + 10) * Math.sin(ang);
      return { chainId: sel.chainId, pos: idx, x: tx, y: ty };
    }
    return null;
  };

  svg.onmousemove = (event) => {
    const found = pickRegionTarget(event);
    if (!found) {
      regionTargetInfo = null;
      regionTarget.style.display = 'none';
      return;
    }
    regionTargetInfo = found;
    regionTarget.setAttribute('cx', found.x);
    regionTarget.setAttribute('cy', found.y);
    regionTarget.style.display = '';
  };
  svg.onmouseleave = () => {
    regionTargetInfo = null;
    regionTarget.style.display = 'none';
  };

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
    let overrideColor = null;
    for (const sel of arcSelections) {
      if (!sel.linkColor) continue;
      if (linkMatchesSelection(link, sel)) overrideColor = sel.linkColor;
    }
    const color = overrideColor || mixColor(a.color, b.color);

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

  svg.appendChild(labelDefs);
  svg.appendChild(linksGroup);
  svg.appendChild(arcsGroup);
  svg.appendChild(selectionFillGroup);
  svg.appendChild(selectionStrokeGroup);
  svg.appendChild(overlayGroup);
  svg.appendChild(calloutGroup);
  svg.appendChild(labelsGroup);
  svg.appendChild(targetGroup);
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
    const lockNote = lockedBottomChain ? ` Locked bottom: ${lockedBottomChain}.` : '';
    renderStatus.textContent = `Rendered ${stats.linkCount} links across ${stats.chainCount} chains.${filterNote}${lockNote}`;
  } catch (err) {
    renderStatus.textContent = `Render error: ${err.message}`;
    console.error(err);
  }
}

function buildSessionState() {
  return {
    schema: 'contacts-circos-session-v2',
    controls: {
      plotTitle: plotTitleInput.value,
      threshold: Number(threshold.value),
      maxWidth: Number(maxWidth.value),
      angle: Number(angle.value),
      labelSize: Number(labelSize.value),
      countMode,
      orderText: orderInput.value,
      lockedBottomChain
    },
    chain: {
      contactVisibility: Object.fromEntries(contactVisibility.entries()),
      flip: Object.fromEntries(chainFlip.entries()),
      displayName: Object.fromEntries(chainDisplayName.entries())
    },
    selections: arcSelections.map((sel) => ({
      id: sel.id,
      chainId: sel.chainId,
      startPos: sel.startPos,
      endPos: sel.endPos,
      color: sel.color ?? null,
      linkColor: sel.linkColor ?? null,
      comment: sel.comment ?? '',
      callout: sel.callout
        ? {
            kind: sel.callout.kind || 'comment',
            title: sel.callout.title || '',
            text: sel.callout.text || '',
            html: sel.callout.html || '',
            x: sel.callout.x ?? null,
            y: sel.callout.y ?? null,
            w: sel.callout.w ?? null,
            h: sel.callout.h ?? null
          }
        : null
    }))
  };
}

function applySessionState(state) {
  if (!state || !['contacts-circos-session-v1', 'contacts-circos-session-v2'].includes(state.schema)) {
    throw new Error('Unsupported session file');
  }
  const controls = state.controls || {};
  plotTitleInput.value = typeof controls.plotTitle === 'string' ? controls.plotTitle : (DATA.title || '');
  threshold.value = clamp(Number(controls.threshold ?? threshold.value), Number(threshold.min), Number(threshold.max));
  thresholdInput.value = threshold.value;
  maxWidth.value = clamp(Number(controls.maxWidth ?? maxWidth.value), Number(maxWidth.min), Number(maxWidth.max));
  maxWidthInput.value = maxWidth.value;
  angle.value = clamp(Number(controls.angle ?? angle.value), Number(angle.min), Number(angle.max));
  angleInput.value = angle.value;
  labelSize.value = clamp(Number(controls.labelSize ?? labelSize.value), Number(labelSize.min), Number(labelSize.max));
  labelSizeInput.value = labelSize.value;
  if (controls.countMode === 'atom' || controls.countMode === 'residue') {
    countMode = controls.countMode;
  }
  if (typeof controls.lockedBottomChain === 'string' && chainMap.has(controls.lockedBottomChain)) {
    lockedBottomChain = controls.lockedBottomChain;
  } else {
    lockedBottomChain = null;
  }
  document.querySelectorAll('input[name="countMode"]').forEach((radio) => {
    radio.checked = radio.value === countMode;
  });
  orderInput.value = typeof controls.orderText === 'string' ? controls.orderText : defaultOrder.join(', ');

  const chainState = state.chain || {};
  const vis = chainState.contactVisibility || {};
  const flip = chainState.flip || {};
  const names = chainState.displayName || {};
  for (const chain of DATA.chains) {
    const id = chain.id;
    contactVisibility.set(id, vis[id] !== false);
    chainFlip.set(id, flip[id] === true);
    chainDisplayName.set(id, typeof names[id] === 'string' && names[id].trim() ? names[id] : id);
  }

  const loadedSelections = [];
  const inputSelections = Array.isArray(state.selections) ? state.selections : [];
  for (const raw of inputSelections) {
    if (!raw || !chainMap.has(raw.chainId)) continue;
    const maxLen = chainMap.get(raw.chainId).length;
    const s = clamp(Number(raw.startPos), 1, maxLen);
    const e = clamp(Number(raw.endPos), 1, maxLen);
    if (!Number.isFinite(s) || !Number.isFinite(e)) continue;
    const sel = {
      id: Number(raw.id) || nextSelectionId++,
      chainId: raw.chainId,
      startPos: Math.min(s, e),
      endPos: Math.max(s, e),
      color: raw.color || null,
      linkColor: raw.linkColor || null,
      comment: typeof raw.comment === 'string' ? raw.comment : '',
      callout: null
    };
    if (raw.callout && typeof raw.callout === 'object') {
      sel.callout = {
        kind: raw.callout.kind === 'sequence' ? 'sequence' : 'comment',
        title: typeof raw.callout.title === 'string' ? raw.callout.title : '',
        text: typeof raw.callout.text === 'string' ? raw.callout.text : '',
        html: typeof raw.callout.html === 'string' ? raw.callout.html : '',
        x: Number.isFinite(raw.callout.x) ? Number(raw.callout.x) : null,
        y: Number.isFinite(raw.callout.y) ? Number(raw.callout.y) : null,
        w: Number.isFinite(raw.callout.w) ? Number(raw.callout.w) : null,
        h: Number.isFinite(raw.callout.h) ? Number(raw.callout.h) : null
      };
    }
    loadedSelections.push(sel);
  }
  arcSelections = loadedSelections;
  nextSelectionId = Math.max(1, ...arcSelections.map((x) => Number(x.id) || 0)) + 1;
  closeSelectionMenu();
  filterChain = null;
  filterResidue = null;
  safeRender();
}

function syncRange(range, input, onUserChange = null) {
  range.addEventListener('input', () => {
    input.value = range.value;
    if (typeof onUserChange === 'function') onUserChange();
    safeRender();
  });
  input.addEventListener('change', () => {
    range.value = clamp(Number(input.value), Number(range.min), Number(range.max));
    input.value = range.value;
    if (typeof onUserChange === 'function') onUserChange();
    safeRender();
  });
}

syncRange(threshold, thresholdInput);
syncRange(maxWidth, maxWidthInput);
syncRange(angle, angleInput, () => {
  if (lockedBottomChain) lockedBottomChain = null;
});
syncRange(labelSize, labelSizeInput);
plotTitleInput.addEventListener('input', () => safeRender());

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
  let createdSelection = false;
  if (activeCalloutDrag) {
    activeCalloutDrag = null;
  }
  if (dragSelection) {
    if (dragSelection.moved) {
      const [s, e] = normalizeRange(dragSelection.startPos, dragSelection.currentPos);
      arcSelections.push({
        id: nextSelectionId++,
        chainId: dragSelection.chainId,
        startPos: s,
        endPos: e,
        color: null,
        callout: false
      });
      createdSelection = true;
      suppressNextClick = true;
    }
    dragSelection = null;
  }
  if (filterChain || filterResidue) {
    filterChain = null;
    filterResidue = null;
    safeRender();
    return;
  }
  if (createdSelection) {
    safeRender();
    return;
  }
  if (suppressNextClick) {
    suppressNextClick = false;
    return;
  }
});

window.addEventListener('pointermove', (event) => {
  if (!activeCalloutDrag) return;
  const p = svgPoint(event);
  const sel = arcSelections.find((x) => x.id === activeCalloutDrag.id);
  if (!sel || !sel.callout) return;
  const dx = p.x - activeCalloutDrag.sx;
  const dy = p.y - activeCalloutDrag.sy;
  const c = sel.callout;
  if (activeCalloutDrag.mode === 'move') {
    c.x = activeCalloutDrag.ox + dx;
    c.y = activeCalloutDrag.oy + dy;
  } else {
    let x = activeCalloutDrag.ox;
    let y = activeCalloutDrag.oy;
    let w = activeCalloutDrag.ow;
    let h = activeCalloutDrag.oh;
    if (activeCalloutDrag.mode.includes('e')) w = activeCalloutDrag.ow + dx;
    if (activeCalloutDrag.mode.includes('s')) h = activeCalloutDrag.oh + dy;
    if (activeCalloutDrag.mode.includes('w')) { w = activeCalloutDrag.ow - dx; x = activeCalloutDrag.ox + dx; }
    if (activeCalloutDrag.mode.includes('n')) { h = activeCalloutDrag.oh - dy; y = activeCalloutDrag.oy + dy; }
    c.x = x; c.y = y; c.w = w; c.h = h;
  }
  safeRender();
});

plotWrap.addEventListener('click', (event) => {
  if (suppressNextClick) {
    suppressNextClick = false;
    return;
  }
  if (!selectionMenu.contains(event.target)) {
    closeSelectionMenu();
  }
});

saveSessionBtn.addEventListener('click', () => {
  const blob = new Blob([JSON.stringify(buildSessionState(), null, 2)], {type: 'application/json'});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'contacts_circos_session.json';
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
});

loadSessionBtn.addEventListener('click', () => {
  loadSessionFile.value = '';
  loadSessionFile.click();
});

loadSessionFile.addEventListener('change', async () => {
  const file = loadSessionFile.files && loadSessionFile.files[0];
  if (!file) return;
  try {
    const text = await file.text();
    const obj = JSON.parse(text);
    applySessionState(obj);
    renderStatus.textContent = 'Session loaded.';
  } catch (err) {
    renderStatus.textContent = `Session load error: ${err.message}`;
  }
});

document.getElementById('downloadSvg').addEventListener('click', () => {
  const ns = 'http://www.w3.org/2000/svg';
  const exportSvg = svg.cloneNode(true);
  exportSvg.setAttribute('xmlns', ns);
  exportSvg.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
  exportSvg.querySelectorAll('[data-export-ignore="true"]').forEach((node) => node.remove());
  exportSvg.querySelectorAll('text[data-export-label="arc"]').forEach((label) => {
    const textPath = label.querySelector('textPath');
    if (!textPath) return;
    const flat = document.createElementNS(ns, 'text');
    flat.setAttribute('x', label.getAttribute('data-export-x') || '0');
    flat.setAttribute('y', label.getAttribute('data-export-y') || '0');
    flat.setAttribute('fill', label.getAttribute('fill') || '#0f172a');
    flat.setAttribute('font-size', label.getAttribute('font-size') || '12');
    flat.setAttribute('font-weight', label.getAttribute('font-weight') || '600');
    flat.setAttribute('font-family', label.getAttribute('font-family') || 'IBM Plex Sans, sans-serif');
    flat.setAttribute('text-anchor', 'middle');
    flat.setAttribute('dominant-baseline', 'middle');
    const x = label.getAttribute('data-export-x') || '0';
    const y = label.getAttribute('data-export-y') || '0';
    const rot = label.getAttribute('data-export-rotation') || '0';
    flat.setAttribute('transform', `rotate(${rot} ${x} ${y})`);
    flat.textContent = textPath.textContent || '';
    label.parentNode.insertBefore(flat, label.nextSibling);
    label.remove();
  });
  exportSvg.querySelectorAll('foreignObject').forEach((fo) => {
    const x = Number(fo.getAttribute('x') || 0);
    const y = Number(fo.getAttribute('y') || 0);
    const w = Number(fo.getAttribute('width') || 120);
    const encoded = fo.getAttribute('data-export-text') || '';
    const textVal = decodeURIComponent(encoded);
    const lines = wrapLines(textVal, Math.max(10, w - 4));
    const boxId = fo.getAttribute('data-callout-fo');
    const box = boxId ? exportSvg.querySelector(`rect[data-callout-box="${boxId}"]`) : null;
    const title = boxId ? exportSvg.querySelector(`text[data-callout-title="${boxId}"]`) : null;
    const titleBarH = title ? 16 : 4;
    const neededHeight = Math.max(40, titleBarH + 4 + Math.max(1, lines.length) * 14);
    if (box) {
      box.setAttribute('height', `${neededHeight}`);
    }
    const text = document.createElementNS(ns, 'text');
    text.setAttribute('fill', '#0f172a');
    text.setAttribute('font-size', '12');
    text.setAttribute('font-family', 'IBM Plex Sans, sans-serif');
    lines.forEach((line, idx) => {
      const t = document.createElementNS(ns, 'tspan');
      t.setAttribute('x', `${x + 2}`);
      t.setAttribute('y', `${y + 12 + idx * 14}`);
      t.textContent = line;
      text.appendChild(t);
    });
    fo.parentNode.insertBefore(text, fo.nextSibling);
    fo.remove();
  });
  const serializer = new XMLSerializer();
  const svgString = serializer.serializeToString(exportSvg);
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

saveHtml.addEventListener('click', () => {
  const state = buildSessionState();
  const stateJson = JSON.stringify(state);
  const docClone = document.documentElement.cloneNode(true);
  docClone.setAttribute('data-saved-session', encodeURIComponent(stateJson));
  const html = '<!doctype html>\\n' + docClone.outerHTML;
  const blob = new Blob([html], {type: 'text/html;charset=utf-8'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'contacts_circos.html';
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(a.href);
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

function applyEmbeddedSavedSessionIfPresent() {
  const attr = document.documentElement.getAttribute('data-saved-session');
  if (attr) {
    try {
      const obj = JSON.parse(decodeURIComponent(attr));
      applySessionState(obj);
      renderStatus.textContent = 'Session restored from saved HTML.';
      return true;
    } catch (err) {
      renderStatus.textContent = `Saved HTML load error: ${err.message}`;
      console.error(err);
      return false;
    }
  }
  const node = document.getElementById('savedSessionData');
  if (!node) return false;
  try {
    const obj = JSON.parse(node.textContent || '{}');
    applySessionState(obj);
    renderStatus.textContent = 'Session restored from saved HTML.';
    return true;
  } catch (err) {
    renderStatus.textContent = `Saved HTML load error: ${err.message}`;
    console.error(err);
    return false;
  }
}

if (!applyEmbeddedSavedSessionIfPresent()) {
  safeRender();
}
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
    output_path.write_text(html, encoding="utf-8")


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
    parser.add_argument(
        "--dna-mode",
        choices=["auto", "merge", "split"],
        default="auto",
        help="DNA display mode: auto/merge collapses to one DNA arc, split keeps individual chains.",
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
        chain_res,
        chain_resnums,
        dna_chains,
        args.dna_reverse,
        args.dna_mismatches,
        dna_mode=args.dna_mode,
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
        {},
        output_path,
    )
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
