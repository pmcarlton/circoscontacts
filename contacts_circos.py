#!/usr/bin/env python3
"""
Generate an interactive circos-style contact plot (HTML + SVG export)
from ChimeraX contact files and corresponding mmCIFs.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import shlex
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
TOOL_VERSION = "0.4.20"


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


def parse_contacts(
    contact_paths: Iterable[Path],
) -> Tuple[
    Counter,
    List[set[Tuple[str, int, str, int]]],
    Dict[Tuple[str, int, str, int], set[str]],
    Dict[Tuple[str, int, str, int], Dict[str, set[Tuple[str, str, int]]]],
]:
    def extract_chain(token: str) -> str:
        if token.startswith("/"):
            return token[1:]
        if "/" in token:
            return token.rsplit("/", 1)[-1]
        return token

    counts_atom: Counter = Counter()
    residue_pairs_by_file: List[set[Tuple[str, int, str, int]]] = []
    contact_models: Dict[Tuple[str, int, str, int], set[str]] = defaultdict(set)
    contact_sources: Dict[Tuple[str, int, str, int], Dict[str, set[Tuple[str, str, int]]]] = defaultdict(
        lambda: {"a": set(), "b": set()}
    )
    for path in contact_paths:
        per_file_pairs = set()
        with path.open() as fh:
            for line in fh:
                parts = line.split()
                chain1 = chain2 = None
                res1 = res2 = None
                model1 = model2 = ""

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
                    parsed_model1 = parts[1].split("/", 1)[0]
                    parsed_model2 = parts[6].split("/", 1)[0]
                    model1 = parsed_model1 if parsed_model1.startswith("#") else ""
                    model2 = parsed_model2 if parsed_model2.startswith("#") else ""
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
                if model1 and model2:
                    if model1 == model2:
                        contact_models[key].add(model1)
                    else:
                        contact_models[key].update((model1, model2))
                if key == (chain1, res1, chain2, res2):
                    contact_sources[key]["a"].add((model1, chain1, res1))
                    contact_sources[key]["b"].add((model2, chain2, res2))
                else:
                    contact_sources[key]["a"].add((model2, chain2, res2))
                    contact_sources[key]["b"].add((model1, chain1, res1))
        residue_pairs_by_file.append(per_file_pairs)
    return counts_atom, residue_pairs_by_file, contact_models, contact_sources


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
    residue_pairs_by_file: List[set[Tuple[str, int, str, int]]],
    contact_models: Dict[Tuple[str, int, str, int], set[str]],
    contact_sources: Dict[Tuple[str, int, str, int], Dict[str, set[Tuple[str, str, int]]]],
    pos_map: Dict[str, Dict[int, int]],
    display_chain_of: Dict[str, str],
) -> Tuple[List[Dict[str, int | str | list]], int, int]:
    def model_sort_key(label: str) -> Tuple[int, object]:
        if label.startswith("#") and label[1:].isdigit():
            return (0, int(label[1:]))
        return (1, label)

    agg_atom: Counter = Counter()
    agg_residue: Counter = Counter()
    agg_models: Dict[Tuple[str, int, str, int], set[str]] = defaultdict(set)
    agg_source_entries: Dict[Tuple[str, int, str, int], Dict[str, set[Tuple[str, str, int]]]] = defaultdict(
        lambda: {"a": set(), "b": set()}
    )
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
        agg_models[key].update(contact_models.get((chain1, res1, chain2, res2), set()))
        for side in ("a", "b"):
            agg_source_entries[key][side].update(contact_sources.get((chain1, res1, chain2, res2), {}).get(side, set()))
        agg_atom[key] += counts_atom[(chain1, res1, chain2, res2)]

    for per_file_pairs in residue_pairs_by_file:
        mapped_keys = set()
        for chain1, res1, chain2, res2 in per_file_pairs:
            if chain1 not in pos_map or chain2 not in pos_map:
                continue
            if res1 not in pos_map[chain1] or res2 not in pos_map[chain2]:
                continue
            a_chain = display_chain_of.get(chain1, chain1)
            b_chain = display_chain_of.get(chain2, chain2)
            a_pos = pos_map[chain1][res1]
            b_pos = pos_map[chain2][res2]
            mapped_keys.add(canonical_contact(a_chain, a_pos, b_chain, b_pos))
        for key in mapped_keys:
            agg_residue[key] += 1
    if skipped:
        sys.stderr.write(f"Warning: skipped {skipped} contacts without CIF mapping\n")

    max_count_atom = max(agg_atom.values()) if agg_atom else 1
    max_count_residue = max(agg_residue.values()) if agg_residue else 1
    contacts = []
    for (a, pa, b, pb), c in agg_atom.items():
        src = sources.get((a, pa, b, pb), {"a": set(), "b": set()})
        source_entries = agg_source_entries.get((a, pa, b, pb), {"a": set(), "b": set()})
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
                "source_entries_a": [
                    {"model": model, "chain": chain, "resnum": resnum}
                    for model, chain, resnum in sorted(source_entries["a"], key=lambda item: model_sort_key(item[0]) + (item[1], item[2]))
                ],
                "source_entries_b": [
                    {"model": model, "chain": chain, "resnum": resnum}
                    for model, chain, resnum in sorted(source_entries["b"], key=lambda item: model_sort_key(item[0]) + (item[1], item[2]))
                ],
                "models": sorted(
                    agg_models.get((a, pa, b, pb), set()),
                    key=model_sort_key,
                ),
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


def map_residue_colors_to_display(
    pos_info: Dict[str, List[Dict[str, str | int]]],
    residue_colors: Dict[str, Dict[int, str]],
) -> Dict[str, List[str | None]]:
    display_colors: Dict[str, List[str | None]] = {}
    for chain_id, infos in pos_info.items():
        colors: List[str | None] = []
        for info in infos:
            color = None
            pair = info.get("pair") if isinstance(info, dict) else None
            if pair:
                f_chain = str(pair.get("f_chain", ""))
                f_resnum = int(pair.get("f_resnum", 0))
                r_chain = str(pair.get("r_chain", ""))
                r_resnum = int(pair.get("r_resnum", 0))
                color = residue_colors.get(f_chain, {}).get(f_resnum)
                if color is None:
                    color = residue_colors.get(r_chain, {}).get(r_resnum)
            else:
                source_chain = str(info.get("source_chain", ""))
                resnum = int(info.get("resnum", 0))
                color = residue_colors.get(source_chain, {}).get(resnum)
            colors.append(color)
        display_colors[chain_id] = colors
    return display_colors


def build_html_metadata(
    invocation: str,
    models: List[Dict[str, str]],
) -> Dict[str, object]:
    return {
        "tool_version": TOOL_VERSION,
        "created_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "invocation": invocation,
        "models": models,
    }


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
    model_ids: List[str],
    display_residue_colors: Dict[str, List[str | None]],
    metadata: Dict[str, object],
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
                "chimerax_colors": display_residue_colors.get(chain, []),
            }
            for chain in chains
        ],
        "contacts": contacts,
        "max_count_atom": max_count_atom,
        "max_count_residue": max_count_residue,
        "default_order": default_order,
        "model_ids": model_ids,
        "has_chimerax_colors": any(
            any(color for color in display_residue_colors.get(chain, []))
            for chain in chains
        ),
        "info": metadata,
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
  box-sizing: border-box;
  padding: 8px 10px;
  border-radius: 8px;
  border: 1px solid #cbd5f5;
  font-size: 13px;
}
button {
  background: #1d4ed8;
  color: #fff;
  border: none;
  border-radius: 8px;
  padding: 8px 12px;
  cursor: pointer;
  font-size: 13px;
}
button.secondary {
  background: #3b82f6;
}
button:hover {
  background: #1e40af;
}
button.secondary:hover {
  background: #2563eb;
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
.mode-row {
  align-items: flex-start;
  justify-content: space-between;
}
.mode-options,
.display-options {
  display: flex;
  flex-direction: column;
  gap: 6px;
}
.display-header {
  font-size: 12px;
  color: #334155;
  font-weight: 600;
  line-height: 1.2;
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
.selection-menu .menu-group-label {
  font-size: 11px;
  color: #64748b;
  letter-spacing: 0.03em;
  text-transform: uppercase;
  margin-top: 8px;
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
.selection-menu .menu-row.wrap {
  flex-wrap: wrap;
}
.selection-menu .menu-button {
  flex: 1 1 0;
  min-width: 0;
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
.version-footer {
  margin-top: 8px;
  font-size: 11px;
  color: #64748b;
}
.button-stack {
  display: grid;
  gap: 10px;
}
.button-group-title {
  font-size: 11px;
  color: #64748b;
  letter-spacing: 0.03em;
  text-transform: uppercase;
}
.button-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 8px;
}
.button-grid button,
.button-stack > button {
  width: 100%;
}
.info-panel {
  position: absolute;
  right: 18px;
  top: 18px;
  width: min(560px, calc(100% - 36px));
  max-height: calc(100% - 36px);
  overflow: auto;
  display: none;
  z-index: 11;
  background: rgba(255, 255, 255, 0.98);
  border: 1px solid #cbd5e1;
  border-radius: 12px;
  box-shadow: 0 18px 40px rgba(15, 23, 42, 0.22);
  padding: 14px;
}
.info-panel.visible {
  display: block;
}
.info-panel-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 12px;
  margin-bottom: 10px;
}
.info-panel-title {
  font-size: 14px;
  font-weight: 700;
  color: #0f172a;
}
.info-panel-section {
  margin-top: 10px;
}
.info-panel-label {
  font-size: 11px;
  color: #64748b;
  text-transform: uppercase;
  letter-spacing: 0.03em;
  margin-bottom: 4px;
}
.info-panel pre,
.info-panel .info-box {
  margin: 0;
  padding: 10px;
  background: #f8fafc;
  border: 1px solid #e2e8f0;
  border-radius: 8px;
  white-space: pre-wrap;
  word-break: break-word;
  font-size: 12px;
  color: #0f172a;
}
.info-model {
  margin-top: 8px;
}
.info-model:first-child {
  margin-top: 0;
}
.copy-toast {
  position: absolute;
  left: 12px;
  top: 12px;
  z-index: 8;
  max-width: min(360px, calc(100% - 24px));
  padding: 8px 10px;
  border-radius: 8px;
  border: 1px solid #bfdbfe;
  background: #eff6ff;
  color: #1e3a8a;
  font-size: 12px;
  box-shadow: 0 10px 24px rgba(15, 23, 42, 0.14);
  opacity: 0;
  pointer-events: none;
  transform: translateY(-4px);
  transition: opacity 180ms ease, transform 180ms ease;
}
.copy-toast.visible {
  opacity: 1;
  transform: translateY(0);
}
.copy-toast.error {
  border-color: #fecaca;
  background: #fef2f2;
  color: #991b1b;
}
.copy-toast pre {
  margin: 6px 0 0 0;
  white-space: pre-wrap;
  word-break: break-word;
  font-family: "SFMono-Regular", Consolas, monospace;
  font-size: 11px;
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
    <div class="control">
      <label for="plotTitle">Plot title</label>
      <input id="plotTitle" type="text" value="__TITLE__">
    </div>
    <div class="control">
      <label for="threshold">Minimum contact visibility threshold</label>
      <div class="row">
        <input id="threshold" type="range" min="1" max="__MAX_COUNT__" step="1" value="2">
        <input id="thresholdInput" type="number" min="1" max="__MAX_COUNT__" step="1" value="2">
      </div>
    </div>
    <div class="control">
      <label>Contact mode</label>
      <div class="row mode-row">
        <div class="mode-options">
          <label class="chain-toggle"><input type="radio" name="countMode" value="atom" checked>Atom</label>
          <label class="chain-toggle"><input type="radio" name="countMode" value="residue">Residue</label>
        </div>
        <div class="display-options">
          <div class="display-header">Display</div>
          <label class="chain-toggle"><input id="useTransparency" type="checkbox" checked>Transparency</label>
          <label class="chain-toggle"><input id="useDefaultColors" type="checkbox" checked>Default colors</label>
        </div>
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
    </div>
    <div class="control">
      <div class="button-stack">
        <div>
          <div class="button-group-title">Export</div>
          <div class="button-grid" style="margin-top:6px;">
            <button id="downloadSvg">Download SVG</button>
            <button id="saveHtml" class="secondary">Save HTML</button>
          </div>
          <button id="exportCxc" class="secondary" style="margin-top:8px;">ChimeraX Colors</button>
        </div>
        <div>
          <div class="button-group-title">Session</div>
          <div class="button-grid" style="margin-top:6px;">
            <button id="saveSession" class="secondary">Save Session</button>
            <button id="loadSession" class="secondary">Load Session</button>
          </div>
        </div>
        <button id="infoButton" class="secondary">Info</button>
        <input id="loadSessionFile" type="file" accept=".json,application/json" style="display:none;">
      </div>
      <div class="status" id="renderStatus">Render status: pending</div>
      <div class="version-footer">CircosContacts v__VERSION__</div>
    </div>
  </div>
  <div class="plot-wrap panel" id="plotWrap" style="position:relative;">
    <svg id="circos" viewBox="0 0 1400 900" role="img" aria-label="Circos plot"></svg>
    <div class="copy-toast" id="copyToast"></div>
    <div class="tooltip" id="tooltip"></div>
    <div class="selection-menu" id="selectionMenu"></div>
    <div class="info-panel" id="infoPanel"></div>
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
const orderInput = document.getElementById('order');
const orderHint = document.getElementById('orderHint');
const renderStatus = document.getElementById('renderStatus');
const copyToast = document.getElementById('copyToast');
const plotWrap = document.getElementById('plotWrap');
const saveHtml = document.getElementById('saveHtml');
const exportCxc = document.getElementById('exportCxc');
const selectionMenu = document.getElementById('selectionMenu');
const saveSessionBtn = document.getElementById('saveSession');
const loadSessionBtn = document.getElementById('loadSession');
const loadSessionFile = document.getElementById('loadSessionFile');
const infoButton = document.getElementById('infoButton');
const infoPanel = document.getElementById('infoPanel');
const plotTitleInput = document.getElementById('plotTitle');
const useTransparencyInput = document.getElementById('useTransparency');
const useDefaultColorsInput = document.getElementById('useDefaultColors');

const chainMap = new Map(DATA.chains.map(c => [c.id, c]));
const defaultOrder = DATA.default_order.slice();
let currentOrder = defaultOrder.slice();
const contactVisibility = new Map(DATA.chains.map(c => [c.id, true]));
const chainFlip = new Map(DATA.chains.map(c => [c.id, false]));
const chainDisplayName = new Map(DATA.chains.map(c => [c.id, c.id]));
let lockedBottomChain = null;
let countMode = 'atom';
let useTransparency = true;
let useDefaultColors = true;

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
let pendingNameFocusChain = null;
const dnaSeqSideCache = new Map();
const CANVAS_W = 1400;
const CANVAS_H = 900;

const AA1 = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G',
  HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
  THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V', SEC: 'U', PYL: 'O'
};
const hasChimeraXColors = DATA.has_chimerax_colors === true;
useDefaultColorsInput.disabled = !hasChimeraXColors;
useDefaultColorsInput.checked = true;

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

function linkCount(link) {
  return countMode === 'atom' ? link.count_atom : link.count_residue;
}

function countMeetsThreshold(count, thresholdValue) {
  return count >= thresholdValue;
}

function thresholdIntensity(count, thresholdValue, maxCount) {
  return clamp((count - thresholdValue + 1) / Math.max(1, maxCount - thresholdValue + 1), 0, 1);
}

function escapeHtml(value) {
  return String(value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/\"/g, '&quot;');
}

function showCopyToast(message, command, isError = false) {
  if (!copyToast) return;
  const safeMessage = escapeHtml(message || 'Copied ChimeraX command');
  const safeCommand = escapeHtml(command || '');
  copyToast.innerHTML = `<div>${safeMessage}</div>${safeCommand ? `<pre style="margin:6px 0 0;white-space:pre-wrap;font:11px IBM Plex Mono, monospace;">${safeCommand}</pre>` : ''}`;
  copyToast.style.borderColor = isError ? '#dc2626' : '#1d4ed8';
  copyToast.classList.add('visible');
  window.clearTimeout(showCopyToast._timer);
  showCopyToast._timer = window.setTimeout(() => {
    copyToast.classList.remove('visible');
  }, 3200);
}

async function copyTextToClipboard(text) {
  const value = String(text || '');
  if (navigator.clipboard && navigator.clipboard.writeText) {
    try {
      await navigator.clipboard.writeText(value);
      return true;
    } catch (err) {
      // fall through
    }
  }
  const ta = document.createElement('textarea');
  ta.value = value;
  ta.setAttribute('readonly', 'readonly');
  ta.style.position = 'fixed';
  ta.style.opacity = '0';
  ta.style.left = '-9999px';
  document.body.appendChild(ta);
  ta.select();
  let ok = false;
  try {
    ok = document.execCommand('copy');
  } catch (err) {
    ok = false;
  }
  document.body.removeChild(ta);
  return ok;
}

function compareModelId(a, b) {
  const sa = String(a || '');
  const sb = String(b || '');
  const ma = /^#(\\d+)$/.exec(sa);
  const mb = /^#(\\d+)$/.exec(sb);
  if (ma && mb) return Number(ma[1]) - Number(mb[1]);
  if (ma) return -1;
  if (mb) return 1;
  return sa.localeCompare(sb);
}

function compressNumbers(values) {
  const nums = Array.from(new Set(values.map((value) => Number(value)).filter(Number.isFinite))).sort((a, b) => a - b);
  const out = [];
  for (let i = 0; i < nums.length; i += 1) {
    const start = nums[i];
    let end = start;
    while (i + 1 < nums.length && nums[i + 1] === end + 1) {
      i += 1;
      end = nums[i];
    }
    out.push(start === end ? `${start}` : `${start}-${end}`);
  }
  return out.join(',');
}

function buildSelectSpec(records) {
  const normalized = [];
  const seen = new Set();
  for (const rec of records || []) {
    if (!rec || !rec.chain || !Number.isFinite(Number(rec.resnum))) continue;
    const model = String(rec.model || '');
    const chain = String(rec.chain);
    const resnum = Number(rec.resnum);
    const key = `${model}|${chain}|${resnum}`;
    if (seen.has(key)) continue;
    seen.add(key);
    normalized.push({ model, chain, resnum });
  }
  if (!normalized.length) return '';

  const groups = new Map();
  for (const rec of normalized) {
    const groupKey = rec.model || '__nomodel__';
    if (!groups.has(groupKey)) groups.set(groupKey, new Map());
    const chainMapForModel = groups.get(groupKey);
    if (!chainMapForModel.has(rec.chain)) chainMapForModel.set(rec.chain, []);
    chainMapForModel.get(rec.chain).push(rec.resnum);
  }

  const modelKeys = Array.from(groups.keys()).sort((a, b) => compareModelId(a === '__nomodel__' ? '' : a, b === '__nomodel__' ? '' : b));
  const signatureBuckets = new Map();
  for (const modelKey of modelKeys) {
    const chainEntries = Array.from(groups.get(modelKey).entries())
      .sort((a, b) => a[0].localeCompare(b[0]))
      .map(([chain, nums]) => `${chain}:${compressNumbers(nums)}`);
    const signature = chainEntries.join('/');
    if (!signatureBuckets.has(signature)) signatureBuckets.set(signature, []);
    signatureBuckets.get(signature).push(modelKey === '__nomodel__' ? '' : modelKey);
  }

  const parts = [];
  for (const [signature, models] of signatureBuckets.entries()) {
    const chainSpec = '/' + signature;
    const compactModels = models.filter(Boolean);
    if (compactModels.length) {
      const modelSpec = '#' + compactModels.map((m) => m.replace(/^#/, '')).sort((a, b) => Number(a) - Number(b)).join(',');
      parts.push(`${modelSpec}${chainSpec}`);
    } else {
      parts.push(chainSpec);
    }
  }
  if (!parts.length) return '';
  if (parts.length === 1) return `select ${parts[0]}`;
  return ['select clear', ...parts.map((part) => `select add ${part}`)].join('\\n');
}

function mergedModelRecords(sourceEntries, fallbackSources, models) {
  const records = [];
  if (Array.isArray(sourceEntries) && sourceEntries.length) {
    for (const rec of sourceEntries) {
      records.push({ model: rec.model || '', chain: rec.chain, resnum: Number(rec.resnum) });
    }
    return records;
  }
  const modelList = Array.isArray(models) && models.length ? models : [''];
  for (const src of fallbackSources || []) {
    for (const model of modelList) {
      records.push({ model, chain: src.chain, resnum: Number(src.resnum) });
    }
  }
  return records;
}

function linkSelectionRecords(link) {
  const records = [];
  records.push(...mergedModelRecords(link.source_entries_a, link.sources_a, link.models));
  records.push(...mergedModelRecords(link.source_entries_b, link.sources_b, link.models));
  return records;
}

function regionSelectionRecords(sel) {
  const chainMeta = chainMap.get(sel.chainId);
  if (!chainMeta) return [];
  const [startPos, endPos] = normalizeRange(sel.startPos, sel.endPos);
  const models = Array.isArray(DATA.model_ids) && DATA.model_ids.length ? DATA.model_ids : [''];
  const records = [];
  for (let pos = startPos; pos <= endPos; pos += 1) {
    const info = chainMeta.resinfo[pos - 1];
    if (!info) continue;
    if (info.pair) {
      for (const model of models) {
        records.push({ model, chain: info.pair.f_chain, resnum: Number(info.pair.f_resnum) });
        records.push({ model, chain: info.pair.r_chain, resnum: Number(info.pair.r_resnum) });
      }
    } else {
      for (const model of models) {
        records.push({ model, chain: info.source_chain, resnum: Number(info.resnum) });
      }
    }
  }
  return records;
}

function regionContactSelectionRecords(sel, thresholdValue, visibleChains) {
  const records = [];
  for (const link of DATA.contacts) {
    if (!staticLinkVisible(link, thresholdValue, visibleChains)) continue;
    if (!linkMatchesSelection(link, sel)) continue;
    records.push(...linkSelectionRecords(link));
  }
  return records;
}

async function copyChimeraXSelectCommand(records, message) {
  const command = buildSelectSpec(records);
  if (!command) {
    showCopyToast('No ChimeraX selection could be generated', '', true);
    return;
  }
  const ok = await copyTextToClipboard(command);
  if (ok) showCopyToast(message || 'Copied ChimeraX select command', command, false);
  else showCopyToast('Copy failed — command shown below', command, true);
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
        event.preventDefault();
        const inputs = Array.from(orderHint.querySelectorAll('input[data-role="name"]:not([disabled])'));
        const idx = inputs.indexOf(target);
        const nextInput = inputs.length ? inputs[(idx + 1) % inputs.length] : null;
        pendingNameFocusChain = nextInput ? nextInput.getAttribute('data-chain') : null;
        if (commitName(target)) {
          safeRender();
        } else {
          target.value = chainDisplayName.get(target.getAttribute('data-chain')) || target.getAttribute('data-chain');
          if (nextInput) {
            pendingNameFocusChain = null;
            nextInput.focus();
            nextInput.select();
          }
        }
      }
    });

    node.addEventListener('blur', (event) => {
      const target = event.currentTarget;
      const chainId = target.getAttribute('data-chain');
      target.value = chainDisplayName.get(chainId) || chainId;
    });
  });
  if (pendingNameFocusChain) {
    const target = orderHint.querySelector(`input[data-role="name"][data-chain="${pendingNameFocusChain}"]`);
    pendingNameFocusChain = null;
    if (target && !target.disabled) {
      target.focus();
      target.select();
    }
  }
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

function selectionSequenceHtml(sel, thresholdValue, visibleChains) {
  const chainMeta = chainMap.get(sel.chainId);
  if (!chainMeta) return '';
  const activePositions = new Set();
  for (const link of DATA.contacts) {
    if (!staticLinkVisible(link, thresholdValue, visibleChains)) continue;
    if (link.a === sel.chainId) activePositions.add(link.a_pos);
    if (link.b === sel.chainId) activePositions.add(link.b_pos);
  }
  const dnaSide = dnaPreferredSide(sel.chainId, chainMeta);
  const [s, e] = normalizeRange(sel.startPos, sel.endPos);
  const tokens = [];
  for (let pos = s; pos <= e; pos += 1) {
    const info = chainMeta.resinfo[pos - 1];
    if (!info) continue;
    const token = info.pair
      ? (dnaSide === 'r' ? baseCode(info.pair.r_resname) : baseCode(info.pair.f_resname))
      : (AA1[String(info.resname || '').toUpperCase()] || 'X');
    const safe = escapeHtml(token);
    tokens.push(activePositions.has(pos) ? `<strong style="color:#0000FF;">${safe}</strong>` : safe);
  }
  return tokens.join('');
}

function staticLinkVisible(link, thresholdValue, visibleChains) {
  const count = linkCount(link);
  if (!countMeetsThreshold(count, thresholdValue)) return false;
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

function escapeAttr(value) {
  return escapeHtml(value).replace(/'/g, '&#39;');
}

function closeInfoPanel() {
  infoPanel.classList.remove('visible');
  infoPanel.innerHTML = '';
}

function openInfoPanel() {
  const info = DATA.info || {};
  const models = Array.isArray(info.models) ? info.models : [];
  const modelBlocks = models.length
    ? models.map((model, idx) => {
        const parts = [];
        if (model.name) parts.push(`<div><strong>Name:</strong> ${escapeHtml(model.name)}</div>`);
        if (model.filename) parts.push(`<div><strong>File:</strong> ${escapeHtml(model.filename)}</div>`);
        if (model.path) parts.push(`<div><strong>Path:</strong> ${escapeHtml(model.path)}</div>`);
        if (model.atomspec) parts.push(`<div><strong>Spec:</strong> ${escapeHtml(model.atomspec)}</div>`);
        return `<div class="info-box info-model">${parts.join('')}</div>`;
      }).join('')
    : '<div class="info-box">No model metadata stored.</div>';
  infoPanel.innerHTML = `
    <div class="info-panel-header">
      <div class="info-panel-title">Plot Information</div>
      <button id="closeInfoPanel" class="secondary">Close</button>
    </div>
    <div class="info-panel-section">
      <div class="info-panel-label">Version</div>
      <div class="info-box">${escapeHtml(info.tool_version || '')}</div>
    </div>
    <div class="info-panel-section">
      <div class="info-panel-label">Created</div>
      <div class="info-box">${escapeHtml(info.created_at || '')}</div>
    </div>
    <div class="info-panel-section">
      <div class="info-panel-label">Invocation</div>
      <pre>${escapeHtml(info.invocation || '')}</pre>
    </div>
    <div class="info-panel-section">
      <div class="info-panel-label">Models</div>
      ${modelBlocks}
    </div>
  `;
  infoPanel.classList.add('visible');
  const closeBtn = document.getElementById('closeInfoPanel');
  if (closeBtn) {
    closeBtn.addEventListener('click', (event) => {
      event.stopPropagation();
      closeInfoPanel();
    });
  }
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

function parseInlineStyle(styleText) {
  const style = {};
  for (const chunk of String(styleText || '').split(';')) {
    const idx = chunk.indexOf(':');
    if (idx === -1) continue;
    const key = chunk.slice(0, idx).trim().toLowerCase();
    const value = chunk.slice(idx + 1).trim();
    if (key) style[key] = value;
  }
  return style;
}

function extractStyledTextSegments(root, inherited = null) {
  const base = inherited || { fill: '#0f172a', weight: '400', fontStyle: 'normal' };
  const out = [];
  const walk = (node, current) => {
    if (!node) return;
    if (node.nodeType === Node.TEXT_NODE) {
      if (node.nodeValue) {
        out.push({ text: node.nodeValue, fill: current.fill, weight: current.weight, fontStyle: current.fontStyle });
      }
      return;
    }
    if (node.nodeType !== Node.ELEMENT_NODE) return;
    const tag = node.tagName.toLowerCase();
    if (tag === 'br') {
      out.push({ text: '\\n', fill: current.fill, weight: current.weight, fontStyle: current.fontStyle });
      return;
    }
    const next = { ...current };
    if (tag === 'strong' || tag === 'b') next.weight = '700';
    if (tag === 'em' || tag === 'i') next.fontStyle = 'italic';
    const inline = parseInlineStyle(node.getAttribute('style'));
    if (inline.color) next.fill = inline.color;
    if (inline['font-weight']) next.weight = inline['font-weight'];
    if (inline['font-style']) next.fontStyle = inline['font-style'];
    Array.from(node.childNodes).forEach((child) => walk(child, next));
  };
  walk(root, base);
  return out;
}

function wrapStyledSegments(segments, maxWidthPx) {
  const charPx = 7.2;
  const maxChars = Math.max(1, Math.floor(maxWidthPx / charPx));
  const lines = [[]];
  let lineLen = 0;
  const pushChar = (style, ch) => {
    const line = lines[lines.length - 1];
    const prev = line[line.length - 1];
    if (prev && prev.fill === style.fill && prev.weight === style.weight && prev.fontStyle === style.fontStyle) {
      prev.text += ch;
    } else {
      line.push({ text: ch, fill: style.fill, weight: style.weight, fontStyle: style.fontStyle });
    }
    lineLen += 1;
  };
  for (const seg of segments) {
    const style = { fill: seg.fill || '#0f172a', weight: seg.weight || '400', fontStyle: seg.fontStyle || 'normal' };
    for (const ch of String(seg.text || '')) {
      if (ch === '\\n') {
        lines.push([]);
        lineLen = 0;
        continue;
      }
      if (lineLen >= maxChars) {
        lines.push([]);
        lineLen = 0;
      }
      pushChar(style, ch);
    }
  }
  return lines;
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
  threshold.value = clamp(Number(threshold.value), 1, maxCount);
  thresholdInput.value = threshold.value;
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

  const linkVisible = (link) => {
    if (!countMeetsThreshold(linkCount(link), thresholdValue)) return false;
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

  const chainPathData = (chain, chainMeta) => {
    if (chain.length >= 10) {
      const indicatorResidues = Math.max(2, Math.min(6, Math.round(chain.length * 0.03)));
      const segmentAngle = (indicatorResidues / chain.length) * (chain.end - chain.start);
      const startAtStart = ((chainMeta.start_pos || 1) === 1) !== (chain.flip === true);
      const flare = 6;
      return styledArcPath(cx, cy, outerR, innerR, chain.start, chain.end, startAtStart, flare, segmentAngle);
    }
    return arcPath(cx, cy, outerR, innerR, chain.start, chain.end);
  };

  const appendColoredChainSegments = (chainId, chain, chainMeta, opacityValue, ranges = null) => {
    const clipId = `chain-clip-${chainId.replace(/[^A-Za-z0-9_-]/g, '_')}-${Math.round(chain.start * 1000)}-${String(opacityValue).replace(/[^A-Za-z0-9_-]/g, '_')}`;
    const clipPath = document.createElementNS('http://www.w3.org/2000/svg', 'clipPath');
    clipPath.setAttribute('id', clipId);
    const clipShape = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    clipShape.setAttribute('d', chainPathData(chain, chainMeta));
    clipPath.appendChild(clipShape);
    labelDefs.appendChild(clipPath);

    const segmentGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    segmentGroup.setAttribute('clip-path', `url(#${clipId})`);
    segmentGroup.setAttribute('opacity', String(opacityValue));
    const allowed = ranges ? new Set(ranges.flatMap(([a, b]) => Array.from({ length: b - a + 1 }, (_, idx) => a + idx))) : null;
    const span = chain.end - chain.start;
    const colors = Array.isArray(chainMeta.chimerax_colors) ? chainMeta.chimerax_colors : [];
    for (let displayIdx = 1; displayIdx <= chain.length; displayIdx += 1) {
      if (allowed && !allowed.has(displayIdx)) continue;
      const sourceIdx = chain.flip ? (chain.length - displayIdx + 1) : displayIdx;
      const color = colors[sourceIdx - 1] || chain.color;
      const segStart = chain.start + ((displayIdx - 1) / chain.length) * span;
      const segEnd = chain.start + (displayIdx / chain.length) * span;
      const segPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      segPath.setAttribute('d', arcPath(cx, cy, outerR, innerR, segStart, segEnd));
      segPath.setAttribute('fill', color);
      segmentGroup.appendChild(segPath);
    }
    arcsGroup.appendChild(segmentGroup);
  };

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
    const hasChimeraXChainColors =
      !useDefaultColors &&
      Array.isArray(chainMeta.chimerax_colors) &&
      chainMeta.chimerax_colors.some((color) => color);

    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    const baseOpacity = (hasActiveMask && !fullyActive) || (anyActiveMask && !hasActiveMask) ? '0.45' : '0.9';
    path.setAttribute('d', chainPathData(chain, chainMeta));
    if (hasChimeraXChainColors) {
      appendColoredChainSegments(id, chain, chainMeta, baseOpacity);
      if (hasActiveMask && !fullyActive) {
        appendColoredChainSegments(id, chain, chainMeta, 0.9, activeRanges);
      }
      path.setAttribute('fill', '#ffffff');
      path.setAttribute('opacity', '0.001');
    } else {
      path.setAttribute('fill', chain.color);
      path.setAttribute('opacity', baseOpacity);
    }
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
    if (hasActiveMask && !fullyActive && !hasChimeraXChainColors) {
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
    label.setAttribute('font-size', '24');
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
      const contactSelectEnabled = regionContactSelectionRecords(item, thresholdValue, visibleChains).length > 0;
      const maxLen = chainMap.get(item.chainId).length;
      return `
        <div class="menu-block">
          <div class="menu-title">${item.chainId}:${a}-${b}</div>
          <div class="menu-group-label">Annotations</div>
          <div class="menu-row wrap">
            <button class="menu-button" data-action="callout" data-id="${item.id}" ${seqEnabled ? '' : 'disabled'}>Sequence callout</button>
            <button class="menu-button" data-action="comment" data-id="${item.id}">Comment</button>
          </div>
          <div class="menu-group-label">ChimeraX Select</div>
          <div class="menu-row wrap">
            <button class="menu-button" data-action="chimeraSelectAll" data-id="${item.id}">All models</button>
            <button class="menu-button" data-action="chimeraSelectContacts" data-id="${item.id}" ${contactSelectEnabled ? '' : 'disabled'}>Contacts @ threshold</button>
          </div>
          <div class="menu-group-label">Range</div>
          <div class="menu-row wrap">
            <button class="menu-button" data-action="autotrim" data-id="${item.id}" ${canTrim ? '' : 'disabled'}>Autotrim</button>
            <button class="menu-button" data-action="clear" data-id="${item.id}">Clear</button>
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
          <div class="menu-group-label">Link Display</div>
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
        } else if (act === 'chimeraSelectAll') {
          copyChimeraXSelectCommand(regionSelectionRecords(selObj), `Copied ChimeraX select command for ${selObj.chainId}:${selObj.startPos}-${selObj.endPos}`);
        } else if (act === 'chimeraSelectContacts') {
          copyChimeraXSelectCommand(
            regionContactSelectionRecords(selObj, thresholdValue, visibleChains),
            `Copied ChimeraX contact select command for ${selObj.chainId}:${selObj.startPos}-${selObj.endPos}`
          );
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
        if (act !== 'autotrim' && act !== 'applyEdit' && act !== 'callout') {
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
      const htmlVal = callout.kind === 'sequence'
        ? selectionSequenceHtml(sel, thresholdValue, visibleChains)
        : (callout.html && callout.html.length
          ? callout.html
          : escapeHtml(textVal).replace(/\\n/g, '<br>'));
      const titleText = callout.kind === 'sequence'
        ? (callout.title || calloutRangeTitle(sel))
        : '';
      const titleBarH = 16;
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
      fo.setAttribute('data-callout-bar-height', `${titleBarH}`);
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

    const count = linkCount(link);
    const t = thresholdIntensity(count, thresholdValue, maxCount);
    const eased = Math.pow(t, 0.65);
    const alpha = useTransparency ? (minAlpha + (maxAlpha - minAlpha) * eased) : 1;
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
      const formatModels = (models) => {
        if (!Array.isArray(models) || !models.length) return '';
        const flat = models.flatMap((item) => String(item).split(',')).filter(Boolean);
        const allNumeric = flat.every((item) => /^#\\d+$/.test(item));
        if (allNumeric) {
          return `#${flat.map((item) => item.slice(1)).join(',')}`;
        }
        return flat.join(', ');
      };
      const modelText = formatModels(link.models);
      tooltip.innerHTML = `
        <div><strong>${formatInfo(aInfo)}</strong> ↔ <strong>${formatInfo(bInfo)}</strong></div>
        <div>Count (${countMode}): ${count}</div>
        ${modelText ? `<div>Models: ${modelText}</div>` : ''}
      `;
      tooltip.style.opacity = '1';
      tooltip.classList.remove('clickable');
      positionTooltip(event);
    });
    path.style.cursor = 'copy';
    path.addEventListener('click', (event) => {
      event.stopPropagation();
      copyChimeraXSelectCommand(
        linkSelectionRecords(link),
        `Copied ChimeraX select command for ${link.a}:${link.a_pos} ↔ ${link.b}:${link.b_pos}`
      );
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
      countMode,
      useTransparency,
      useDefaultColors,
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
  if (controls.countMode === 'atom' || controls.countMode === 'residue') {
    countMode = controls.countMode;
  }
  useTransparency = controls.useTransparency !== false;
  useTransparencyInput.checked = useTransparency;
  useDefaultColors = hasChimeraXColors ? (controls.useDefaultColors !== false) : true;
  useDefaultColorsInput.checked = useDefaultColors;
  useDefaultColorsInput.disabled = !hasChimeraXColors;
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
useTransparencyInput.addEventListener('change', () => {
  useTransparency = useTransparencyInput.checked;
  safeRender();
});
useDefaultColorsInput.addEventListener('change', () => {
  useDefaultColors = useDefaultColorsInput.checked || !hasChimeraXColors;
  safeRender();
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
  if (infoPanel.classList.contains('visible') && !infoPanel.contains(event.target) && event.target !== infoButton) {
    closeInfoPanel();
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

infoButton.addEventListener('click', (event) => {
  event.stopPropagation();
  if (infoPanel.classList.contains('visible')) closeInfoPanel();
  else openInfoPanel();
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
    const div = fo.querySelector('div');
    const styledSegments = div
      ? extractStyledTextSegments(div)
      : [{ text: textVal, fill: '#0f172a', weight: '400', fontStyle: 'normal' }];
    const lines = wrapStyledSegments(styledSegments, Math.max(10, w - 4));
    const boxId = fo.getAttribute('data-callout-fo');
    const box = boxId ? exportSvg.querySelector(`rect[data-callout-box="${boxId}"]`) : null;
    const title = boxId ? exportSvg.querySelector(`text[data-callout-title="${boxId}"]`) : null;
    const titleBarH = Number(fo.getAttribute('data-callout-bar-height') || (title ? 16 : 4));
    const neededHeight = Math.max(40, titleBarH + 4 + Math.max(1, lines.length) * 14);
    if (box) {
      box.setAttribute('height', `${neededHeight}`);
    }
    const text = document.createElementNS(ns, 'text');
    text.setAttribute('font-size', '12');
    text.setAttribute('font-family', 'IBM Plex Sans, sans-serif');
    lines.forEach((line, idx) => {
      if (!line.length) {
        const t = document.createElementNS(ns, 'tspan');
        t.setAttribute('x', `${x + 2}`);
        t.setAttribute('y', `${y + 12 + idx * 14}`);
        t.textContent = ' ';
        text.appendChild(t);
        return;
      }
      line.forEach((seg, segIdx) => {
        const t = document.createElementNS(ns, 'tspan');
        if (segIdx === 0) {
          t.setAttribute('x', `${x + 2}`);
          t.setAttribute('y', `${y + 12 + idx * 14}`);
        }
        t.setAttribute('fill', seg.fill || '#0f172a');
        if (seg.weight && seg.weight !== '400') t.setAttribute('font-weight', seg.weight);
        if (seg.fontStyle && seg.fontStyle !== 'normal') t.setAttribute('font-style', seg.fontStyle);
        t.textContent = seg.text;
        text.appendChild(t);
      });
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
    const count = linkCount(link);
    if (!countMeetsThreshold(count, thresholdValue)) continue;
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
        .replace("__VERSION__", TOOL_VERSION)
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

    (
        raw_counts_atom,
        residue_pairs_by_file,
        contact_models,
        contact_sources,
    ) = parse_contacts(contact_paths)
    contacts, max_count_atom, max_count_residue = aggregate_contacts(
        raw_counts_atom,
        residue_pairs_by_file,
        contact_models,
        contact_sources,
        pos_map,
        display_chain_of,
    )

    chain_lengths = {
        chain: len(pos_info.get(chain, [])) for chain in display_chains
    }
    chain_colors = build_colors(display_chains)
    default_order = [chain for chain in display_chains if chain_lengths.get(chain, 0) > 1] or display_chains[:]
    model_ids = sorted(
        {
            model
            for models in contact_models.values()
            for model in models
            if model
        },
        key=lambda label: (0, int(label[1:])) if label.startswith("#") and label[1:].isdigit() else (1, label),
    )
    metadata = build_html_metadata(
        " ".join(shlex.quote(arg) for arg in sys.argv),
        [
            {
                "name": path.name,
                "filename": path.name,
                "path": str(path.resolve()),
            }
            for path in cif_paths
        ],
    )

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
        model_ids,
        {},
        metadata,
        output_path,
    )
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
