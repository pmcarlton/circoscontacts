from __future__ import annotations

import tempfile
import re
import webbrowser
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

from chimerax.core.commands import BoolArg, CmdDesc, IntArg, ObjectsArg, StringArg, register
from chimerax.core.errors import UserError

from . import contacts_circos as circos


def _atoms_from_arg(arg):
    if arg is None:
        return None
    atoms = getattr(arg, "atoms", None)
    if atoms is not None:
        return atoms
    return arg


def _structures_from_objects(session, objects):
    from chimerax.atomic import AtomicStructure

    if objects is None:
        models = session.models.list(type=AtomicStructure)
        return [m for m in models if hasattr(m, "atomspec")]

    atoms = _atoms_from_arg(objects)
    if atoms is None or len(atoms) == 0:
        return []

    seen = set()
    structures = []
    for atom in atoms:
        structure = atom.structure
        sid = id(structure)
        if sid in seen:
            continue
        seen.add(sid)
        structures.append(structure)
    return structures


def _chain_residue_maps(structures: Iterable) -> Tuple[Dict[str, Dict[int, str]], Dict[str, List[int]]]:
    chain_res: Dict[str, Dict[int, str]] = {}
    chain_resnums: Dict[str, set[int]] = {}
    for structure in structures:
        for residue in structure.residues:
            chain = residue.chain_id
            try:
                resnum = int(residue.number)
            except Exception:
                continue
            resname = residue.name
            if not chain:
                continue
            chain_res.setdefault(chain, {})
            chain_resnums.setdefault(chain, set())
            prev = chain_res[chain].get(resnum)
            if prev is None:
                chain_res[chain][resnum] = resname
            chain_resnums[chain].add(resnum)
    chain_resnums_sorted = {k: sorted(v) for k, v in chain_resnums.items()}
    return chain_res, chain_resnums_sorted


def _objects_to_chain_resnums(objects) -> Dict[str, Set[int]]:
    chain_resnums: Dict[str, Set[int]] = {}
    if objects is None:
        return chain_resnums
    atoms = _atoms_from_arg(objects)
    if atoms is None:
        return chain_resnums
    for atom in atoms:
        residue = atom.residue
        chain = residue.chain_id
        try:
            resnum = int(residue.number)
        except Exception:
            continue
        chain_resnums.setdefault(chain, set()).add(resnum)
    return chain_resnums


def _to_ranges(values: Iterable[int]) -> List[Tuple[int, int]]:
    nums = sorted(set(values))
    if not nums:
        return []
    ranges: List[Tuple[int, int]] = []
    start = prev = nums[0]
    for value in nums[1:]:
        if value == prev + 1:
            prev = value
            continue
        ranges.append((start, prev))
        start = prev = value
    ranges.append((start, prev))
    return ranges


def _active_ranges_from_maps(
    active_chain_resnums: Dict[str, Set[int]],
    pos_map: Dict[str, Dict[int, int]],
    display_chain_of: Dict[str, str],
) -> Dict[str, List[Tuple[int, int]]]:
    active_display: Dict[str, Set[int]] = {}
    for source_chain, source_resnums in active_chain_resnums.items():
        mapped = pos_map.get(source_chain, {})
        display_chain = display_chain_of.get(source_chain, source_chain)
        for resnum in source_resnums:
            pos = mapped.get(resnum)
            if pos is None:
                continue
            active_display.setdefault(display_chain, set()).add(pos)
    return {chain: _to_ranges(pos_set) for chain, pos_set in active_display.items()}


def _spec_chain_hints(spec: str | None) -> Tuple[Set[str], Dict[str, Set[int]]]:
    chains: Set[str] = set()
    chain_resnums: Dict[str, Set[int]] = {}
    if not spec:
        return chains, chain_resnums

    for match in re.finditer(r"/([A-Za-z0-9,]+)(?::([0-9]+)(?:-([0-9]+))?)?", spec):
        raw_chains = match.group(1)
        start = match.group(2)
        end = match.group(3)
        chain_ids = [c.strip() for c in raw_chains.split(",") if c.strip()]
        for chain_id in chain_ids:
            chains.add(chain_id)
            if start is None:
                continue
            a = int(start)
            b = int(end) if end is not None else a
            if b < a:
                a, b = b, a
            chain_resnums.setdefault(chain_id, set()).update(range(a, b + 1))
    return chains, chain_resnums


def _run_contacts_for_structures(
    session,
    structures: Iterable,
    work_dir: Path,
    base_spec: str | None = None,
    restrict_spec: str | None = None,
    inter_model: bool = False,
    intra_mol: bool = False,
) -> List[Path]:
    from chimerax.core.commands import run

    contact_paths: List[Path] = []

    for idx, structure in enumerate(structures, start=1):
        spec = structure.atomspec
        if base_spec:
            contact_spec = f"({base_spec}) & ({spec})"
        else:
            contact_spec = f"{spec} & protein"
        command = (
            f'contacts {contact_spec} interModel {str(inter_model).lower()} '
            f'intramol {str(intra_mol).lower()}'
        )
        if restrict_spec:
            command += f' restrict (({restrict_spec}) & ({spec}))'

        contacts_path = work_dir / f"model_{idx:03d}.contacts"
        command += f' saveFile "{contacts_path}"'
        run(session, command, log=False)

        contact_paths.append(contacts_path)

    return contact_paths


def _open_html(session, html_path: Path):
    # Always open in system browser for reliable visibility in GUI sessions.
    # (Embedded HtmlToolInstance exists in some versions but is not consistently visible.)
    opened = False
    try:
        opened = bool(webbrowser.open(html_path.as_uri()))
    except Exception as err:
        session.logger.error(f"circoscontacts: failed to request browser open: {err}")
        return
    if opened:
        session.logger.info(f"circoscontacts: opened HTML in default browser ({html_path})")
    else:
        session.logger.warning(
            f"circoscontacts: browser open request did not confirm success ({html_path})"
        )


def circoscontacts(
    session,
    objects=None,
    restrict=None,
    output_dir=None,
    title="ChimeraX Contact Circos",
    inter_model=False,
    intramol=False,
    dna_mismatches=0,
    dna_mode="auto",
    open_html=True,
    keep_temp=False,
):
    """
    Generate a circos contact plot from currently open (or selected) models.
    """
    structures = _structures_from_objects(session, objects)
    if not structures:
        raise UserError("No atomic structures found in the requested model specification.")
    if dna_mode not in {"auto", "merge", "split"}:
        raise UserError("dna_mode must be one of: auto, merge, split")

    if output_dir:
        run_root = Path(output_dir).expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=True)
    else:
        run_root = Path(tempfile.mkdtemp(prefix="circoscontacts_"))

    session.logger.info(f"circoscontacts: writing intermediate files to {run_root}")

    base_spec = None
    if objects is not None and hasattr(objects, "atomspec"):
        base_spec = objects.atomspec
    elif objects is not None and hasattr(objects, "spec"):
        try:
            base_spec = str(objects.spec)
        except Exception:
            base_spec = None
    restrict_spec = None
    if restrict is not None and hasattr(restrict, "atomspec"):
        restrict_spec = restrict.atomspec
    elif restrict is not None and hasattr(restrict, "spec"):
        try:
            restrict_spec = str(restrict.spec)
        except Exception:
            restrict_spec = None
    if base_spec:
        session.logger.info(f"circoscontacts: source spec = {base_spec}")
    if restrict_spec:
        session.logger.info(f"circoscontacts: restrict spec = {restrict_spec}")

    contact_paths = _run_contacts_for_structures(
        session,
        structures,
        run_root,
        base_spec=base_spec,
        restrict_spec=restrict_spec,
        inter_model=inter_model,
        intra_mol=intramol,
    )
    if not contact_paths:
        raise UserError("No contact files produced.")

    source_chain_resnums = _objects_to_chain_resnums(objects)
    restrict_chain_resnums = _objects_to_chain_resnums(restrict)
    source_chain_hints, source_range_hints = _spec_chain_hints(base_spec)
    restrict_chain_hints, restrict_range_hints = _spec_chain_hints(restrict_spec)
    if source_range_hints:
        source_chain_resnums = source_range_hints
    if restrict_range_hints:
        restrict_chain_resnums = restrict_range_hints

    included_chain_ids: Set[str] = set()
    if source_chain_hints:
        included_chain_ids.update(source_chain_hints)
        included_chain_ids.update(restrict_chain_hints)
    if source_chain_resnums:
        included_chain_ids.update(source_chain_resnums.keys())
        included_chain_ids.update(restrict_chain_resnums.keys())

    chain_res, chain_resnums = _chain_residue_maps(structures)
    if included_chain_ids:
        chain_res = {c: r for c, r in chain_res.items() if c in included_chain_ids}
        chain_resnums = {c: n for c, n in chain_resnums.items() if c in included_chain_ids}
    dna_chains = circos.detect_dna_chains(chain_res)
    (
        pos_map,
        pos_info,
        display_chain_of,
        display_chains,
        start_pos_map,
    ) = circos.build_chain_maps(
        chain_res,
        chain_resnums,
        dna_chains,
        dna_reverse=None,
        dna_mismatches=dna_mismatches,
        dna_mode=dna_mode,
    )
    active_chain_resnums: Dict[str, Set[int]] = {}
    if source_chain_resnums:
        active_chain_resnums = source_chain_resnums
    elif restrict_chain_resnums:
        active_chain_resnums = restrict_chain_resnums
    active_ranges = _active_ranges_from_maps(active_chain_resnums, pos_map, display_chain_of)

    raw_counts_atom, raw_counts_residue = circos.parse_contacts(contact_paths)
    contacts, max_count_atom, max_count_residue = circos.aggregate_contacts(
        raw_counts_atom, raw_counts_residue, pos_map, display_chain_of
    )

    chain_lengths = {chain: len(pos_info.get(chain, [])) for chain in display_chains}
    chain_colors = circos.build_colors(display_chains)
    default_order = display_chains[:]

    html_path = run_root / "contacts_circos.html"
    try:
        circos.generate_html(
            title,
            display_chains,
            chain_lengths,
            chain_colors,
            pos_info,
            contacts,
            max_count_atom,
            max_count_residue,
            default_order,
            start_pos_map,
            active_ranges,
            html_path,
        )
    except Exception as err:
        session.logger.error(f"circoscontacts: HTML generation failed: {err}")
        raise UserError(f"Failed to generate circos HTML: {err}") from err

    if not html_path.exists():
        session.logger.error(
            f"circoscontacts: expected HTML was not created at {html_path}"
        )
        raise UserError(
            f"Circos HTML was not created. Expected file: {html_path}"
        )

    session.logger.info(f"circoscontacts: wrote {html_path}")
    session.logger.info(
        "circoscontacts: use the Save HTML button in the page to copy output elsewhere."
    )

    if open_html:
        _open_html(session, html_path)

    if not keep_temp and output_dir is None:
        session.logger.info(
            "circoscontacts: temporary directory retained for this run to keep artifacts accessible."
        )

    return str(html_path)


def register_command(command_name, logger):
    desc = CmdDesc(
        optional=[("objects", ObjectsArg)],
        keyword=[
            ("restrict", ObjectsArg),
            ("output_dir", StringArg),
            ("title", StringArg),
            ("inter_model", BoolArg),
            ("intramol", BoolArg),
            ("dna_mismatches", IntArg),
            ("dna_mode", StringArg),
            ("open_html", BoolArg),
            ("keep_temp", BoolArg),
        ],
        synopsis=(
            "Generate aggregated contact circos plot from open/selected models "
            "(default contacts behavior: protein, interModel false, intramol false)"
        ),
    )
    register(command_name, desc, circoscontacts, logger=logger)
