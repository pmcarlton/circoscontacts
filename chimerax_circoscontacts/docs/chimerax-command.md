# ChimeraX Command

## Syntax

```chimerax
circoscontacts [objects] [restrict <spec>] [outputDir <path>] [title <text>] \
               [interModel true|false] [intramol true|false] \
               [dnaMismatches <N>] [dnaMode auto|merge|split] \
               [openHtml true|false] [keepTemp true|false]
```

Run `usage circoscontacts` in ChimeraX for the live installed signature.

## Contact Semantics (matches `contacts` workflow)

The plugin internally runs `contacts` per structure and aggregates outputs.

### Default source set

If `objects` is omitted:

- all open atomic models are used,
- per model, source spec defaults to `<model-spec> & protein`,
- defaults are `interModel false` and `intramol false`.

Equivalent conceptual command per model:

```chimerax
contacts <model-spec> & protein interModel false intramol false saveFile <...>
```

### When `objects` is provided

For each selected structure, source is intersected with that structure:

```text
(<objects>) & (<that model>)
```

This mirrors normal ChimeraX object-spec restriction behavior.

### `restrict` behavior

`restrict` is passed through to `contacts`, also intersected per model:

```text
restrict ((<restrict spec>) & (<that model>))
```

This gives the same “source vs restrict” logic as direct `contacts` usage.

### Direction and counting

- Contact direction is ignored for aggregation (A→B equals B→A).
- Counts are maintained for both atom-level and residue-collapsed representations.

## DNA Handling

`circoscontacts` can treat dsDNA as one logical arc.

### `dnaMode`

- `auto`: detect and merge recognizable dsDNA layouts
- `merge`: force merged DNA behavior when possible
- `split`: keep chains independent (no DNA merge)

### `dnaMismatches`

Used in split/nicked DNA detection where reverse-complement matching is tested.

- `0` (default): strict reverse-complement match
- `N > 0`: allow up to `N` mismatches during merge detection

## Output Control

- `outputDir <path>`: write HTML/intermediate files to a stable directory
- `openHtml true|false`: auto-open final HTML in browser
- `keepTemp true|false`: keep temporary run artifacts when `outputDir` is omitted

## Examples

### 1) Default run on all open structures

```chimerax
circoscontacts
```

### 2) Selected models, protein source only

```chimerax
circoscontacts #1,3,4 & protein
```

### 3) Contact-style source/restrict windowing

```chimerax
circoscontacts #1,2/S,T:100-120 restrict /S,T/C interModel false intramol false
```

Interpretation:

- source: residues 100–120 in chains `S,T` (within models #1,#2),
- restrict target: chains `S,T,C`,
- self-chain exclusions remain governed by `contacts` settings.

### 4) Batch output to explicit directory

```chimerax
circoscontacts #1,2 outputDir /tmp/cc_run title "SPO11 DNA contacts"
```

## Notes for Reproducibility

- Use explicit `objects`, `restrict`, and `outputDir` in saved scripts.
- Export a `.json` session from the web UI to preserve all display and annotation state.
