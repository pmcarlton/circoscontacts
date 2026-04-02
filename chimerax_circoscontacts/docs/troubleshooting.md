# Troubleshooting

## `Missing dependency on ChimeraX version` during upload

Use wheel `0.3.3+` (or newer), which declares explicit ChimeraX dependencies in bundle metadata.

## `HTML generation failed: 'ascii' codec can't encode character ...`

Use wheel `0.3.4+` (or newer).

Cause: older builds wrote HTML with locale default encoding in some ChimeraX runtimes.
Fix: HTML write now forces UTF-8.

## Output HTML opens but no links are visible

Check:

- threshold is not too high,
- chain set/order still includes linked chains,
- per-chain contact checkboxes are enabled,
- `Contact mode` matches what you expect (`Atom` vs `Residue`).

## Shift-click isolate feels inconsistent

Shift-click isolate only acts on residues with links visible under current threshold and chain visibility state.
Hover tooltip blue background indicates residue isolate is currently actionable.

## Lock-to-bottom not staying active

Bottom lock is intentionally cleared when:

- another chain label is locked,
- locked chain is removed from `Chain Presence/Order`,
- angle offset is changed manually.

## `.cxc` color script does not appear to color expected residues

Check:

- exported with desired threshold and chain visibility,
- script run on matching models/chains,
- residue numbering/chain IDs match the originating structures.

## Session file fails to load

- Ensure file is a plugin session JSON (`contacts-circos-session-v1`).
- Invalid or older schema files are rejected.
