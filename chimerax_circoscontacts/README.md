# ChimeraX CircosContacts Plugin (MVP)

This bundle adds a command:

`circoscontacts [objects] [restrict <spec>] [interModel true|false] [intramol true|false] [outputDir <path>] [title <text>] [dnaMismatches <N>] [dnaMode auto|merge|split] [openHtml true|false] [keepTemp true|false]`

## Default behavior

- If `objects` is omitted: uses all open atomic models.
- Contact generation defaults match your current workflow:
  - `contacts <spec> interModel false intramol false saveFile ...`
  - the default `<spec>` per model is `<model-spec> & protein`
- `restrict` is passed through to the underlying `contacts` command (per model).
- If a residue subset is used in `objects` (e.g. `#1/A,B:100-120`), circos arcs show
  full chain length with active source residues opaque and other regions dimmed.
- Outputs an interactive `contacts_circos.html` and opens it.
- Run output is written to a temp folder unless `output_dir` is provided.

## Example usage

- All open models:
  - `circoscontacts`
- Selected models/chains:
  - `circoscontacts #1,3,4 & protein`
- Contacts-style restricted run:
  - `circoscontacts #1,2/S,T:100-120 restrict /S,T/C interModel false intramol false`
- Save artifacts in an explicit directory:
  - `circoscontacts #1,3,4 & protein outputDir /tmp/cc_run`

## Install (development)

In ChimeraX command line:

`devel install /Users/pcarlton/WORK/Projects/KeitaPhD/tools/chimerax_circoscontacts`

Then run:

`circoscontacts`
