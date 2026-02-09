# ChimeraX Contacts Circos Tool

Generate an interactive, presentation-quality circos-style plot from ChimeraX `.contacts` files and their corresponding AlphaFold/AF3 `.cif` models.

This tool:
- Aggregates contacts across multiple `.contacts` files.
- Treats DNA as a single duplex chromosome (base-pair positions), reversing one strand.
- Provides an interactive HTML with controls for visibility threshold, line width, rotation, and chain ordering.
- Exports an SVG snapshot from the browser.

## Requirements

- Python 3.10+
- `gemmi`

Install dependency:

```bash
pip install gemmi
```

## Usage

From the project root:

```bash
python tools/contacts_circos.py \
  --contacts-dir data/chimerax/SPORM/spo-dsb-full-cifs/correct \
  --title "SPO-11 DNA Contacts"
```

The output HTML is written to `contacts_circos.html` in the contacts directory.

### DNA orientation

By default, when two DNA chains are detected the tool:
- Reverses the first detected DNA chain, and
- Uses the second as the forward strand.

To override which chain is reversed:

```bash
python tools/contacts_circos.py \
  --contacts-dir data/chimerax/SPORM/spo-dsb-full-cifs/correct \
  --dna-reverse J
```

## Output

Open the generated HTML in a browser. Controls allow:
- Thresholding contacts by count
- Rotating the plot
- Adjusting max line width
- Reordering chains
- Downloading SVG

## Notes

This repository is intentionally scoped to the `tools/` directory only.
