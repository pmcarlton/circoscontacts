# ChimeraX Contacts Circos Tool

Generate an interactive, presentation-quality circos-style plot from ChimeraX `.contacts` files and their corresponding AlphaFold/AF3 `.cif` models.

[Homepage and documentation](https://circoscontacts.carltonlab.org)

This tool:
- Aggregates contacts across multiple `.contacts` files.
- Treats DNA as a single duplex chromosome (base-pair positions), reversing one strand.
- Provides an interactive HTML with controls for visibility threshold, line width, rotation, zoom, and chain ordering.
- Supports click-and-hold filtering by chain, and shift-click-and-hold filtering by residue.
- Shows tooltip residue details on hover and subtle arc start/end styling.
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
python tools/contacts_circos.py --title "SPO-11 DNA Contacts"
```

To target a specific directory of `.contacts` and `.cif` files:

```bash
python tools/contacts_circos.py \
  --contacts-dir data/chimerax/SPORM/spo-dsb-full-cifs/correct \
  --title "SPO-11 DNA Contacts"
```

The output HTML is written to `contacts_circos.html` in the contacts directory (default: current directory).

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
- Zooming the plot (panel scrolls when zoomed in)
- Reordering chains
- Downloading SVG
 - Click-and-hold on an arc: show only contacts involving that chain
 - Shift-click-and-hold on an arc: show only contacts involving that residue
