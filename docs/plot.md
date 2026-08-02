# pgr plot

Plotting tools for various biological data visualizations.

`pgr plot` generates figures in two styles:

* LaTeX source files (TikZ/PGFPlots) for the `hh`, `nrps`, and `venn` subcommands,
  compiled into PDFs with [Tectonic](https://tectonic-typesetting.github.io/).
* Standalone SVG for the `dot` subcommand, rendered with no external dependencies.

## Subcommands

| Subcommand | Description |
| :--- | :--- |
| `dot` | Dot plot (collinear plot) of PAF alignments, output as SVG |
| `hh` | Histo-heatmap showing distribution of values across groups |
| `nrps` | NRPS (Non-Ribosomal Peptide Synthetase) structure diagram |
| `venn` | Venn diagram for 2-4 sets |

---

## dot

Dot plot (collinear plot) of PAF alignments, output as an SVG file.

The two axes are the target (x) and query (y) sequences, laid out contig by
contig in first-appearance order. Each alignment is drawn as a line segment
colored by identity on a blue scale: `--min-identity` (lightest) to
`--max-identity` (deepest), with a legend at the bottom right. Axes carry bp
tick marks with real genomic coordinates, light grid lines, and contig
separators. All line widths and font sizes scale with `--width`.

### Usage

```bash
pgr plot dot [OPTIONS] <infile>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infile` | | | File | Input PAF file (".paf.gz" supported, "stdin" for standard input) |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `min_len` | | `--min-len` | Int | Minimum alignment block length (default: 100) |
| `min_identity` | | `--min-identity` | Float | Minimum identity to plot, matches / block length; lightest color (default: 0.7) |
| `identity_max` | | `--max-identity` | Float | Identity at which the color scale saturates (default: 1.0) |
| `max_align` | | `--max-align` | Int | Maximum number of alignments to plot, longest first (default: 100000; 0 = all) |
| `width` | | `--width` | Int | Plot width in pixels; height is scaled automatically (default: 1200) |
| `range` | | `--range` | String | Target-side region to zoom into, `chr:start-end` (1-based); the query axis auto-focuses on the main aligned band |

### Notes

* The SVG is self-contained and needs no external libraries to view.
* Convert to PDF or PNG with `rsvg-convert`, `inkscape`, or `cairosvg` when needed.
* Line widths and label/tick font sizes are derived from `--width`, so the
  plot scales consistently at any size.
* Alignments below `--min-identity` or shorter than `--min-len` are skipped;
  when the input exceeds `--max-align`, only the longest alignments are drawn
  to keep the file size manageable.
* With `--range`, only alignments overlapping the target-side region are kept
  and clipped to it; the query axis shows every significant aligned cluster
  (aligned bases >= 1% of the largest cluster), each as its own segment with
  true genomic coordinates. Matches far away on the same chromosome or on
  other chromosomes stay visible; only tiny noise fragments are dropped.

### Examples

```bash
pgr plot dot input.paf -o dot.svg

pgr plot dot input.paf.gz | rsvg-convert -f pdf -o dot.pdf

pgr plot dot input.paf --range chr1:100000-200000 -o zoom.svg
```

---

## hh

Histo-heatmap. This visualization combines a histogram and a heatmap to show the distribution of a numeric variable (X) across different groups (Y).

### Usage

```bash
pgr plot hh [OPTIONS] <infile>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infile` | | | File | Input filename ("stdin" for standard input) |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `column` | `-c` | `--column` | Int | Column index to count (1-based, default: 1) |
| `group` | `-g` | `--group` | Int | Group column index (1-based) |
| `bins` | | `--bins` | Int | Number of bins (default: 40) |
| `xl` | | `--xlabel` | String | X axis label (default: column name) |
| `yl` | | `--ylabel` | String | Y axis label (default: group column name) |
| `xmm` | | `--xmin-max` | F,F | X axis range min,max (e.g., "0,100") |
| `unit` | | `--unit` | F,F | Cell width,height (default: "0.5,1.5") |

### Input Format

A tab-separated file with a header line.
*   **Column 1 (or specified by `--column`)**: Numeric values.
*   **Column 2 (or specified by `--group`)**: Group names (optional).

### Examples

```bash
# Basic usage
pgr plot hh input.tsv -o output.tex

# Compile directly with tectonic
pgr plot hh input.tsv | tectonic - && mv texput.pdf hh.pdf

# Specify columns and labels
pgr plot hh data.tsv -c 2 -g 1 --xlabel "Length" --ylabel "Species" -o plot.tex
```

---

## nrps

Generates a structural diagram for Non-Ribosomal Peptide Synthetase (NRPS) modules and domains.

### Usage

```bash
pgr plot nrps [OPTIONS] <infile>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infile` | | | File | Input filename ("stdin" for standard input) |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `legend` | | `--legend` | Flag | Include legend in the output |
| `color` | | `--color` | String | Default color (default: "grey") |

### Input Format

A tab-separated file defining modules and domains.

*   **Module Definition**: `Module <Name> <Color>`
    *   Starts a new module.
    *   Color is optional.
*   **Domain Definition**: `<Type> <Text> <Color>`
    *   **Type**: A, C, E, CE, T, Te, R, M.
    *   **Text**: Optional label inside the domain.
    *   **Color**: Optional override.

**Supported Colors**: black, grey, red, brown, green, purple, blue.

### Examples

```bash
# Generate diagram
pgr plot nrps input.tsv -o nrps.tex

# With legend and custom default color
pgr plot nrps input.tsv --legend --color blue | tectonic -
```

---

## venn

Generates a Venn diagram for 2, 3, or 4 sets.

### Usage

```bash
pgr plot venn [OPTIONS] <infiles>...
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infiles` | | | Files | 2 to 4 input list files |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |

### Input Format

Plain text files, each containing a list of unique items (one per line). The filename (without extension) is used as the set label.

### Examples

```bash
# 2 sets
pgr plot venn list1.txt list2.txt -o venn2.tex

# 3 sets
pgr plot venn A.txt B.txt C.txt -o venn3.tex

# 4 sets
pgr plot venn A.list B.list C.list D.list | tectonic -
```
