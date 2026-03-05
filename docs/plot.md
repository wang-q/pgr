# pgr plot

Plotting tools for various biological data visualizations.

`pgr plot` generates LaTeX source files (mostly using TikZ/PGFPlots) that can be compiled into high-quality PDFs. The recommended compiler is [Tectonic](https://tectonic-typesetting.github.io/), which automatically handles package downloads.

## Subcommands

| Subcommand | Description |
| :--- | :--- |
| `hh` | Histo-heatmap showing distribution of values across groups |
| `nrps` | NRPS (Non-Ribosomal Peptide Synthetase) structure diagram |
| `venn` | Venn diagram for 2-4 sets |

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
| `infile` | | | File | Input filename (`-` for stdin) |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `col` | `-c` | `--col` | Int | Column index to count (1-based, default: 1) |
| `group` | `-g` | `--group` | Int | Group column index (1-based) |
| `bins` | | `--bins` | Int | Number of bins (default: 40) |
| `xl` | | `--xl` | String | X axis label (default: column name) |
| `yl` | | `--yl` | String | Y axis label (default: group column name) |
| `xmm` | | `--xmm` | F,F | X axis range min,max (e.g., "0,100") |
| `unit` | | `--unit` | F,F | Cell width,height (default: "0.5,1.5") |

### Input Format

A tab-separated file with a header line.
*   **Column 1 (or specified by `--col`)**: Numeric values.
*   **Column 2 (or specified by `--group`)**: Group names (optional).

### Examples

```bash
# Basic usage
pgr plot hh input.tsv -o output.tex

# Compile directly with tectonic
pgr plot hh input.tsv | tectonic - && mv texput.pdf hh.pdf

# Specify columns and labels
pgr plot hh data.tsv -c 2 -g 1 --xl "Length" --yl "Species" -o plot.tex
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
| `infile` | | | File | Input filename (`-` for stdin) |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `legend` | | `--legend` | Flag | Include legend in the output |
| `color` | `-c` | `--color` | String | Default color (default: "grey") |

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
pgr plot nrps input.tsv --legend -c blue | tectonic -
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
