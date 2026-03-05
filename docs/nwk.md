# pgr nwk

`pgr nwk` provides a comprehensive suite of tools for manipulating, analyzing, and visualizing phylogenetic trees in **Newick** format.

## Subcommands

The subcommands are organized into the following categories:

*   **Info**: Extract information or statistics.
    *   `cmp`: Compare trees (RF, WRF, KF distances).
    *   `distance`: Calculate distances between nodes.
    *   `label`: Extract labels/names from the tree.
    *   `stat`: Print tree statistics (nodes, leaves, balance indices).
    *   `support`: Attribute bootstrap support values.
*   **Ops**: Manipulate tree structure.
    *   `cut`: Cut the tree into clusters.
    *   `order`: Reorder nodes (ladderize, alphanumeric).
    *   `prune`: Remove nodes (and their descendants/ancestors).
    *   `rename`: Rename specific nodes.
    *   `replace`: Batch replace names or comments using a file.
    *   `reroot`: Reroot the tree.
    *   `subtree`: Extract a specific clade/subtree.
    *   `topo`: Modify topology (remove branch lengths/labels).
*   **Viz**: Visualization and formatting.
    *   `comment`: Add annotations/comments for visualization.
    *   `indent`: Reformat Newick with indentation.
    *   `to-dot`: Convert to Graphviz DOT format.
    *   `to-forest`: Convert to LaTeX Forest code.
    *   `to-tex`: Convert to a full LaTeX document.

---

## Info Commands

### cmp

Compares trees using Robinson-Foulds (RF) distance and its variants.

```bash
pgr nwk cmp [OPTIONS] <infile> [compare_file]
```

*   `[compare_file]`: Optional second file. If omitted, compares all trees in `<infile>` pairwise.
*   Output columns: `Tree1_ID`, `Tree2_ID`, `RF_Dist`, `WRF_Dist`, `KF_Dist`.

### distance

Calculates distances between nodes or generates distance matrices.

```bash
pgr nwk distance [OPTIONS] <infile>
```

*   `-m, --mode <mode>`: Calculation mode.
    *   `root` (default): Distance to root.
    *   `parent`: Distance to parent.
    *   `pairwise`: All pairwise distances.
    *   `lca`: Distances to the Lowest Common Ancestor.
    *   `phylip`: Output Phylip-formatted distance matrix.
*   `-I`: Ignore internal nodes.
*   `-L`: Ignore leaf nodes.

### label

Extracts labels (names) from the tree.

```bash
pgr nwk label [OPTIONS] <infile>
```

*   `-I`: Skip internal labels.
*   `-L`: Skip leaf labels.
*   `-n <name>` / `-f <file>` / `-r <regex>`: Filter nodes.
*   `-D`: Include descendants of selected nodes.
*   `-M`: Monophyly check (only print if selected nodes form a monophyletic group).
*   `-t`: Tab-separated output (single line).
*   `-c <col>`: Add extra columns (`dup`, `taxid`, `species`, `full`).

### stat

Prints statistics about the trees.

```bash
pgr nwk stat [OPTIONS] <infile>
```

*   `--style <col|line>`: Output format (key-value pairs or TSV).
*   Statistics include: Node count, Leaf count, Rooted status, Dichotomies, Label counts, Cherries, Sackin index, Colless index.

### support

Attributes bootstrap support values to a target tree based on replicate trees.

```bash
pgr nwk support [OPTIONS] <target> <replicates>
```

*   `-p, --percent`: Output support values as percentages (0-100).

---

## Ops Commands

### cut

Cuts the tree into clusters based on various criteria.

```bash
pgr nwk cut [OPTIONS] <infile>
```

*   `--k <N>`: Cut into K clusters.
*   `--height <H>`: Cut at specific height (max distance to leaves).
*   `--root-dist <D>`: Cut at specific distance from root.
*   `--max-clade <T>`: Max pairwise distance in cluster <= T.
*   `--avg-clade <T>`: Avg pairwise distance in cluster <= T.
*   `--inconsistent <T>`: Inconsistent coefficient <= T.
*   `--rep <root|first|medoid>`: Representative selection method.
*   `--format <cluster|pair>`: Output format.

### order

Sorts the children of each node (rotates branches) without changing topology.

```bash
pgr nwk order [OPTIONS] <infile>
```

*   `--nd` / `--ndr`: Sort by number of descendants (Ladderize).
*   `--an` / `--anr`: Sort by label (Alphanumeric).
*   `--list <file>`: Sort by a list of names.
*   `--deladderize`: Alternate sort direction at each level.

### prune

Removes nodes from the tree.

```bash
pgr nwk prune [OPTIONS] <infile>
```

*   `-n <name>` / `-f <file>` / `-r <regex>`: Select nodes to remove.
*   `-v, --invert`: Invert selection (Keep selected nodes, remove others).
*   `-D`: Include descendants.

### rename

Renames specific nodes.

```bash
pgr nwk rename [OPTIONS] <infile>
```

*   `-n <name>`: Select node by name.
*   `-l <name1,name2>`: Select node by LCA of two names.
*   `-r <new_name>`: New name (must correspond to `-n` or `-l` arguments).

### replace

Batch replaces node names or annotations using a TSV file.

```bash
pgr nwk replace [OPTIONS] <infile> <replace.tsv>
```

*   `<replace.tsv>`: Tab-separated file: `Original <TAB> Replacement [TAB Extra...]`.
*   `--mode <mode>`:
    *   `label` (default): Replace node name.
    *   `taxid`: Add as NCBI TaxID (`:T=`).
    *   `species`: Add as species name (`:S=`).
    *   `asis`: Append verbatim to comments.

### reroot

Reroots the tree.

```bash
pgr nwk reroot [OPTIONS] <infile>
```

*   (Default): Reroot at the midpoint of the longest branch.
*   `-n <node>`: Reroot at the edge leading to the LCA of specified nodes (Ingroup).
*   `-l, --lax`: Lax mode (use complement if LCA is already root).
*   `-d, --deroot`: Deroot (create multifurcating root).
*   `-s`: Support values as labels (shift labels during reroot).

### subtree

Extracts a subtree (clade) rooted at the LCA of selected nodes.

```bash
pgr nwk subtree [OPTIONS] <infile>
```

*   `-n` / `-f` / `-r`: Select nodes.
*   `-M`: Monophyly check (only output if clade matches selection exactly).
*   `-c <N>`: Context (extend N levels up).
*   `-C, --condense <name>`: Condense the subtree into a single node.

### topo

Modifies tree topology and attributes.

```bash
pgr nwk topo [OPTIONS] <infile>
```

*   By default, removes branch lengths and comments (topology only).
*   `--bl`: Keep branch lengths.
*   `--comment`: Keep comments.
*   `-I`: Remove internal labels.
*   `-L`: Remove leaf labels.

---

## Viz Commands

### comment

Adds annotations/comments to nodes for visualization.

```bash
pgr nwk comment [OPTIONS] <infile>
```

*   `-n` / `-l`: Select nodes.
*   `--color`, `--label`, `--comment`: Add text attributes.
*   `--dot`, `--bar`, `--rec`, `--tri`: Add shape attributes (for `to-tex` / `to-forest`).
*   `-r, --remove <regex>`: Remove matching comments.

### indent

Formats Newick trees with indentation.

```bash
pgr nwk indent [OPTIONS] <infile>
```

*   `-t, --text <str>`: Indentation string (default: "  ").
*   `-c, --compact`: Compact output (single line).

### to-dot

Converts Newick trees to Graphviz DOT format.

```bash
pgr nwk to-dot [OPTIONS] <infile>
```

### to-forest

Converts Newick trees to raw LaTeX Forest code.

```bash
pgr nwk to-forest [OPTIONS] <infile>
```

*   `--bl`: Include branch lengths.

### to-tex

Converts Newick trees to a full LaTeX document (wrapper around `to-forest`).

```bash
pgr nwk to-tex [OPTIONS] <infile>
```

*   `--bl`: Draw phylogram (with branch lengths).
*   `--forest`: Input is already Forest code.
*   `--style`: Skip default style definitions.
