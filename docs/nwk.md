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
    *   `to-svg`: Convert to SVG format.
    *   `to-tex`: Convert to a full LaTeX document.

---

## Info Commands

### cmp

Compares trees using Robinson-Foulds (RF) distance and its variants.

```bash
pgr nwk cmp [OPTIONS] <infile> [compare_file]
```

*   `[compare_file]`: Optional second file. If omitted, compares all trees in `<infile>` pairwise.
*   `--include-trivial`: Include trivial splits (single-leaf branches) in WRF/KF calculations.
*   Output columns: `Tree1`, `Tree2`, `RF_Dist`, `WRF_Dist`, `KF_Dist`.

### distance

Calculates distances between nodes or generates distance matrices.

```bash
pgr nwk distance [OPTIONS] <infile>
```

*   `--mode <mode>`: Calculation mode.
    *   `root` (default): Distance to root.
    *   `parent`: Distance to parent.
    *   `pairwise`: All pairwise distances.
    *   `lca`: Distances to the Lowest Common Ancestor.
    *   `phylip`: Phylip-formatted distance matrix for the selected nodes.
*   `-I`: Ignore internal nodes.
*   `-L`: Ignore leaf nodes.

### label

Extracts labels (names) from the tree.

```bash
pgr nwk label [OPTIONS] <infile>
```

*   `-I`: Skip internal labels.
*   `-L`: Skip leaf labels.
*   `-n <name>` / `-l <file>` / `-x <regex>`: Filter nodes.
*   `-D`: Include descendants of selected nodes.
*   `-M`: Monophyly check (only print if selected nodes form a monophyletic group).
*   `--tab`: Tab-separated output (single line).
*   `-c <col>` / `--extra-column <col>`: Add extra columns (`dup`, `taxid`, `species`, `full`).
*   `--root`: Only print the root label.

### stat

Prints statistics about the trees.

```bash
pgr nwk stat [OPTIONS] <infile>
```

*   `--style <col|line>`: Output format (key-value pairs or TSV).
*   Statistics include: Type (cladogram/phylogram/neither), Node count, Leaf count, Rooted status, Dichotomies, Leaf labels, Internal labels, Cherries, Sackin index, Colless index.

### support

Attributes bootstrap support values to a target tree based on replicate trees.

```bash
pgr nwk support [OPTIONS] <target> <replicates>
```

*   `-p, --percent`: Output support values as percentages (0-100).

---

## Ops Commands

### order

Sorts the children of each node (rotates branches) without changing topology.

```bash
pgr nwk order [OPTIONS] <infile>
```

*   `--num-descendants` / `--num-descendants-rev`: Sort by number of descendants (Ladderize).
*   `--alphanumeric` / `--alphanumeric-rev`: Sort by label (Alphanumeric).
*   `--name-list <file>`: Sort by a list of names.
*   `--deladderize` (`--dl`): Alternate sort direction at each level.

### prune

Removes nodes from the tree.

```bash
pgr nwk prune [OPTIONS] <infile>
```

*   `-n <name>` / `-l <file>` / `-x <regex>`: Select nodes to remove.
*   `-i, --invert`: Invert selection (Keep selected nodes, remove others).
*   `-D`: Include descendants.

### rename

Renames specific nodes.

```bash
pgr nwk rename [OPTIONS] <infile>
```

*   `-n <name>`: Select node by name.
*   `-l <name1,name2>`: Select node by LCA of two names.
*   `--rename <new_name>`: New name (must correspond to `-n` or `-l` arguments).

### replace

Batch replaces node names or annotations using a TSV file.

```bash
pgr nwk replace [OPTIONS] --replace-tsv <replace.tsv> <infile>
```

*   `--replace-tsv <replace.tsv>`: Tab-separated file: `Original <TAB> Replacement [TAB Extra...]`.
*   `-I, --internal`: Skip internal labels.
*   `-L, --leaf`: Skip leaf labels.
*   `--mode <mode>`:
    *   `label` (default): Replace node name.
    *   `taxid`: Add as NCBI TaxID (`:T=`).
    *   `species`: Add as species name (`:S=`).
    *   `asis`: Append as comments/properties. Values containing `=` are parsed as `key=value` pairs; bare values are stored as keys with empty values.

### reroot

Reroots the tree.

```bash
pgr nwk reroot [OPTIONS] <infile>
```

*   (Default): Reroot at the midpoint of the longest branch.
*   `-n <node>`: Reroot at the edge leading to the LCA of specified nodes (Ingroup).
*   `-l, --lax`: Lax mode (use complement if LCA is already root).
*   `-d, --deroot`: Deroot (create multifurcating root).
*   `--support-as-labels`: Support values as labels (shift labels during reroot).

### subtree

Extracts a subtree (clade) rooted at the LCA of selected nodes.

```bash
pgr nwk subtree [OPTIONS] <infile>
```

*   `-n` / `-l` / `-x`: Select nodes.
*   `-D, --descendants`: Include all descendants of selected internal nodes.
*   `-M`: Monophyly check (only output if clade matches selection exactly).
*   `-c <N>`: Context (extend N levels up).
*   `-C, --condense <name>`: Condense the subtree into a single node.
    * The new node is annotated with `member=<count>` and `tri=white`.
    * `<count>` is the number of named nodes matched by `-n/-l/-x` (including descendants expanded by `-D`).

### topo

Modifies tree topology and attributes.

```bash
pgr nwk topo [OPTIONS] <infile>
```

*   By default, removes branch lengths and comments (topology only).
*   `-b, --bl`: Keep branch lengths.
*   `-c, --comment`: Keep comments.
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
*   `--string <str>`: Add free-form string annotations.
*   `--color`, `--label`, `--comment-text`: Add text attributes.
*   `--dot`, `--bar`, `--rec`, `--tri`: Add shape attributes (for `to-tex` / `to-forest`).
*   `--remove <regex>`: Remove matching comments.
*   `-o, --outfile <file>`: Output filename. `[stdout]` for screen.

### indent

Formats Newick trees with indentation.

```bash
pgr nwk indent [OPTIONS] <infile>
```

*   `--text <str>`: Indentation string (default: "  ").
*   `-c, --compact`: Compact output (single line).
*   `-o, --outfile <file>`: Output filename. `[stdout]` for screen.

### to-dot

Converts Newick trees to Graphviz DOT format.

```bash
pgr nwk to-dot [OPTIONS] <infile>
```

*   `-o, --outfile <file>`: Output filename. `[stdout]` for screen.

### to-forest

Converts Newick trees to raw LaTeX Forest code.

```bash
pgr nwk to-forest [OPTIONS] <infile>
```

*   `-b, --bl`: Include branch lengths.
*   `-o, --outfile <file>`: Output filename. `[stdout]` for screen.

### to-svg

Converts Newick trees to SVG format for visualization.

```bash
pgr nwk to-svg [OPTIONS] <infile>
```

*   Automatically draws a phylogram if branch lengths are present, otherwise a cladogram.
*   `-w, --width <N>`: SVG width in pixels (default: 800).
*   `-v, --vskip <N>`: Vertical spacing between leaf nodes in pixels (default: 20).
*   `-o, --outfile <file>`: Output filename. `[stdout]` for screen.

### to-tex

Converts Newick trees to a full LaTeX document (wrapper around `to-forest`).

```bash
pgr nwk to-tex [OPTIONS] <infile>
```

*   `-b, --bl`: Draw phylogram (with branch lengths).
*   `--forest`: Input is already Forest code.
*   `--no-default-style`: Skip default style definitions.
*   `-o, --outfile <file>`: Output filename. `[stdout]` for screen.

---

## Planned Subcommands

*   `eval` [Planned]: Multi-dimensional tree evaluation framework (geometric, taxonomic, phylogenetic, trait consistency). Referenced by `pgr clust eval` and `pgr clust cut` for tree-based metrics.
