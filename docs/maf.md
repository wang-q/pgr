# pgr maf

`pgr maf` provides tools for manipulating **MAF** (Multiple Alignment Format) files.

MAF is treated as an **exchange format** with external tools: upstream aligners
(lastz/FastGA) produce MAF, and downstream ecosystems consume it. Inside `pgr`,
multiple-sequence alignment and core-genome work happens in **Block FA** (the
`fas` command family, e.g. `pgr pl p2m`), and the pangenome implicit graph works
in PAF. Therefore `pgr maf` intentionally stays thin — conversions only, no
filter/subset commands (use `pgr fas` for those).

## Subcommands

*   `to-fas`: Convert MAF files to block FA format.
*   `to-paf`: Convert two-sequence MAF files to PAF format.

---

## to-fas

Converts MAF files into block FA format.

MAF files typically contain multiple sequence alignments. This command extracts each alignment block and converts it into the block FA format used by other `pgr` tools (like `pgr fas`).

```bash
pgr maf to-fas [OPTIONS] <infiles>...
```

### Options

*   `-o, --outfile <file>`: Output filename (default: stdout).

### Notes

*   Supports both plain text and gzipped (`.gz`) input files.
*   Reads from stdin if the input file is `stdin`.
*   The output preserves the alignment structure, with each block separated by a blank line.

### Examples

1.  **Convert a MAF file to block FA format**:
    ```bash
    pgr maf to-fas tests/maf/example.maf
    ```

---

## to-paf

Converts MAF (Multiple Alignment Format) files containing pairwise alignments into PAF (Pairwise mApping Format).

Only blocks with exactly two `s` lines are converted. Multi-sequence blocks are skipped with a warning.

```bash
pgr maf to-paf [OPTIONS] <infiles>...
```

### Options

*   `-o, --outfile <file>`: Output filename (default: stdout).

### Custom PAF Tags

*   `cg:Z:` – CIGAR string derived from the MAF alignment strings.
*   `gi:f:` – Gap-compressed identity.
*   `bi:f:` – Block identity.
*   `ms:i:` – MAF score (from the `a` line `score=` field).

### Notes

*   Supports both plain text and gzipped (`.gz`) input files.
*   Reads from stdin if the input file is `stdin`.

### Examples

1.  **Convert a MAF file to PAF**:
    ```bash
    pgr maf to-paf ref_vs_query.maf -o ref_vs_query.paf
    ```
