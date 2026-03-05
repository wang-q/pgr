# pgr maf

`pgr maf` provides tools for manipulating **MAF** (Multiple Alignment Format) files.

## Subcommands

*   `to-fas`: Convert MAF files to block FASTA format.

---

## to-fas

Converts MAF files into block FASTA format.

MAF files typically contain multiple sequence alignments. This command extracts each alignment block and converts it into the block FASTA format used by other `pgr` tools (like `pgr fas`).

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

1.  **Convert a MAF file to block FASTA format**:
    ```bash
    pgr maf to-fas tests/maf/example.maf
    ```
