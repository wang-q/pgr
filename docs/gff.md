# pgr gff

`pgr gff` provides tools for manipulating **GFF** (General Feature Format) files.

## Subcommands

*   `rg`: Extract ranges and features from GFF files.

---

## rg

Extracts ranges and features from GFF files, optionally simplifying identifiers.

This command is useful for converting GFF annotations into a simpler tab-separated format (name, location) or extracting specific feature types.

```bash
pgr gff rg [OPTIONS] <infile>
```

### Options

*   `--tag <string>`: Feature type to retain (default: "gene").
*   `--asm <string>`: Assembly name (default: inferred from filename).
*   `--key <string>`: GFF attribute to use as the feature identifier (default: "ID").
    *   Choices: `ID`, `Name`, `Parent`, `gene`, `locus_tag`, `protein_id`, `product`.
*   `-s, --simplify`: Simplify sequence names (identifiers) by truncating at the first space, dot, comma, or dash.
*   `--ss`: Simplify reference sequence names (chromosome names) by truncating at the first space, dot, comma, or dash.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### Output Format

The output is a tab-separated list with two columns:
1.  Feature Key (e.g., Gene ID)
2.  Location: `Assembly.Chromosome(Strand):Start-End`

### Examples

1.  **Extract default 'gene' features**:
    ```bash
    pgr gff rg tests/gff/test.gff
    ```

2.  **Extract 'mRNA' features with a custom assembly name**:
    ```bash
    pgr gff rg tests/gff/test.gff --tag mRNA --asm Human -o output.tsv
    ```

3.  **Use 'Name' attribute as the identifier**:
    ```bash
    pgr gff rg tests/gff/test.gff --key Name
    ```

4.  **Simplify identifiers and sequence names**:
    ```bash
    pgr gff rg tests/gff/test.gff --simplify --ss
    ```
