# pgr sam

`pgr sam` provides tools for manipulating **SAM** alignment files.

## Subcommands

*   `to-rg`: Extract alignment coordinates as ranges (.rg).

---

## to-rg

Extracts alignment coordinates from SAM files and outputs them in `.rg`
format (`chr:start-end`, 1-based inclusive). This is the bridge for
deriving per-base coverage from the mapped SAM of `pgr asm map` (anchr
`anchors` step): each mapped record becomes one range, and `pgr rg
coverage` computes the depth over all ranges.

```bash
pgr sam to-rg [OPTIONS] <infile>
```

### Options

*   `-o, --outfile <file>`: Output filename. `[stdout]` for screen.
*   `--strict`: Fail on parse errors instead of skipping malformed lines.

### Notes

*   Header lines (`@...`) and unmapped records (FLAG 0x4, RNAME `*`, or
    POS 0) are skipped.
*   The output range spans the full reference-consuming CIGAR
    (`M`/`D`/`N`/`=`/`X`); insertions, soft/hard clips, and padding do not
    extend the range.
*   Supports both plain text and gzipped (.gz) files.
*   Reads from stdin if input file is `stdin`.

### Examples

1.  **Convert mapped reads to ranges**:
    ```bash
    pgr sam to-rg mapped.sam > mapped.rg
    ```

2.  **Derive per-base coverage from an `asm map` SAM (anchr anchors step)**:
    ```bash
    pgr asm map UT.fasta R1.fq.gz R2.fq.gz --outm mapped.sam --outu unmapped.sam
    pgr sam to-rg mapped.sam | pgr rg coverage stdin -m 2 -o cov.json
    ```
