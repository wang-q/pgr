# pgr sam

`pgr sam` provides tools for manipulating **SAM** alignment files.

## Subcommands

*   `ihist`: Compute an insert-size histogram from a paired SAM.
*   `to-rg`: Extract alignment coordinates as ranges (.rg).

---

## ihist

Reads a paired SAM (e.g. the `--paired` output of `anchr asm map`) and writes
the insert-size histogram in the BBTools `reformat.sh ihist` text format:
`#Mean`/`#Median`/`#Mode`/`#STDev`/`#PercentOfPairs` lines followed by
`#InsertSize  Count` rows. This replaces the `reformat.sh ihist` call of
the anchr `2_insert_size` step.

```bash
pgr sam ihist [OPTIONS] <infile>
```

### Options

*   `-o, --outfile <file>`: Output filename. `[stdout]` for screen.

### Notes

*   Pairs are grouped by read name (first whitespace token, trailing
    `/1`/`/2` stripped).
*   Only proper FR pairs — both ends mapped, same reference, opposite
    strands, pointing inward — contribute an insert size.
*   `#PercentOfPairs` is the fraction of pairs contributing to the
    histogram.
*   The median is the lower median; the mode is the most frequent size
    (ties -> the smallest); the standard deviation is population.
*   Supports both plain text and gzipped (.gz) files.
*   Reads from stdin if input file is `stdin`.

### Examples

1.  **Insert-size histogram from a paired mapping (anchr 2_insert_size
    step)**:
    ```bash
    anchr asm map UT.fasta R1.fq.gz R2.fq.gz --paired \
        --outm mapped.sam --outu unmapped.sam --max-reads 1000000
    pgr sam ihist mapped.sam -o insert_size.ihist.txt
    ```

---

## to-rg

Extracts alignment coordinates from SAM files and outputs them in `.rg`
format (`chr:start-end`, 1-based inclusive). This is the bridge for
deriving per-base coverage from the mapped SAM of `anchr asm map` (anchr
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
    anchr asm map UT.fasta R1.fq.gz R2.fq.gz --outm mapped.sam --outu unmapped.sam
    pgr sam to-rg mapped.sam | pgr rg coverage stdin -m 2 -o cov.json
    ```
