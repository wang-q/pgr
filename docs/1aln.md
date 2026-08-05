# pgr 1aln

`pgr 1aln` reads **.1aln** files — the compact binary alignment container from
FastGA, built on the ONEcode container format. Alignments are stored as sampled
**trace points** rather than base-level columns, so the file is much smaller
than PAF while preserving the aligned path. Each alignment is expanded back to
base-level columns on read.

The `.1aln` header stores only the source genome file *references* (not the
sequences), so the read-side commands require the two source genomes.

## Subcommands

*   `stat`: Report header and per-record statistics.
*   `to-paf`: Expand alignments to PAF format.
*   `to-psl`: Expand alignments to PSL format.

The write-side mirrors live on the source-format commands:
`pgr paf to-1aln` and `pgr maf to-1aln`.

---

## stat

Reports a TSV (`key<TAB>value`) summary of the `.1aln` header and record
statistics. No source sequences are needed.

```bash
pgr 1aln stat [OPTIONS] <infile>
```

### Options

*   `-o, --outfile <file>`: Output filename (default: stdout).

### Reported fields

*   `tspace` – trace point spacing (the `t` line).
*   `records` – number of `A` alignment records (from the footer counts).
*   `max_trace_points` – largest trace-point count across records.
*   `total_diffs` – summed differences across records.
*   `skeletons` – number of GDB skeleton objects.
*   `scaffolds`, `contigs` – across all skeletons.
*   `refs` – number of reference (`<`) entries, each with `count` 1/2/3.
*   `ref.N.filename` / `ref.N.count` – per-reference source file.
*   `prov.N.program` / `prov.N.version` – provenance.

### Notes

*   Reads a single `.1aln` file; does not support gzip or stdin (the ONEcode
    container requires random access to the footer offset at EOF).

### Examples

1.  **Report header stats**:
    ```bash
    pgr 1aln stat mg1655-sakai.1aln -o stats.tsv
    ```

---

## to-paf

Expands each alignment record back into base-level aligned columns and emits a
PAF record per alignment.

```bash
pgr 1aln to-paf [OPTIONS] <infile>
```

### Options

*   `--ref-seq <file>`: Reference (a side) genome FASTA / gzipped FASTA.
*   `--query-seq <file>`: Query (b side) genome FASTA / gzipped FASTA.
*   `--cigar`: Emit the `cg:Z` CIGAR tag.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### Custom PAF tags

*   `dv:f:` – identity (matches / aligned bases).
*   `df:i:` – number of differences (substitutions + indels).
*   `cg:Z:` – X-CIGAR (only with `--cigar`).

### Notes

*   Requires `--ref-seq` and `--query-seq` (the source genomes).
*   A `-` strand means the `b`-side sequence was stored reverse-complemented;
    its PAF coordinates are given in forward orientation.
*   Reads a single `.1aln` file; does not support gzip or stdin.

### Examples

1.  **Convert with default tags (no CIGAR)**:
    ```bash
    pgr 1aln to-paf mg1655-sakai.1aln \
        --ref-seq mg1655.fa.gz --query-seq sakai.fa.gz -o out.paf
    ```
2.  **Include the cg:Z CIGAR tag**:
    ```bash
    pgr 1aln to-paf mg1655-sakai.1aln \
        --ref-seq mg1655.fa.gz --query-seq sakai.fa.gz --cigar -o out.paf
    ```

---

## to-psl

Expands each alignment record back into base-level aligned columns and emits a
PSL record per alignment.

```bash
pgr 1aln to-psl [OPTIONS] <infile>
```

### Options

*   `--ref-seq <file>`: Reference (a side) genome FASTA / gzipped FASTA.
*   `--query-seq <file>`: Query (b side) genome FASTA / gzipped FASTA.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### Notes

*   Requires `--ref-seq` and `--query-seq` (the source genomes).
*   A `+-` strand means the target (`b`) sequence was stored
    reverse-complemented; the emitted target coordinates are in forward
    orientation.
*   Reads a single `.1aln` file; does not support gzip or stdin.

### Examples

1.  **Convert**:
    ```bash
    pgr 1aln to-psl mg1655-sakai.1aln \
        --ref-seq mg1655.fa.gz --query-seq sakai.fa.gz -o out.psl
    ```