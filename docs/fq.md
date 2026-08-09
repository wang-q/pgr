# pgr fq

`pgr fq` provides tools for manipulating **FASTQ** files.

## Subcommands

*   `interleave`: Interleave paired-end sequences from one or two files.
*   `to-fa`: Convert FASTQ files to FASTA format.

---

## interleave (il)

Interleaves paired-end sequences from one or two files.

It can combine separate R1 and R2 files into a single interleaved file, or generate dummy R2 sequences (N's) from a single R1 file.

```bash
pgr fq interleave [OPTIONS] <infiles>...
```

### Options

*   `--fq`: Write output in FASTQ format (default is FASTA).
    *   For FASTQ output, quality scores are preserved from input FASTQ files.
    *   If input is FASTA, quality scores are set to '!' (ASCII 33).
*   `--name-prefix <string>`: Prefix for record names (default: "read").
*   `--start-index <int>`: Starting index for record numbering (default: 0).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### Examples

1.  **Interleave two FASTQ files into one**:
    ```bash
    pgr fq interleave R1.fq R2.fq -o interleaved.fq --fq
    ```

2.  **Generate dummy pairs from a single FASTA file**:
    ```bash
    pgr fq interleave R1.fa --name-prefix sample --start-index 1
    ```

3.  **Convert separate FASTA files to interleaved FASTQ**:
    ```bash
    pgr fq interleave R1.fa R2.fa --fq -o out.fq
    ```

---

## to-fa

Converts FASTQ files to FASTA format.

This command preserves sequence names and supports multiple input files.

```bash
pgr fq to-fa [OPTIONS] <infiles>...
```

### Options

*   `-o, --outfile <file>`: Output filename (default: stdout).

### Examples

1.  **Convert a FASTQ file to FASTA**:
    ```bash
    pgr fq to-fa input.fq -o output.fa
    ```

2.  **Convert multiple FASTQ files to a single FASTA**:
    ```bash
    pgr fq to-fa input1.fq input2.fq -o output.fa
    ```

## range

Extracts FASTQ records by read name (or a region within a read) using a
`.loc` index that is created automatically.

### Options

| Argument | Description |
|----------|-------------|
| `infile` | Input FASTQ file (plain text or BGZF `.gz`) |
| `--mate <FILE>` | Second mate file for paired-end extraction |
| `ranges` | Read names and/or `name:start-end` regions (or `-r` list file) |
| `-r`, `--rgfile <FILE>` | Read names/regions from a file |
| `-c`, `--cache <N>` | LRU cache capacity for extracted records (default 1) |
| `-o`, `--outfile <FILE>` | Output filename (default: stdout) |
| `--outfile-2 <FILE>` | Output for the second mate (requires `--mate`) |
| `-u`, `--update` | Force rebuild the `.loc` index |

Read names with `/1` `/2` suffixes are matched by their pair name;
interleaved reads with identical names are both returned in order.
`name:start-end` returns the subsequence of both sequence and quality
(1-based inclusive). The index is rebuilt when the input is newer.

### Examples

```bash
pgr fq range reads.fq read1 read2 -o out.fq
pgr fq range reads.fq "read1:10-100" -o out.fq
pgr fq range reads.fq -r names.txt -c 10 -o out.fq
pgr fq range R1.fq --mate R2.fq read1 -o r1.out.fq --outfile-2 r2.out.fq
```
