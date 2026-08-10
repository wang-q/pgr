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

## clump

Sorts paired reads by the pivot k-mer of R1, reproducing the BBTools
`clumpify.sh` default output order byte for byte. Reads sharing k-mers end up
adjacent, which speeds up the k-mer steps that follow in a read-cleaning
pipeline.

```bash
pgr fq clump [OPTIONS] <infiles>...
```

### Options

*   `-o, --outfile <file>`: Output filename (default: stdout).
*   `-k, --kmer <int>`: K-mer size (default: 31).
*   `--seed <int>`: Comparator seed (default: 1).
*   `--dedupe`: Remove duplicate read pairs (R1 and R2 both exact within
    `--dupesubs`, N wildcard; keeps the higher-quality copy).
*   `--dupesubs <int>`: Maximum substitutions allowed in a duplicate
    (default: 0).
*   `--mem <size>`: In-memory sort budget (KMG, default 2g). Data estimated
    to exceed `min(--mem, physical/2)` is sorted via external hash buckets
    (deterministic, bucket-concatenated order).
*   `--buckets <int>`: External-path hash bucket count (default: derived from
    the memory budget).
*   `--sort-mode <auto|global|bucket>`: Force the sorting path. `auto`
    (default) picks by the memory budget; `global` always sorts in memory;
    `bucket` always uses the external hash-bucket path (implied when
    `--buckets` is given).

### Examples

1.  **Sort two paired FASTQ files**:
    ```bash
    pgr fq clump R1.fq.gz R2.fq.gz -o clumped.fq
    ```

2.  **Sort an interleaved file with a different seed**:
    ```bash
    pgr fq clump in.fq.gz -o out.fq --seed 2
    ```

3.  **Remove exact duplicate pairs**:
    ```bash
    pgr fq clump R1.fq.gz R2.fq.gz -o out.fq --dedupe --dupesubs 0
    ```

4.  **Bound memory to 1 GiB** (larger data goes through external buckets):
    ```bash
    pgr fq clump R1.fq.gz R2.fq.gz -o out.fq --mem 1g
    ```

5.  **Force the external bucket path**:
    ```bash
    pgr fq clump R1.fq.gz R2.fq.gz -o out.fq --sort-mode bucket
    ```

---

## to-fa

Converts FASTQ files to FASTA format.

This command preserves sequence names and supports multiple input files.

```bash
pgr fq to-fa [OPTIONS] <infiles>...
```

---

## split

Splits an interleaved FASTQ file into paired-end R1/R2 outputs and a singles
file. It is the inverse of `pgr fq interleave` and matches BBTools
`repair.sh` in `rp` mode.

```bash
pgr fq split [OPTIONS] <infile>
```

### Options

*   `-o, --outfile <file>`: R1 output filename (default: stdout).
*   `--outfile-2 <file>`: R2 output filename (required).
*   `--outfile-single <file>`: Output file for unpaired reads (optional).

### Examples

1.  **Split into R1/R2 and singles**:
    ```bash
    pgr fq split interleaved.fq -o r1.fq --outfile-2 r2.fq --outfile-single s.fq
    ```

2.  **Split a gzipped file into R1/R2 only**:
    ```bash
    pgr fq split interleaved.fq.gz -o r1.fq --outfile-2 r2.fq
    ```

---

## sample

Downsamples reads so the output contains approximately the requested number of
bases, preserving input order. Selection is deterministic for a given seed and
matches BBTools `reformat.sh` with `samplebasestarget` and a fixed
`sampleseed`.

```bash
pgr fq sample [OPTIONS] <infile>
```

### Options

*   `-o, --outfile <file>`: Output filename (default: stdout).
*   `--bases <int>`: Target number of output bases.
*   `--seed <int>`: Random seed for deterministic selection (default: 1).

### Examples

1.  **Keep about 1 million bases**:
    ```bash
    pgr fq sample reads.fq -o out.fq --bases 1000000
    ```

2.  **Reproduce a BBTools run with a fixed seed**:
    ```bash
    pgr fq sample reads.fq.gz -o out.fq --bases 1000000 --seed 42
    ```

---

## trim-adapter

Removes adapter/contaminant sequences by matching read k-mers against a
reference, then quality-trims and length-filters the reads. It reproduces
BBTools 39.38 `bbduk.sh` output byte for byte for the anchr trim pipeline
parameters (deterministic `ordered=t` mode).

```bash
pgr fq trim-adapter [OPTIONS] <infiles>...
```

### Options

*   `--ref <file>`: Reference FASTA of adapters/contaminants (required).
*   `-k, --k <int>`: K-mer size (default: 23).
*   `--mink <int>`: Minimum short k-mer size at read ends (default: 11).
*   `--hdist <int>`: Reference hamming distance (default: 1).
*   `--no-ktrim`: Disable k-mer trimming (filtering mode).
*   `--no-tbo`: Disable overlap trimming.
*   `--no-tpe`: Disable even pair trimming.
*   `--no-qtrim`: Disable quality trimming.
*   `--trimq <int>`: Quality threshold for `qtrim=r` (default: 15).
*   `--minlen <int>`: Minimum kept read length (default: 60).
*   `--maxns <int>`: Maximum allowed N bases; negative disables (default: 0).
*   `--ftm <int>`: Right-trim lengths to a multiple (default: 5).
*   `--no-toss-broken-reads`: Keep surviving mates of discarded reads.
*   `-p, --threads <int|auto>`: Worker threads (default: logical CPU count);
    output order is preserved for any thread count.

### Examples

1.  **Adapter trim with the anchr pipeline defaults**:
    ```bash
    pgr fq trim-adapter R1.fq.gz R2.fq.gz --ref illumina_adapters.fa -o out.fq
    ```

2.  **K-mer filtering mode (bbduk filter step)**:
    ```bash
    pgr fq trim-adapter in.fq --ref illumina_adapters.fa --no-ktrim \
        --no-tbo --no-tpe --no-qtrim --k 27 --mink 0 --minlen 0 \
        --maxns=-1 --ftm 0 -o out.fq
    ```

3.  **Run with a specific number of threads**:
    ```bash
    pgr fq trim-adapter R1.fq.gz R2.fq.gz --ref illumina_adapters.fa \
        -o out.fq --threads 8
    ```

---

## norm

Removes reads whose k-mer coverage is below a minimum depth, following the
BBTools 39.38 `bbnorm.sh passes=1 bits=16 min=<n> target=9999999` read
decision logic. Counts are exact (canonical table); bbnorm's `bits=16`
approximate hash counts can differ on reads near the depth boundary.

```bash
pgr fq norm [OPTIONS] <infiles>...
```

### Options

*   `-k, --kmer <int>`: K-mer size (default: 31).
*   `--min <int>`: Minimum k-mer depth cutoff (default: 3).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### Examples

1.  **Keep reads with at least one k-mer at depth 3**:
    ```bash
    pgr fq norm reads.fq.gz -k 31 --min 3 -o out.fq
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
