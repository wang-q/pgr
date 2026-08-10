# pgr fq

`pgr fq` provides tools for manipulating **FASTQ** files.

## Subcommands

*   `interleave` (`il`): Interleave paired-end sequences from one or two files.
*   `to-fa`: Convert FASTQ files to FASTA format.
*   `clump`: Sort reads by k-mer signature (clumpify-compatible).
*   `split`: Split interleaved FASTQ into R1/R2/singles files.
*   `sample`: Subsample reads to a target base count.
*   `clean`: Adapter k-mer trimming, quality and composition filtering (bbduk).
*   `filter`: Discard reads matching reference k-mers (bbduk kfilter).
*   `norm`: Filter reads by k-mer depth (bbnorm-style cutoff).
*   `range`: Extract FASTQ records by name or region.
*   `trim-qual`: Trim reads by quality score (sickle/cutadapt-style).

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
*   `-p, --parallel <int|auto>`: Worker threads for the parallel sort and
    bucket processing (default: logical CPU count).

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

6.  **Limit parallelism**:
    ```bash
    pgr fq clump R1.fq.gz R2.fq.gz -o out.fq --parallel 4
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

---

## split

Splits an interleaved FASTQ file into paired-end R1/R2 outputs and a singles
file. It is the inverse of `pgr fq interleave` and matches BBTools
`repair.sh` in `rp` mode. By default reads are paired by position (every two
records); `--repair` instead matches mates by read-name prefix, recovering
disordered pairs and routing orphaned reads to singles.

```bash
pgr fq split [OPTIONS] <infile>
```

### Options

*   `-o, --outfile <file>`: R1 output filename (default: stdout).
*   `--outfile-2 <file>`: R2 output filename (required).
*   `--outfile-single <file>`: Output file for unpaired reads (optional).
*   `--repair`: Pair mates by read-name prefix (`/1` `/2` or `1:` `2:`
    suffixes, with the repair.sh fallback) instead of position; buffers
    unpaired reads in memory like `repair.sh`.

### Examples

1.  **Split into R1/R2 and singles**:
    ```bash
    pgr fq split interleaved.fq -o r1.fq --outfile-2 r2.fq --outfile-single s.fq
    ```

2.  **Split a gzipped file into R1/R2 only**:
    ```bash
    pgr fq split interleaved.fq.gz -o r1.fq --outfile-2 r2.fq
    ```

3.  **Repair disordered pairs and orphaned reads**:
    ```bash
    pgr fq split disordered.fq.gz --repair -o r1.fq --outfile-2 r2.fq \
        --outfile-single s.fq
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

## clean

Cleans reads in one pass: adapter/contaminant k-mer trimming, quality
trimming, polymer and GC filtering, and masking. Reproduces the first
BBTools 39.38 `bbduk.sh` call of the anchr trim pipeline (the `ktrim` pass)
byte for byte (`ordered=t`, deterministic). For the second bbduk call
(k-mer contaminant filtering) use `pgr fq filter`; for sickle-style pure
quality trimming use `pgr fq trim-qual`.

```bash
pgr fq clean [OPTIONS] <infiles>...
```

### Options

Options identical to their bbduk counterparts are not annotated; renamed
options show the bbduk name in parentheses.

*   `--ref <file>`: Reference FASTA of adapters/contaminants. Omit it to
    skip all k-mer operations and only quality-trim/filter.
*   `-k, --k <int>`: K-mer size (default 23).
*   `--min-k <int>`: Minimum short k-mer size at read ends;
    default 11).
*   `--hamming-distance <int>`: Reference hamming distance (bbduk: `hdist`;
    default 1).
*   `--no-trim-by-overlap`: Disable overlap trimming (bbduk: `tbo=f`).
*   `--no-trim-pair-evenly`: Disable even pair trimming (bbduk: `tpe=f`).
*   `--no-qtrim`: Disable quality trimming.
*   `--qtrim <r|l|rl|w|f>`: Quality trim mode (default `r`;
    `w` uses a sliding window of `--qtrim-window`).
*   `--qtrim-window <int>`: Window size for `--qtrim w` (bbduk: `qtrim=w,N`;
    default 4).
*   `--trim-quality <int>`: Quality threshold for qtrim (bbduk: `trimq`;
    default 15).
*   `--minlen <int>`: Minimum kept read length (default 60).
*   `--minlen-fraction <float>`: Minimum read length as a fraction of the
    original (bbduk: `mlf`; default 0).
*   `--max-ns <int>`: Maximum allowed N bases; negative disables
    (default 0).
*   `--max-n-rate <float>`: Discard reads with more than this fraction of Ns
    (default 1, disabled).
*   `--force-trim-mod <int>`: Right-trim lengths to a multiple (bbduk: `ftm`;
    default 5).
*   `--force-trim-left <int>`: Trim bases left of this position.
*   `--force-trim-right <int>`: Trim bases right of this position.
*   `--force-trim-right2 <int>`: Trim this many bases on the right end.
*   `--trim-poly-a <int>`: Trim poly-A/T tails.
*   `--trim-poly-g-left/--trim-poly-g-right <int>`: Trim poly-G
    prefixes/tails.
*   `--filter-poly-g <int>`: Discard reads with a poly-G prefix; poly-C
    equivalents: `--trim-poly-c-left/right`, `--filter-poly-c`.
*   `--max-non-poly <int>`: Allowed non-polymer bases in a polymer run
    (default 1).
*   `--min-avg-quality <float>`: Discard reads with average quality below
    this (bbduk: `maq`).
*   `--min-avg-quality-bases <int>`: Use only this many leading bases for
    `--min-avg-quality` (bbduk: `maqb`).
*   `--min-base-quality <int>`: Discard reads with any base below this
    quality (bbduk: `mbq`).
*   `--min-consecutive-bases <int>`: Discard reads without this many
    consecutive ACGT bases (bbduk: `mcb`).
*   `--maxlength <int>`: Discard reads longer than this (default 0 = off).
*   `--min-gc <float>` / `--max-gc <float>`: Discard reads with GC content
    below/above these bounds.
*   `--no-pair-gc`: Check GC per read instead of the pair average (bbduk:
    `gcpairs=f`).
*   `--mask-kmers <symbol|lc|t>`: Mask matching k-mers with a symbol
    (`t` = `N`) or lowercase them (`lc`) instead of trimming (bbduk:
    `kmask`); requires `--ref`.
*   `--mask-fully-covered`: Only mask bases fully covered by matching k-mers.
*   `--trim-pad <int>`: Extra bases to mask around matching k-mers.
*   `--no-toss-broken-reads`: Keep surviving mates of discarded reads (bbduk:
    `removeifeitherbad=f`).
*   `-p, --parallel <int|auto>`: Worker threads (bbduk: `threads`; default
    logical CPU count); output order is preserved for any thread count.
*   `--stats <file>`: Write per-reference match statistics in the bbduk
    `stats=` format (tab-separated, `#Matched`/`#Name` header lines).

### Examples

1.  **Clean with the anchr pipeline defaults**:
    ```bash
    pgr fq clean R1.fq.gz R2.fq.gz --ref illumina_adapters.fa -o out.fq
    ```

2.  **Clean without a reference (bbduk `qtrim=r minlen=...` style)**:
    ```bash
    pgr fq clean unmerged.raw.fq.gz -o unmerged.trim.fq \
        --no-trim-by-overlap --no-trim-pair-evenly \
        --max-ns=-1 --force-trim-mod 0 --trim-quality 25 --minlen 60
    ```

3.  **Mask matching k-mers instead of trimming**:
    ```bash
    pgr fq clean in.fq --ref illumina_adapters.fa --mask-kmers N -o out.fq
    ```

---

## filter

Discards reads containing k-mers matching a reference (adapters,
contaminants, spike-ins). Reproduces the second BBTools 39.38 `bbduk.sh`
call of the anchr trim pipeline (`k=<matchk> cardinality`) byte for byte
(`ordered=t`, deterministic). A read is discarded when more than zero
k-mers match (bbduk `minkmerhits=1`); surviving mates follow
`--toss-broken-reads`.

```bash
pgr fq filter [OPTIONS] <infiles>...
```

### Options

*   `--ref <file>`: Reference FASTA of contaminants/adapters (required).
*   `-k, --k <int>`: K-mer size (default 27, the anchr `matchk`).
*   `--min-k <int>`: Minimum short k-mer size at read ends;
    default 0).
*   `--hamming-distance <int>`: Reference hamming distance (bbduk: `hdist`;
    default 0).
*   `--minlen <int>`: Minimum kept read length (default 10).
*   `--max-ns <int>`: Maximum allowed N bases; negative disables
    (default -1).
*   `--no-toss-broken-reads`: Keep surviving mates of discarded reads (bbduk:
    `removeifeitherbad=f`).
*   `-p, --parallel <int|auto>`: Worker threads (bbduk: `threads`).
*   `--stats <file>`: Write per-reference match statistics in the bbduk
    `stats=` format.

### Examples

1.  **Filter adapter/artifact matches (anchr filter step)**:
    ```bash
    pgr fq filter trim.fq.gz --ref illumina_adapters.fa -k 27 -o filter.fq
    ```

---

## s-filter

Discards reads whose own k-mer counts look erroneous, using the input reads
themselves as the reference (quorum's error-correction signals, no external
contaminant database needed). The `s-` prefix marks this as a self/internal
check, in contrast to `pgr fq filter`, which matches reads against an
external reference. A quality-weighted k-mer table is built from the input
first, then each read is checked for quorum's signals: no high-quality
anchor, a k-mer with no continuation (truncation), or a base that quorum
would substitute (including the Poisson collision test). No corrected
sequence is produced — the read is kept as-is or discarded.

```bash
pgr fq s-filter [OPTIONS] <infile>
```

### Options

| Argument | Description |
|----------|-------------|
| `infile` | Input FASTQ file to process (FASTA is rejected) |
| `-k`, `--kmer` | K-mer size (default: 17) |
| `-q`, `--qual-thresh` | Table quality threshold (default: detected Phred offset + 5) |
| `-b`, `--bits` | Table count bits (default: 7, max count 127) |
| `--skip` / `--good` / `--anchor-count` | Anchor search parameters (defaults 0 / 1 / 1) |
| `--min-count` / `--cutoff` | Trusted-count parameters (defaults 1 / 4) |
| `--apriori-error-rate` / `--poisson-threshold` | Poisson collision test (defaults 0.01 / 1e-6) |
| `-o`, `--outfile` | Output FASTQ of kept reads |
| `--discard-file` | Optional FASTQ of flagged reads |

### Examples

1.  **Self-check and discard error-prone reads**:
    ```bash
    pgr fq s-filter reads.fq.gz -k 21 -o kept.fq.gz
    ```

2.  **Also write the discarded reads**:
    ```bash
    pgr fq s-filter reads.fq.gz -k 21 -o kept.fq --discard-file bad.fq
    ```

The error signals mirror quorum's `find_starting_mer` / `extend` semantics
(high-quality anchors only, substitution/truncation events), so `s-filter`
reproduces which reads quorum would touch without producing corrected
sequences. End-to-end agreement with `quorum_error_correct_reads` is verified
on the Lambda golden data.

---

## norm

Removes reads whose k-mer coverage is below a minimum depth, following the
BBTools 39.38 `bbnorm.sh passes=1 bits=16 min=<n> target=9999999` read
decision logic. Counts are exact (canonical table); bbnorm's `bits=16`
approximate hash counts can differ on reads near the depth boundary.
The bbnorm defaults `changequality` (N bases get quality 0, ACGT bases a
minimum of 2) and `minq=6` (k-mers containing bases below quality 6 are
excluded from the count table) are applied; bbnorm's nominal `minprob=0.5`
is not, matching the KmerCount table bbnorm actually uses for `bits=16`.

```bash
pgr fq norm [OPTIONS] <infiles>...
```

### Options

*   `-k, --kmer <int>`: K-mer size (default: 31).
*   `--min <int>`: Minimum k-mer depth cutoff (default: 3).
*   `-p, --parallel <int>`: Worker threads (default: logical CPU count).
*   `--mem <size>`: In-memory count budget, KMG units (default: 2g). Data
    estimated to exceed it is counted via external hash buckets.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### Examples

1.  **Keep reads with at least one k-mer at depth 3**:
    ```bash
    pgr fq norm reads.fq.gz -k 31 --min 3 -o out.fq
    ```

---

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

---

## trim-qual

Trims low-quality bases from read ends using a sliding window (default) or the
Mott cumulative-quality algorithm. Quality trimming only; adapters are not
removed (use `pgr fq clean`/`pgr fq filter` for adapter/contaminant trimming).

```bash
pgr fq trim-qual [OPTIONS] <infiles>...
```

### Options

*   `-o, --outfile <file>`: Output filename (default: stdout).
*   `--outfile-2 <file>`: R2 output file (paired-end; omit for interleaved
    output).
*   `--outfile-single <file>`: Output file for surviving single-end reads.
*   `-q, --qual-threshold <float>`: Quality threshold (default: 20).
*   `-l, --length-threshold <int>`: Minimum kept length; shorter reads are
    discarded (default: 20).
*   `--method <sliding|mott>`: Trimming algorithm (default: sliding).
*   `--no-fiveprime`: Disable 5' trimming.
*   `--quality-base <33|64|auto>`: Input quality encoding (default: auto).
*   `--polyg-right <int>`: Trim 3' poly-G tails of at least this length
    (0 disables).

### Examples

1.  **Single-end trimming**:
    ```bash
    pgr fq trim-qual in.fq -o out.fq
    ```

2.  **Paired-end with separate outputs and singles**:
    ```bash
    pgr fq trim-qual R1.fq R2.fq -o r1.fq --outfile-2 r2.fq --outfile-single s.fq
    ```

3.  **Paired-end interleaved output**:
    ```bash
    pgr fq trim-qual R1.fq R2.fq -o interleaved.fq
    ```
