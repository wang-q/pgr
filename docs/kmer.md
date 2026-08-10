# pgr kmer

K-mer counting, per-sequence profiling, and frequency histograms. `pgr kmer`
is the general-purpose k-mer analysis entry point; `pgr rept s-kmer` /
`e-kmer` build on the same library for repeat detection.

All sequence inputs accept FASTA or FASTQ (plain or gzipped) and read from
stdin when the input file is `stdin`.

## Files

Three self-contained single-file formats are produced:

| File | Content |
|------|---------|
| `.pkt` | canonical k-mer count table (sorted, deduplicated, compact encoding) |
| `.pkp` | per-sequence k-mer profiles (raw `u16` counts per position) |
| `.hist` | k-mer frequency histogram in the FastK `.hist` binary layout |
| `.kgc` | GC-content × coverage matrix in the KatGC `.kgc` format |

`.pkt` and `.pkp` are pgr-native formats. `.hist` is byte-compatible with
FastK, so external tooling (Histex, KatGC, GenomeScope) can read it directly.

## Commands

### table

Builds a canonical k-mer count table from one or more sequence files.
Counts accumulate across all inputs; every k-mer is kept, including
singletons (FastK `-t1` semantics).

| Argument | Description |
|----------|-------------|
| `infiles` | Input FASTA/FASTQ file(s) to process |
| `-k`, `--kmer` | K-mer size (default: 17) |
| `-o`, `--outfile` | Output `.pkt` filename |

```bash
pgr kmer table reads.fq.gz -k 21 -o reads.pkt
pgr kmer table a.fa b.fa.gz -k 17 -o all.pkt
```

### hist

Builds a k-mer frequency histogram from sequences or an existing table and
writes it in the FastK `.hist` binary layout (readable by Histex, KatGC, and
GenomeScope tooling). With `--khist-text` and `--peaks`, it also writes the
BBTools `kmercountexact.sh` text outputs (`#Depth Count logScale` histogram
and the peaks summary), byte-identical to kmercountexact with default
parameters.

```bash
pgr kmer hist [OPTIONS] --outfile <outfile> [infile]
```

| Argument | Description |
|----------|-------------|
| `infile` | Input FASTA/FASTQ file to process (unless `--table` is given) |
| `-t, --table <file>` | Reuse a k-mer table (`.pkt`); k is read from the table |
| `-k, --kmer <int>` | K-mer size (required unless `--table` is given) |
| `-o, --outfile <file>` | Output `.hist` filename |
| `--khist-text <file>` | Also write the kmercountexact text histogram |
| `--peaks <file>` | Also write the kmercountexact peaks summary |

```bash
pgr kmer hist reads.fq.gz -k 31 -o reads.hist \
    --khist-text reads.khist.txt --peaks reads.peaks.txt
```


### profile

Generates one k-mer count profile per sequence (read or chromosome) and
writes them to a `.pkp` file. For every k-mer position of every input
sequence, the profile records one count; the counts are looked up from a
k-mer table, either built on the fly from the input itself or reused via
`--table`:

* Without `--table` (self): the input sequences are counted first, and each
  position reports how many times its k-mer occurs in the input dataset.
  Repeated regions therefore show high values (FastK `-p` semantics).
* With `--table` (relative): each position reports the count of its k-mer
  *in the given table*; positions whose k-mer is absent from the table
  report 0 (FastK `-p:<table>` semantics). This is a lookup against an
  external table — not a comparison between profiles.

Both modes write the same `.pkp` format; only the source of the counts
differs. `k` is read from the table unless `--kmer` is given explicitly.

| Argument | Description |
|----------|-------------|
| `infile` | Input FASTA/FASTQ file to process |
| `-t`, `--table` | Reuse a k-mer table (`.pkt`); k is read from the table |
| `-k`, `--kmer` | K-mer size (required unless `--table` is given) |
| `-o`, `--outfile` | Output `.pkp` filename |

```bash
pgr kmer profile genome.fa -k 17 -o genome.pkp
pgr kmer profile reads.fq.gz -t lib.pkt -o reads.pkp
```

### hist

Builds a k-mer frequency histogram and writes it in the FastK `.hist`
binary layout (fixed bins `1..=32767`; counts above the top bin are folded
into it, matching FastK semantics).

Give either a sequence file (histogram computed on the fly) or `--table` to
reuse an existing `.pkt` table.

| Argument | Description |
|----------|-------------|
| `infile` | Input FASTA/FASTQ file (unless `--table` is given) |
| `-t`, `--table` | Reuse a k-mer table (`.pkt`); k is read from the table |
| `-k`, `--kmer` | K-mer size (required unless `--table` is given) |
| `-o`, `--outfile` | Output `.hist` filename |

```bash
pgr kmer hist reads.fq.gz -k 21 -o reads.hist
pgr kmer hist -t reads.pkt -o reads.hist
```

The `.hist` output can be fed directly to external FastK tooling, e.g.:

```bash
Histex -G reads.hist            # GenomeScope ASCII format
Histex -A reads.hist            # tab-separated ASCII histogram
```

### gc

Builds the two-dimensional GC-content × k-mer coverage matrix and writes it
in the KatGC `.kgc` format (`GCP KF Count` rows, 2×2 neighbor average,
values clamped to the peak). Rows are GC counts, columns are count bins, so
the matrix shows how k-mer coverage varies with GC content — a typical
quality diagnostic for read sets.

Give either a sequence file (table built on the fly) or `--table` to reuse
an existing `.pkt` table. Counts above the count cap are folded into the top
bin; the output x-range is the peak coverage times `--xrel` (default 2.1),
unless `--xmax` pins it absolutely (which also sets the count cap, as in
KatGC).

| Argument | Description |
|----------|-------------|
| `infile` | Input FASTA/FASTQ file (unless `--table` is given) |
| `-t`, `--table` | Reuse a k-mer table (`.pkt`); k is read from the table |
| `-k`, `--kmer` | K-mer size (required unless `--table` is given) |
| `-X`, `--xmax` | Absolute x max (also caps the count axis; default: auto) |
| `-x`, `--xrel` | Max x as a multiple of the peak coverage (default: 2.1) |
| `--tex` | Render a LaTeX heatmap instead of the `.kgc` matrix |
| `-o`, `--outfile` | Output `.kgc` filename |

```bash
pgr kmer gc reads.fq.gz -k 21 -o reads.kgc
pgr kmer gc -t reads.pkt -x 1.9 -o reads.kgc
pgr kmer gc reads.fq.gz -k 21 --tex -o reads.tex
```

The matrix is byte-identical to MerquryFK KatGC output on the same input
(verified against a locally compiled KatGC). `--tex` renders the KatGC
heat-map equivalent as a LaTeX figure (pgfplots; compile with tectonic),
via the shared `pgr plot heat` renderer (`pgr plot heat reads.kgc -o reads.tex`
is equivalent).

### qhist

Builds the quality-weighted k-mer frequency histogram from FASTQ reads and
writes it in quorum's `histo_mer_database` format: `count n_lowq n_highq`
per non-empty count bin.

A k-mer counts as high quality iff all `k` bases of its window score at
least the quality threshold; per k-mer the final count is the number of
high-quality occurrences when any exist, otherwise the number of
low-quality occurrences (quorum `hash_with_quality` semantics:
low-quality evidence never raises a high-quality count). Counts are capped
at 1000, and additionally by `--bits` (default 7, quorum's
`create_database` default: max count 127).

The threshold defaults to the detected Phred offset (+33/+64) plus 5
(quorum's default min-quality offset) and can be pinned with
`--qual-thresh`.

| Argument | Description |
|----------|-------------|
| `infile` | Input FASTQ file to process (FASTA is rejected) |
| `-k`, `--kmer` | K-mer size (default: 17) |
| `-q`, `--qual-thresh` | Quality ASCII threshold (default: detected Phred offset + 5) |
| `-b`, `--bits` | Count bits (quorum `create_database -b`; max count = 2^bits − 1) |
| `-o`, `--outfile` | Output histogram filename |

```bash
pgr kmer qhist reads.fq.gz -k 21 -o reads.qhist
pgr kmer qhist reads.fq.gz -k 21 -q 43 -o reads.qhist   # Phred+33 Q10
```

The output format matches quorum's `histo_mer_database`, so it can be
compared against quorum/Jellyfish pipelines directly.

### qcheck

Flags reads that quorum would correct or truncate and keeps the rest
untouched. A quality-weighted k-mer table is built from the input reads
first, then each read is checked for quorum's error signals: no high-quality
anchor, a k-mer with no continuation (truncation), or a base that quorum
would substitute (including the Poisson collision test). No corrected
sequence is produced — the read is kept as-is or discarded.

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

```bash
pgr kmer qcheck reads.fq.gz -k 21 -o kept.fq.gz
pgr kmer qcheck reads.fq.gz -k 21 -o kept.fq --discard-file bad.fq
```

The error signals mirror quorum's `find_starting_mer` / `extend` semantics
(high-quality anchors only, substitution/truncation events), so `qcheck`
reproduces which reads quorum would touch without producing corrected
sequences.

### gsize

Estimates the k-mer coverage peak and genome size from a count table.
`peak_coverage` is the frequency carried by the most distinct k-mers (the
main mode); `genome_size` is total k-mer instances / peak coverage — the
cheap haploid estimate that precedes GenomeScope-style model fitting.

Give either a sequence file (table built on the fly) or `--table` to reuse
an existing `.pkt` table.

| Argument | Description |
|----------|-------------|
| `infile` | Input FASTA/FASTQ file (unless `--table` is given) |
| `-t`, `--table` | Reuse a k-mer table (`.pkt`); k is read from the table |
| `-k`, `--kmer` | K-mer size (required unless `--table` is given) |
| `--model` | Fit the GenomeScope model (kmercov/het/genome size) |
| `--plot` | With `--model`, also write `spectra.tex` to the output directory |
| `-p`, `--ploidy` | Ploidy for the model (1 or 2; default 1) |
| `-o`, `--outfile` | Output statistics (default: stdout) |

```bash
pgr kmer gsize reads.fq.gz -k 21
pgr kmer gsize reads.fq.gz -k 21 --model -o gs_out
pgr kmer gsize -t reads.pkt -o stats.tsv
```

Output rows are tab-separated: `k`, `peak_coverage`, `total_distinct`,
`total_kmers`, `genome_size`. On a synthetic 30× 1 kb test the estimate
lands within ~3% of the true size.

With `--model`, a native port of `genescopefk.R` (GenomeScope 2.0) is run:
a negative-binomial mixture (unique + repeat components for p=1; the four
AA/AB/BB classes for p=2) is fitted by Levenberg-Marquardt least squares
across four trimming rounds, with model scoring and the best-round
selection. `-o` is then an output directory holding `summary.txt` and
`model.txt` in the GenomeScope formats consumed by anchr's `2_fastk`
(`grep '^kmercov' model.txt`), and the summary is printed to stdout.
With `--plot`, the spectra figure is also written to
`spectra.tex` in that directory (equivalent to `pgr plot spectra`).

The fit recovers known parameters exactly on noiseless synthetic spectra
(p=1: kmercov, bias, d, length all within rounding of the truth). Real
reads with coverage heterogeneity (read ends) are handled by the
negative-binomial dispersion, though the simplified port (no
polyploidy>2, no error component) can still overestimate the genome size;
the `kmercov` estimate stays close to the true coverage.
