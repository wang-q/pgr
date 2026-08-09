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
