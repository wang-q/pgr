# pgr pgi

`pgr pgi` manages **pgr genome index (.pgi)** files. A `.pgi` is a binary,
syncmer-sparse, sorted k-mer index of a genome: it stores the 2-bit encoded
k-mers that start at closed-syncmer positions, with per-k-mer positions
(contig, offset, strand). The sorted key stream supports linear two-index
merges for distance computation and seed discovery.

## Core positioning

- **Purpose**: build and consume a genome k-mer index.
- **Input**: FASTA (plain or `.gz`) or 2bit genome files.
- **Output**: `.pgi` binary indexes, `.hv` hypervectors, or PSL alignments.
- **Complements**:
  - Upstream: `pgr fa to-2bit` (fastest index input), `pgr dist mini/mash/frac/hv`
    (sketch distances from sequences).
  - Downstream: `pgr dist pgi` (merge distance on the index; secondary
    consumer — see note below), `pgr psl to-chain` and `pgr pl chainnet`
    (chain the PSL blocks from `pgr align pgi`).
- **Storage note (2026-08-08)**: a `.pgi` is ~27× the gzipped FASTA
  (e.g. MG1655: 37.6 MB index vs 1.4 MB FASTA; `.hv` is 16 KB). The index
  exists for `align pgi` (chaining needs positions/strands); **distance
  computation should not build `.pgi`** — use `dist mini`/`dist frac` (FASTA, no index)
  or `dist hv` (`.hv`, 1/87 of FASTA). `dist pgi` is a secondary consumer
  for the "index already built" case and as a calibration reference for
  `.hv`; its merge distance has syncmer sampling bias (see
  `notes/benchmarks/dist-cohort-validation.md`).
- **Design notes**: `notes/design/pbit.md` (index consumers)
  and `notes/design/pgi-align.md` (alignment pipeline).

## Index format

Each `.pgi` file (format v2, magic `PGI1`) stores, in order:

1. Sampling parameters: `k` (k-mer size), syncmer `smer`/`window`,
   record field widths;
2. The contig table (name, length);
3. A sorted per-occurrence record stream: each record is the packed 2-bit
   k-mer (`ceil(k/4)` bytes, big-endian) + the position in minimal
   little-endian bytes + a packed `contig_id | (strand << high_bit)` byte.

The packed layout follows FastGA's GIX (see
`notes/benchmarks/bench-pgi-vs-gix-storage.md`): k-mers are stored at their
true 2-bit width (10 bytes at k=40) and positions use only the bytes
needed for the largest contig, with the strand flag folded into the contig
field. Repeated k-mers are stored per occurrence on disk; readers re-group
them into unique entries + position lists in memory.

Both strands are indexed by default: each syncmer position emits the forward
k-mer and its reverse complement (strand flag 0/1). Sampling parameters must
match for any two-index operation (`dist pgi`, `align`).

## Subcommands

### `pgr pgi build`

Builds a `.pgi` index from FASTA or 2bit.

```
pgr pgi build <infile> -o out.pgi [-k 40] [--smer 8] [--window 5] [--no-rev] [--mask]
```

- `-k`/`--kmer`: k-mer size, at most 64 (default 40, matching FastGA GIX);
- `--smer`/`--window`: closed-syncmer sampling parameters (default 8/5,
  span = smer + window - 1);
- `--no-rev`: index the forward strand only;
- `--mask`: skip soft-masked regions (lowercase FASTA bases / 2bit mask
  blocks), matching FastGA `-M` — masked repeats and low-complexity regions
  produce no seeds.

The default parameters (k=40, syncmer 8/5) align with FastGA's GIX so the
index parameter set is comparable against the C tool. Note the sampling rule
differs: pgr uses Edgar's closed syncmer and anchors each seed k-mer at its
window's minimal s-mer endpoint, whereas GIX's match-mer anchors at the window
start. On divergent sequences the two rules can yield different seed
positions, but the resulting alignment is equivalent (verified end-to-end).
2bit input is fastest (no FASTA decoding).

### `pgr pgi stat`

Prints the index parameters, contig count, unique k-mer and position counts,
and the file size.

### `pgr pgi to-hv`

Projects the index's k-mer set onto a fixed-dimension hypervector for fast
pairwise comparison:

```
pgr pgi to-hv in.pgi -o out.hv [--dim 4096] [--sparse 1]
```

The projection is sparse: each k-mer updates `--sparse` random dimensions,
so the shared-k-mer signal stays dominant for large k-mer sets. The output
`.hv` can be compared with `pgr dist hv`, which recovers the k-mer set
overlap via cosine similarity (approximating `pgr dist pgi` at ~50x speed).
Sampling parameters, sparse-update count, and dimension must match for
comparisons; `.hv` files from different parameter sets are not comparable.

### `pgr align pgi`

The pairwise genome alignment consumer lives under `pgr align`; see
[`align-pgi.md`](align-pgi.md). It accepts genome sequences directly
(indexing them automatically, reusing a same-named `.pgi` when present) or
explicit `.pgi` indexes.

## Notes

* 2bit inputs are preferred for speed and random access.
* `.pgi` files are not gzip-compressed.
* Both indexes in a comparison must use identical `-k/--smer/--window`.

## Examples

1. Compute the merge distance between two indexes (see the bias note above):
   ```
   pgr dist pgi a.pgi b.pgi
   ```

> **Bias note (2026-08-08)**: `dist pgi` merges the *syncmer-sampled*
> k-mer sets, not the full k-mer sets. On closely related genomes (or any
> pair rich in repeats/mobile elements) the syncmer positions drift apart,
> so the intersection is underestimated and Jaccard/containment are biased
> low (up to ~3% ANI; containment is slightly more stable than Jaccard,
> ~34% vs ~39% relative error on a 5-strain test). Rankings stay roughly
> usable. For unbiased numeric ANI use `pgr dist frac`
> (with `--ci`); `dist pgi` is for the "index already built" case and as a
> calibration reference for `.hv`. Details:
> `notes/benchmarks/dist-cohort-validation.md`.
