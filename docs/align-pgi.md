# pgr align pgi

`pgr align pgi` aligns two genomes on the pgr k-mer index pipeline, emitting
one PSL block per chain:

```
pgr align pgi ref query -o out.psl
  [-f 10] [-c 85] [-s 1000] [--band 128] [--merge-gap 5000]
  [--min-shared N] [--workflow greedy|tube] [--parallel 8]
  [-k 40] [--smer 8] [--window 5] [--keep-index]
  [--ref-seq ref.fa|ref.2bit] [--query-seq query.fa|query.2bit]
```

`ref` and `query` may be genome sequences (FASTA, gzipped FASTA or .2bit) or
`.pgi` indexes, mixed freely:

- A genome sequence is indexed automatically. An index named like the input
  (e.g. `ref.fa` → `ref.pgi`) is reused when present; otherwise one is built
  in a temporary directory and removed afterwards. `--keep-index` keeps it
  next to the input. The sequence itself is then used to refine the chains.
  For `.gz` inputs the sibling index is `<name-without-.gz>.pgi`
  (e.g. `ref.fa.gz` → `ref.fa.pgi`).
- A `.pgi` index is used directly; `--ref-seq`/`--query-seq` may then supply
  the sequences for chain refinement, and are validated against the index
  contig table.

With a single input (or `--self` with the same input as query) the genome is
aligned to itself (internal repeats and haplotype-level homology, FastGA's
self mode); exact self-identity hits are dropped.

## Options

- `-f`/`--freq`: drop k-mers occurring more than this many times on either
  side (repeats are not seeds);
- `-c`/`--min-span`: minimum per-axis seed span (bp) for a chain;
- `-s`/`--max-gap`: maximum bp gap between consecutive seeds in a chain;
- `--band`: diagonal band half-width (bp) around the chain mean;
- `--merge-gap`: merge adjacent colinear chains separated by at most this gap
  (bp), stitching blocks split by insertions (IS elements);
- `--min-shared`: minimum shared seed length (bp); default = k for greedy,
  12 for tube;
- `--workflow`: chaining workflow, `greedy` (default) or FastGA-style `tube`;
- `--parallel`: rayon thread count (default 8, FastGA `-T` default);
- `-k`/`--smer`/`--window`: sampling parameters for automatic indexing of
  genome inputs only (default 40/8/5, matching FastGA GIX);
- `--keep-index`: write automatically built indexes next to the inputs;
- `--self`: self-alignment; the query input may be omitted or must be the same
  as the reference;
- `--ref-seq`/`--query-seq`: sequences for chain refinement of `.pgi` inputs
  (FASTA or .2bit).

Without extension sequences, each chain is emitted as a single PSL block.
With sequences, chains are refined by a banded local alignment into scored
PSL records with real blocks; chains longer than 16 kb are split into
overlapping windows. The output feeds directly into
`pgr psl to-chain` / `pgr pl chainnet`.

## Notes

* 2bit inputs are preferred for speed and random access.
* Both sides must use identical sampling parameters (k, syncmer, window);
  `.pgi` inputs carry theirs in the index header.
* The query index is memory-mapped and must be a regular file (`stdin` and
  gzipped indexes are not supported).
* `.pgi` files are not gzip-compressed. The reference index is streamed and
  the query index is memory-mapped (positions are decoded on demand from
  mapped pages), so neither index is materialized in full.

## Examples

1. Align two genomes directly (indexes built automatically):
   ```
   pgr align pgi a.fa.gz b.2bit -o ab.psl
   ```
2. Reuse self-built indexes for repeated comparisons:
   ```
   pgr pgi build a.fa.gz -o a.pgi
   pgr pgi build b.2bit -o b.pgi
   pgr align pgi a.pgi b.pgi --ref-seq a.fa.gz --query-seq b.2bit -o ab.psl
   ```
3. Keep the automatically built indexes:
   ```
   pgr align pgi a.fa b.fa --keep-index -o ab.psl
   ```
4. Tune seed filtering and chaining:
   ```
   pgr align pgi a.fa b.fa -f 20 -c 100 -s 2000 --band 64 -o ab.psl
   ```
5. Detect internal repeats via self-alignment:
   ```
   pgr align pgi genome.fa -o self.psl
   ```
