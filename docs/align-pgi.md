# pgr align pgi

`pgr align pgi` aligns two genomes on the pgr k-mer index pipeline:

```
pgr align pgi ref query -o out.psl
  [-f 10] [--min-shared N] [--parallel 8]
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
  (e.g. `ref.fa.gz` → `ref.fa.pgi`, distinct from `ref.fa` → `ref.pgi`).
  A sibling index whose mtime is older than the genome file is rebuilt
  automatically (same convention as the e-kmer repeat-table cache). The
  current `-k`/`--smer`/`--window` (explicit or the defaults) must match the
  cached index's parameters; a mismatch is an error, never a silent reuse
  with different seeds.
- A `.pgi` index is used directly; `--ref-seq`/`--query-seq` may then supply
  the sequences for chain refinement, and are validated against the index
  contig table (count, names, lengths). The sequences must be the ones the
  index was built from: a same-length, same-name but different sequence
  passes the contig check and yields fragmented low-identity alignments
  instead of an error, so keep index and sequences in sync (the automatic
  sibling-index path does this via the mtime check).

With a single input (or `--self` with the same input as query) the genome is
aligned to itself (internal repeats and haplotype-level homology, FastGA's
self mode); exact self-identity hits are dropped, and no same-contig forward
block is emitted on diagonal 0 (FastGA's self-mode wave boundary).

## Options

- `-f`/`--freq`: drop k-mers occurring at least this many times on either
  side (repeats are not seeds);
- `--min-shared`: minimum shared seed length (bp); default is FastGA's plen
  floor (12);
- `--parallel`: rayon thread count (default 8, FastGA `-T` default);
- `-k`/`--smer`/`--window`: sampling parameters for automatic indexing of
  genome inputs only (default 40/8/5, matching FastGA GIX);
- `--keep-index`: write automatically built indexes next to the inputs;
- `--self`: self-alignment; the query input may be omitted or must be the same
  as the reference;
- `--ref-seq`/`--query-seq`: sequences for chain refinement of `.pgi` inputs
  (FASTA or .2bit).

The alignment uses FastGA's `align_contigs` workflow (diagonal buckets,
mid-line wave extension; chain break 2000 bp, min cover 85 bp, 128 bp slide).
Without extension sequences (a `.pgi` pair without `--ref-seq`/`--query-seq`)
each tube is emitted as one geometric block from its seed span; with
sequences the tubes are refined by the wave into scored multi-block records.
The output feeds directly into
`pgr psl to-chain` / `pgr pl chainnet`.

## Terminology

* **Wave**: FastGA's wavefront local aligner (`forward_wave` /
  `reverse_wave`, from Myers' wavefront algorithm). A wavefront is the
  furthest-reaching point per diagonal at a given edit distance; it expands
  one diagonal per edit. The two opposing waves (forward from the mid-line,
  reverse on the mirrored sequences) frame the alignment span, and the exact
  path is reconstructed by a divide-and-conquer edit script. `pgr` uses this
  to refine each tube when sequences are available.
* **Tube**: a seed chain together with its search box (anti-diagonal range ×
  diagonal band); the wave runs inside this box. The name comes from a
  comment in FastGA's `align_contigs` (`FastGA.c:3160`) and its
  `DEBUG_TUBE` macro; it is not an official FastGA API term. There is no
  "cube" in FastGA — the spelling is sometimes confused with "tube".

## Notes

* 2bit inputs are preferred for speed and random access.
* Both sides must use identical sampling parameters (k, syncmer, window);
  `.pgi` inputs carry theirs in the index header.
* The query index is memory-mapped and must be a regular file (`stdin` and
  gzipped indexes are not supported).
* `.pgi` files are not gzip-compressed. The reference index is streamed and
  the query index is memory-mapped (positions are decoded on demand from
  mapped pages), so neither index is materialized in full.
* Automatic indexing applies FastGA `-M` semantics: soft-masked (lowercase)
  bases are replaced by N and produce no seeds or blocks, so a lowercase
  (soft-masked) copy is not aligned against its uppercase twin. Use
  `pgr pgi build --mask` for the same behavior on explicitly built indexes.

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
4. Tune seed filtering:
   ```
   pgr align pgi a.fa b.fa -f 20 -o ab.psl
   ```
5. Detect internal repeats via self-alignment:
   ```
   pgr align pgi genome.fa -o self.psl
   ```
