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

## Background

`pgr align pgi` is a native reimplementation of the FastGA alignment pipeline
(Myers, Durbin and Zhou, *FastGA: fast genome alignment*, Bioinformatics
Advances 5(1):vbaf238, 2025, DOI 10.1093/bioadv/vbaf238): syncmer-sparse
k-mer indexes, a linear merge of the two sorted index streams to find
adaptive-seed hits, diagonal-bucket seed chaining, and a wavefront local
aligner. FastGA is typically an order of magnitude faster than tools of
comparable sensitivity (e.g. two 2 Gbp bat genomes in 2.1 min with 8 threads
and 5.7 GB RAM), while its sensitivity is close to minimap2's and slightly
below LastZ's, which remains the most sensitive aligner in the paper's
benchmarks.

Like FastGA, pgi deliberately separates *genome alignment* (finding all
statistically significant local alignments, with maximal internal gaps of
about 40 bp) from *homology inference* (chaining those alignments across
larger gaps). The output is a PSL of local alignments; chaining and netting
are left to `pgr psl to-chain` / `pgr pl chainnet`, exactly as FastGA leaves
chaining to a second step. The paper also evaluates a FastGA + LastZ hybrid
("FastGA-gapfill": FastGA anchors, LastZ fills the gaps between them with a
1 kb overlap for seeding), whose sensitivity approaches LastZ's at
19.3×–137.5× its speed; the same hybrid design with `pgr align pgi` anchors
and `pgr align lastz` gap filling is implemented as `pgr align hybrid` (see
[align-hybrid.md](align-hybrid.md)).

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

* Follows the UCSC/LastZ convention: the reference genome comes first,
  the query second (`pgr align pgi <ref> <query>`); in the output PSL the
  reference is the target (`tName`) and the second input the query
  (`qName`). This is the opposite positional order of FastGA's
  `FastGA <query> <target>` but matches `pgr pl chainnet`.
* Seed semantics are always anchored to the reference genome: seeds are
  defined and emitted from the reference (the first input) and matched
  against the query. Swapping the inputs is therefore not symmetric, and
  results differ from FastGA's query-anchored model by design; the same
  reference yields reproducible seeds regardless of which query set is
  compared against it.
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
