# pgr sd

`pgr sd` detects and analyzes **segmental duplications (SDs)** in one or more
genomes, following the BISER pipeline (search → align → cluster →
decompose → cover):

```
sd search（self-alignment）→ putative SD hits (PSL)
  → sd align（chain/net refinement）→ refined hits (PAF)
    → sd cluster（overlap union-find）→ repeat-family FASTA
      → sd decompose（k-mer chaining）→ elementary SD (BED)
        → sd cover（greedy set cover）→ core duplicon marking
sd cross（cross-genome SD mapping）
sd run（whole pipeline in one command）
```

## sd search

Detects putative SDs by self-aligning a genome and keeping hits meeting the
T2T-CHM13 SD standard (> 1 kbp, > 90% identity).

```
pgr sd search <genome.fa> -o hits.psl
  [--engine pgi|lastz] [--min-len 1000] [--min-identity 0.90]
  [--preset set01] [--query-depth 50] [--parallel 4]
```

* `--engine`: alignment engine — `pgi` (default, native `pgr align pgi`
  self-alignment, no external tools) or `lastz` (external `lastz --self`,
  requires lastz in PATH);
* `--min-len` / `--min-identity`: T2T-CHM13 SD filter (default 1000 bp /
  0.90 identity);
* `--preset` / `--query-depth`: lastz-only parameters (set01..set07 presets,
  query-depth coverage cutoff);
* `--parallel`: worker threads (default 4).

The output PSL is **not** chained; feed it to `pgr sd align`.

## sd align

Refines the putative hits via chain/net (without `--syn`, so rearranged SDs
survive) and merges the MAF blocks into a single PAF.

```
pgr sd align <genome.fa> <hits.psl> -o hits.paf
```

## sd cluster

Clusters overlapping SD mates by union-find: both mates of one hit share a
cluster, and intervals overlapping on the same chromosome are unioned.

```
pgr sd cluster <genome.fa> <hits.paf> -o clusters.dir/
```

Each cluster is written as `cluster_N.fa` with headers in BISER form
`{species}#{chrom}{strand}#{start}#{end}` (0-based).

## sd decompose

Decomposes one cluster FASTA into elementary SDs (k=10 shared-k-mer chaining
with a 50 bp gap tolerance; fragments shorter than 100 bp are dropped).

```
pgr sd decompose <cluster_N.fa> -o cluster_N.elem.bed
```

Output rows: `species<TAB>chrom<TAB>begin<TAB>end<TAB>set_id<TAB>length<TAB>score<TAB>strand`
in genome coordinates.

## sd cover

Marks core duplicons: an elementary SD set covers a hit if any copy overlaps
the hit's query or target interval; a greedy set cover selects the smallest
elementary sets covering all hits, marked `CORE`.

```
pgr sd cover <hits.paf> <elems.bed> -o covered.bed
```

`<elems.bed>` is the merged `pgr sd decompose` output.

## sd cross

Maps SD-like homology from one genome to another (cross-genome counterpart of
search+align).

```
pgr sd cross <target.fa> <query.fa> -o cross.paf
  [--engine pgi|lastz] [--min-len 1000] [--min-identity 0.90]
  [--preset set01] [--query-depth 50] [--parallel 4]
```

`--engine` defaults to `pgi` (native two-genome alignment); `--preset` /
`--query-depth` apply to the lastz engine only.

## sd run

Runs the whole pipeline in one command (search → align → cluster →
decompose → cover) and writes the final CORE-annotated elementary BED to
`<outdir>/out.elem.bed`.

```
pgr sd run <genome.fa> -o sd_out/
  [--engine pgi|lastz] [--min-len 1000] [--min-identity 0.90] [--preset set01]
```

## Notes

* The SD filter follows the T2T-CHM13 standard: block length > 1 kbp and
  identity > 90%, computed as `(matches + repeats) / block_length`.
* The default `pgi` engine is fully native (no external tools); the `lastz`
  engine requires lastz in PATH and is kept for comparison.
* Chain/net refinement always runs without `--syn` (`pgr pl chainnet`), so
  SDs associated with rearrangements are kept.
* Both `sd search` and `sd cross` share the same engines and parameters.

## Examples

1. Run the full SD pipeline on a genome:
   ```
   pgr sd run genome.fa -o sd_out/
   ```
2. Step-by-step with intermediate files:
   ```
   pgr sd search genome.fa -o hits.psl
   pgr sd align genome.fa hits.psl -o hits.paf
   pgr sd cluster genome.fa hits.paf -o clusters/
   pgr sd decompose clusters/cluster_1.fa -o cluster_1.elem.bed
   pgr sd cover hits.paf elems.bed -o covered.bed
   ```
3. Cross-genome SD mapping:
   ```
   pgr sd cross genomeA.fa genomeB.fa -o cross.paf
   ```
4. Use the external lastz engine:
   ```
   pgr sd search genome.fa --engine lastz -o hits.psl
   ```
