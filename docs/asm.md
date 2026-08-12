# pgr asm

`pgr asm` provides assembly-related tools: building contigs/unitigs from
reads and mapping reads back to an assembly.

## Subcommands

*   `contig`: Assemble reads into contigs (tadpole-compatible).
*   `unitig`: Assemble reads into maximal unitigs (BCALM-style compaction).
*   `map`: Map reads to a reference requiring perfect matches (bbmap
    perfectmode replacement).
*   `ovlp`: Find exact overlaps between unitigs (OLC stage 1).
*   `layout`: Chain unitigs into layouts from an overlap PAF (OLC stage 2).
*   `cns`: Stitch layouts into consensus contigs (OLC stage 3).
*   `olc`: Assemble reads into contigs via multi-k unitig OLC (full
    pipeline).

---

## contig

Assembles reads into contigs through the k-mer graph, reproducing the
BBTools `tadpole.sh` contig mode (the default mode when no `ecc`/`extend`
flag is set). K-mers are counted with a quality gate, contigs are seeded
from k-mers above a depth threshold and extended greedily in both
directions, stopping at branches and dead ends, then bubbles are resolved
and contigs are sorted longest-first. This replaces the tadpole assembly
steps of the anchr `2_insert_size` and `unitigs` flows. Bubble popping is
enabled by default (tadpole `popbubbles=t`); `--no-bubbles` keeps
parallel-path contigs separate (tadpole `popbubbles=f`), which preserves
both branches of a bubble instead of merging them into one representative
path. For the strict graph-compression counterpart (maximal unitigs, no
seeded extension), see [`pgr asm unitig`](#unitig).

```bash
pgr asm contig [OPTIONS] <infiles>...
```

### Options

*   `-k, --kmer <int>`: K-mer length (default 31; up to 128, the k-mer key
    table limit — k > 64 uses multi-word k-mers).
*   `-o, --outfile <file>`: Output FASTA filename (default: stdout).
*   `--min-contig-len <int>`: Minimum contig length (default:
    `max(124, 2*k)`).
*   `--min-count-seed <int>`: Minimum k-mer depth to seed a contig
    (tadpole `mincountseed`, default 3).
*   `--min-coverage <float>`: Minimum mean k-mer coverage for a contig
    (tadpole `mincoverage`, default 1.0).
*   `--no-bubbles`: Keep parallel-path contigs separate; disable bubble
    popping (default: bubble popping on, matching tadpole `popbubbles=t`).
*   `-p, --parallel <int|auto>`: Accepted for tadpole.sh compatibility;
    ignored (processing is deterministic single-pass).

### Examples

1.  **Assemble contigs from corrected reads (anchr unitigs step)**:
    ```bash
    pgr asm contig pe.cor.fa -o unitigs_K31.fasta --kmer 31
    ```

2.  **Assemble from paired-end reads (anchr 2_insert_size step)**:
    ```bash
    pgr asm contig R1.fq.gz R2.fq.gz -o contigs.fasta
    ```

3.  **Raise the minimum contig length**:
    ```bash
    pgr asm contig in.fq -o out.fasta --min-contig-len 500
    ```

---

## unitig

Assembles reads into maximal unitigs through the k-mer graph, following the
BCALM 2 compaction semantics (GATB `ograph.cpp` `graph3`): every solid
k-mer (count >= 3) extends in both directions only while it has exactly one
solid successor whose own predecessor is also unique, so the assembly stops
at branches, junctions, coverage gaps, and loops. Parallel paths stay
separate (no bubble popping), and the result is independent of the k-mer
scan order.

This is the strict graph-compression counterpart of
[`pgr asm contig`](#contig), whose seeded contig mode keeps extending
through weak branches (tadpole-compatible). Unitigs are best suited to
high-coverage or error-corrected input, such as the anchr `unitigs` step's
`pe.cor.fa`.

```bash
pgr asm unitig [OPTIONS] <infiles>...
```

### Options

*   `-k, --kmer <int>`: K-mer length (default 31; up to 128, the k-mer key
    table limit — k > 64 uses multi-word k-mers).
*   `-o, --outfile <file>`: Output FASTA filename (default: stdout).
*   `--min-contig-len <int>`: Minimum unitig length (default:
    `max(124, 2*k)`).
*   `--min-count-seed <int>`: Solid k-mer count threshold (default 3, like
    bcalm `-abundance-min`).
*   `--links`: Append BCALM-style `L:+:<to>:<sign>` links to unitig FASTA
    headers (links connect unitigs sharing an endpoint (k-1)-mer).
*   `--gfa`: Emit a GFA 1.0 graph (`H`/`S`/`L` lines, overlap `(k-1)M`)
    instead of FASTA.
*   `-p, --parallel <int|auto>`: Accepted for compatibility; ignored
    (processing is deterministic).

### Examples

1.  **Assemble unitigs from corrected reads (anchr unitigs step)**:
    ```bash
    pgr asm unitig pe.cor.fa -o unitigs_K31.fasta --kmer 31
    ```

2.  **Assemble from paired-end reads**:
    ```bash
    pgr asm unitig R1.fq.gz R2.fq.gz -o unitigs.fasta
    ```

3.  **Raise the minimum unitig length**:
    ```bash
    pgr asm unitig in.fq -o out.fasta --min-contig-len 500
    ```

4.  **Raise the solid k-mer threshold (bcalm `-abundance-min` equivalent)**:
    ```bash
    pgr asm unitig in.fq -o out.fasta --min-count-seed 5
    ```

5.  **Emit the unitig graph as GFA**:
    ```bash
    pgr asm unitig in.fq -o unitigs.gfa --gfa
    ```

---

## map

Maps reads to a reference (typically an assembly) requiring every read to
match exactly: no mismatches and no gaps, mirroring BBTools
`bbwrap.sh perfectmode maxindel=0 strictmaxindel`. This replaces the bbwrap
call of the anchr `anchors` flow, whose downstream only needs the
mapped/unmapped counts and the per-base coverage.

Mapping is seed-and-verify: the reference's canonical k-mers are indexed
once (sorted, radix), each read seeds on its first k-mer, and every
candidate position is verified over the full read length (forward or
reverse strand). Reads matching multiple positions are reported at all of
them (`ambiguous=all` semantics).

Per-base coverage is not accumulated in memory here: it is derived from the
mapped SAM with `pgr sam to-rg` and `pgr rg coverage` (see examples), which
is cheaper when the reference is large and keeps the mapping command
single-purpose.

With `--paired`, the two read files are interleaved as R1/R2 pairs: a pair
is mapped only when both ends match perfectly, and the SAM carries pair
flags (0x1/0x2/0x40/0x80), mate coordinates and signed TLEN. This supports
insert-size estimation with `pgr sam ihist` (anchr `2_insert_size` step).

```bash
pgr asm map [OPTIONS] <ref.fa> <reads.fq...>
```

### Options

*   `-k, --kmer <int>`: Seed k-mer length (default 31, range 1..=128).
*   `--outm <file>`: SAM output of perfectly matched reads.
*   `--outu <file>`: SAM output of unmapped reads.
*   `--paired`: Map reads as R1/R2 pairs (exactly 2 read files; pairs with
    an unmapped end go to `--outu`).
*   `--max-reads <int>`: Stop after processing this many read records
    (pairs count as two).
*   `-p, --parallel <int|auto>`: Worker threads (real parallelism via
    rayon; output stays deterministic and in input order).

### Examples

1.  **Map reads back to an assembly (anchr anchors step)**:
    ```bash
    pgr asm map UT.fasta R1.fq.gz R2.fq.gz \
        --outm mapped.sam --outu unmapped.sam
    ```

2.  **Derive per-base coverage from the mapped SAM (anchr anchors step)**:
    ```bash
    pgr sam to-rg mapped.sam | pgr rg coverage stdin -m 2 -o cov.json
    ```

3.  **Map paired reads and estimate the insert size (anchr 2_insert_size
    step)**:
    ```bash
    pgr asm map UT.fasta R1.fq.gz R2.fq.gz --paired \
        --outm mapped.sam --outu unmapped.sam --max-reads 1000000
    pgr sam ihist mapped.sam -o insert_size.ihist.txt
    ```

4.  **Use a longer seed k-mer**:
    ```bash
    pgr asm map ref.fa reads.fq.gz -k 41 --outm mapped.sam
    ```

---

## ovlp

Finds exact overlaps between unitigs by seeding a canonical k-mer index
with the boundary k-mers of every unitig and verifying each candidate by
extension, so overlaps are exact and error-free (unitigs come from the de
Bruijn graph). This is the overlap stage of the OLC assembly pipeline; the
caller is expected to assemble unitigs at several k values first and pass
the FASTA files here.

Overlaps are written as PAF with an `ov:A:D` (dovetail) or `ov:A:C`
(contain) tag. Unitig names are prefixed with the input file stem
(`stem:name`) so identical `unitig_<id>` names across k files stay unique.

```bash
pgr asm ovlp [OPTIONS] <infiles>...
```

### Options

*   `-o, --outfile <file>`: Output PAF filename (default: stdout).
*   `--overlap-k <int>`: Seed k-mer length (default 17; clamped to the
    shortest unitig).
*   `--min-overlap <int>`: Minimum accepted overlap length in bases
    (default 34).

### Examples

1.  **Overlap unitigs from two k values**:
    ```bash
    pgr asm ovlp k21.fa k51.fa -o ovlp.paf
    ```

2.  **Raise the seed and minimum overlap**:
    ```bash
    pgr asm ovlp unitigs.fa -o ovlp.paf --overlap-k 21 --min-overlap 51
    ```

---

## layout

Builds greedy layouts from the exact overlaps produced by `pgr asm ovlp`:
every unitig end gets its best extension edge, unplaced unitigs are seeded
longest-first, and chains grow in both directions through mutual-best
junctions. Ambiguous junctions (two near-equal best partners, e.g. repeats)
and non-reciprocal edges stop the chain, so branches stay separate and no
heuristic picks a bubble path.

The unitig FASTA files must be the same files passed to `pgr asm ovlp` (the
`stem:name` prefixes are re-derived here and must match the PAF names). The
PAF file is the first positional argument.

Output is a layout TSV (no header), one line per step:
`contig_id<TAB>step<TAB>unitig_name<TAB>strand<TAB>q_start<TAB>q_end<TAB>overlap_len`
where `q_start`/`q_end` is the unitig's interval in the contig and
`overlap_len` is the exact overlap with the previous step (0 for the first
step).

```bash
pgr asm layout <paf> <infiles>... -o layout.tsv
```

### Examples

1.  **Layout overlaps from two k values**:
    ```bash
    pgr asm layout ovlp.paf k21.fa k51.fa -o layout.tsv
    ```

---

## cns

Stitches the layouts produced by `pgr asm layout` into consensus contigs.
Overlaps are exact, so each layout is walked in order, every unitig is
oriented by its strand, and only the bases beyond the exact overlap with
the previous step are appended. A layout whose overlapping bases disagree
with the already-stitched contig is reported as an error (exact overlaps
must agree).

The unitig FASTA files must be the same files passed to `pgr asm ovlp` and
`pgr asm layout`. The layout TSV is the first positional argument.

Output is FASTA (`>contig_<id>,len=...,cov=...`, 70-column wrap, longest
first); `cov` is the approximate unitig depth (sum of unitig lengths over
the contig length).

```bash
pgr asm cns <layout.tsv> <infiles>... -o contigs.fa
```

### Options

*   `-o, --outfile <file>`: Output FASTA filename (default: stdout).
*   `--min-contig-len <int>`: Minimum contig length (default 500).

### Examples

1.  **Consensus from a layout**:
    ```bash
    pgr asm cns layout.tsv k21.fa k51.fa -o contigs.fa
    ```

2.  **Drop short contigs**:
    ```bash
    pgr asm cns layout.tsv unitigs.fa -o contigs.fa --min-contig-len 500
    ```

---

## olc

Runs the full OLC pipeline in memory: for every k in `--kmer` the reads are
assembled into maximal unitigs (`pgr asm unitig` semantics), all unitigs are
pooled as pseudo-reads, exact overlaps are found (`pgr asm ovlp`), layouts
are built greedily (`pgr asm layout`), and each layout is stitched into a
consensus contig (`pgr asm cns`). See `notes/design/olc.md`.

Unitigs are named `k<k>:unitig_<id>` so the per-k sets stay distinguishable
and reproducible. Overlaps are exact (error-free unitigs), layouts stop at
ambiguous junctions and non-reciprocal edges, and no bubble heuristics are
applied.

```bash
pgr asm olc [OPTIONS] <infiles>...
```

### Options

*   `-k, --kmer <int,int,...>`: Comma-separated k-mer lengths for the unitig
    sets (default `21,51,81`).
*   `-o, --outfile <file>`: Output FASTA filename (default: stdout).
*   `--min-count-seed <int>`: Solid k-mer count threshold for unitig
    assembly (default 3).
*   `--overlap-k <int>`: Seed k-mer length for overlap detection (default
    17).
*   `--min-overlap <int>`: Minimum accepted overlap length in bases (default
    34).
*   `--min-contig-len <int>`: Minimum output contig length (default 500).
*   `--keep-dir <dir>`: Write the intermediate unitigs/ovlp/layout files
    for debugging or re-running the stage commands separately.

### Examples

1.  **Assemble a small metagenome with three k values**:
    ```bash
    pgr asm olc reads.fq.gz -o contigs.fa --kmer 21,51,81
    ```

2.  **Keep the intermediates and raise the minimum contig length**:
    ```bash
    pgr asm olc R1.fq.gz R2.fq.gz -o contigs.fa \
        --kmer 21,51,81 --min-contig-len 1000 --keep-dir stage/
    ```
