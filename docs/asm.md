# pgr asm

`pgr asm` provides assembly-related tools: building contigs/unitigs from
reads and mapping reads back to an assembly.

## Subcommands

*   `contig`: Assemble reads into contigs (tadpole-compatible).
*   `unitig`: Assemble reads into maximal unitigs (BCALM-style compaction).
*   `map`: Map reads to a reference requiring perfect matches (bbmap
    perfectmode replacement).

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

*   `-k, --kmer <int>`: K-mer length (default 31; no upper bound — k > 64
    uses multi-word k-mers).
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

*   `-k, --kmer <int>`: K-mer length (default 31; no upper bound — k > 64
    uses multi-word k-mers).
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

*   `-k, --kmer <int>`: Seed k-mer length (default 31, range 1..=64).
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
