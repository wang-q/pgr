# CIGAR — Compact Idiosyncratic Gapped Alignment Report

CIGAR strings encode how a query sequence aligns to a reference (target)
sequence.  They originated in the SAM/BAM format and are also used in PAF
(via the `cg:Z:` tag).  pgr uses CIGAR as the internal representation for
all coordinate projection and identity calculation in its PAF-based
implicit pangenome graph pipeline.

## Operators

| Op | Meaning (query‑centric)          | Target delta | Query delta | SAM  | PAF  |
|----|----------------------------------|:------------:|:-----------:|:----:|:----:|
| M  | Match or mismatch                | len          | len         | ✓    | ✓    |
| =  | Exact match                      | len          | len         | ✓    | ✓    |
| X  | Mismatch                         | len          | len         | ✓    | ✓    |
| I  | Insertion (bases in query only)  | 0            | len         | ✓    | ✓    |
| D  | Deletion (bases in target only)  | len          | 0           | ✓    | ✓    |
| N  | Skipped region (intron)          | len          | 0           | ✓    |      |
| S  | Soft clip                        | 0            | len         | ✓    |      |
| H  | Hard clip                        | 0            | 0           | ✓    |      |
| P  | Padding (silent deletion)        | —            | —           | ✓    |      |

pgr supports `M = X I D`; the remaining operators are SAM‑specific and
do not appear in PAF.

> **pgr v1 convention**: all non‑gap columns are reported as `M` (match
> or mismatch, not distinguished).  `M`, `=`, and `X` behave identically
> for coordinate projection — they all advance both target and query by
> the same amount.  The difference matters only for identity calculation.

## Examples

```
ref:  ACGT---ACG
qry:  ACGTAC-ACG
CIGAR: 4=2I1D3=
```

Walk through column‑by‑column:

| ref | qry | action     | CIGAR segment |
|-----|-----|------------|---------------|
| A   | A   | match      | 1=            |
| C   | C   | match      | 2=            |
| G   | G   | match      | 3=            |
| T   | T   | match      | 4=            |
| –   | A   | insertion  | 1I            |
| –   | C   | insertion  | 2I            |
| A   | –   | deletion   | 1D            |
| C   | C   | match      | 1=            |
| G   | G   | match      | 2=            |
|     |     |            | 3=            |

After run‑length merging: `4=2I1D3=`.

## Coordinate projection

Given a target interval `[t_start, t_end)` and a CIGAR, you can project
to query coordinates by walking the ops and accumulating deltas:

| Op | Target axis      | Query axis       |
|----|------------------|------------------|
| M/= | +len            | +len             |
| X   | +len            | +len             |
| I   | 0               | +len             |
| D   | +len            | 0                |

For reverse‑strand alignments the query axis runs backward (negative
deltas).

In pgr this is implemented by `CigarOp::target_delta()` and
`CigarOp::query_delta()`.  See the projection tests in
`src/libs/paf/cigar.rs` for worked examples.

## Identity

Two identity metrics are computed from CIGAR:

| Metric | Formula | What it measures |
|--------|---------|-----------------|
| **gap‑compressed** (gi) | `matches / (matches + mismatches + #indel_events)` | Homology — each indel counts as 1 event regardless of length |
| **block** (bi) | `matches / (matches + mismatches + indel_bp_total)` | Sequence identity — each indel base counts |

```
CIGAR: 10=5I   →  gi = 10/(10+0+1) = 0.909
                   bi = 10/(10+0+5) = 0.667
```

## CIGAR in pgr

pgr uses CIGAR for two distinct purposes:

1. **`pgr maf to-paf`** — extract CIGAR from MAF alignment strings.
   Two `s` lines are compared base‑by‑base; `-` gaps are translated to
   `I`/`D`, everything else to `M`.

2. **`pgr paf query`** — project target intervals through CIGAR to
   obtain query coordinates, and filter results by identity threshold.

The internal representation is a bit‑packed `u32` per operation
(bits 31‑29: op code, bits 28‑0: length), directly modeled after impg's
`CigarOp` (`impg-0.4.1/src/impg.rs:73-138`).  This gives compact storage
(4 bytes per op) and branch‑free projection arithmetic.

## References

- [CIGAR Strings Explained — Replicon Genetics](https://replicongenetics.com/cigar-strings-explained/)
- [CIGAR Strings — timd.one](https://timd.one/blog/genomics/cigar.php)
- [Redefining the CIGAR String — omicstutorials](https://omicstutorials.com/step-by-step-guide-understanding-and-redefining-the-cigar-string-in-sam-bam-format/)
- [CIGAR Processing — impg DeepWiki](https://deepwiki.com/pangenome/impg/9.3-cigar-processing)
- [Structural variants and the SAM format](https://cmdcolin.github.io/posts/2022-02-06-sv-sam/)
- [PAF specification (lh3/miniasm)](https://github.com/lh3/miniasm/blob/master/PAF.md)
- [SAM v1 spec §1.4 — CIGAR](https://samtools.github.io/hts-specs/SAMv1.pdf)
