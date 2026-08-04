# pgr align

`pgr align` performs pairwise genome alignment. It is a top-level wrapper
around two engines:

- **`pgr align pgi`** — native FastGA-style pipeline: syncmer-sparse k-mer
  indexes, two-index merge, chaining and banded extension. See
  [align-pgi.md](align-pgi.md) and `notes/design/pgi-align.md`.
- **`pgr align lastz`** — wrapper around the external `lastz` aligner with
  Cactus-style presets. See [align-lastz.md](align-lastz.md).

`pgr align pgi` emits PSL blocks directly; `pgr align lastz` emits LAV files
(one per target/query pair), which convert to PSL with `pgr lav to-psl`.
Either way the PSL blocks feed the UCSC-style chain pipeline
(`pgr psl to-chain` → `pgr pl chainnet`) or PAF conversion (`pgr psl to-paf`).

## Core positioning

- **Purpose**: produce pairwise genome alignments (PSL, or LAV via lastz).
- **Input**: FASTA (plain or `.gz`) or 2bit genomes; `pgi` accepts `.pgi`
  indexes directly.
- **Output**: PSL (`align pgi`) or LAV (`align lastz`).
- **Complements**: `pgr pgi build/stat` (index lifecycle), `pgr dist`
  (index/sequence distances), `pgr sd search --engine pgi|lastz`
  (self-alignment reuse).
- **External dependency**: only `align lastz` requires the `lastz` binary.

## Examples

1. FastGA-style alignment with default parameters:

   ```bash
   pgr align pgi ref.fa query.fa -o out.psl
   ```

2. LASTZ alignment:

   ```bash
   pgr align lastz ref.fa query.fa -o lastz_out
   ```

   `align lastz` writes one LAV file per target/query pair into the output
   directory; convert them with `pgr lav to-psl` (see
   [align-lastz.md](align-lastz.md)).
