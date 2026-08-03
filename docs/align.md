# pgr align

`pgr align` performs pairwise genome alignment, emitting PSL blocks. It is a
top-level wrapper around two engines:

- **`pgr align pgi`** — native FastGA-style pipeline: syncmer-sparse k-mer
  indexes, two-index merge, chaining and banded extension. See
  [align-pgi.md](align-pgi.md) and `notes/design/pgi-align.md`.
- **`pgr align lastz`** — wrapper around the external `lastz` aligner with
  Cactus-style presets. See [align-lastz.md](align-lastz.md).

Both engines output PSL blocks, which feed the UCSC-style chain pipeline
(`pgr psl to-chain` → `pgr pl chainnet`) or PAF conversion (`pgr psl to-paf`).

## Core positioning

- **Purpose**: produce pairwise genome alignments in PSL form.
- **Input**: FASTA (plain or `.gz`) or 2bit genomes; `pgi` accepts `.pgi`
  indexes directly.
- **Output**: PSL.
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
