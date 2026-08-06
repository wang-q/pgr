# pgr align hybrid

`pgr align hybrid` combines the speed of `pgr align pgi` with the sensitivity
of **LASTZ**, following the FastGA-gapfill idea from the FastGA paper: pgi
produces coarse anchors, LASTZ fills the colinear gaps between consecutive
same-strand anchors, and the two PSL sets are emitted together.

```bash
pgr align hybrid [OPTIONS] <target> <query>
```

Only **syntenic (colinear, same-strand)** searches are supported. Feed the
combined PSL to `pgr pl chainnet --syn` for the final chain/net/axt/maf.

**Note**: `lastz` must be installed and available in your `PATH`.
**Note**: the inputs are converted to 2bit in a tempdir and each box is
extracted with `pgr 2bit range` (random access — no whole-genome sequences are
kept in memory).

## How it works

1.  **Anchors**: `pgr align pgi` runs first, or an existing PSL is reused with
    `--avail-psl` (producible by pgi, FastGA, minimap2...), producing coarse
    PSL blocks.
2.  **Boxes**: for each pair of adjacent, non-overlapping, same-strand anchors
    whose target *and* query gaps both fall in `[--min-gap, --max-gap]`, a
    bounding box is built overlapping the anchors by `--overlap` bp on each
    side (LASTZ seeding buffer).
3.  **Gap fill**: the two sub-sequences of each box are extracted with
    `pgr 2bit range`, written as single-sequence FASTA files, and aligned with
    LASTZ (LAV output).
4.  **Lift**: LASTZ LAV is converted to PSL and the sub-range coordinates are
    lifted back to genomic coordinates.
5.  **Combine**: the anchors and the LASTZ records are written together. No
    dedup happens here — `pgr pl chainnet` handles the overlap/merge downstream.

The default box/overlap parameters follow the paper: 1 kb overlap, and gaps
from 100 bp to 1 Mb. Tune `--overlap`/`--min-gap`/`--max-gap` for your data.

## Sensitivity

On the FastGA-paper-style simulation (§5.1: 6 Mb genomes, target regions of
100–5000 bp at 1–40% divergence, reported as regions recovered/preserved on
both genomes), `pgr align hybrid` recovers **251/600** target regions vs
**256/600** for `lastz` and **186/600** for `pgr align pgi` alone — i.e. hybrid
nearly matches LASTZ's sensitivity while staying close to pgi's speed, exactly
the FastGA-gapfill result. False-positive aligned bases are <1% for all three
engines. Reproduce with `scripts/verify-hybrid-sensitivity.sh`; full table in
`notes/design/pgi-lastz-hybrid.md` §5.1.

## Options

*   `--avail-psl <file>`: Precomputed PSL anchor file (skips the internal
    `align pgi` run). Any aligner's PSL works (pgi, FastGA, minimap2...).
*   `--preset <set01..set07>`: Predefined LASTZ parameter set (see
    `pgr align lastz` for the list). Conservative by default; use a
    close-species preset (e.g. `set01`) for pangenome-scale comparisons.
*   `--overlap <int>`: Box overlap with the anchors in bp (LASTZ seeding
    buffer; default: 1000).
*   `--min-gap <int>`: Shortest gap to fill in bp; smaller gaps are left to
    pgi (default: 100).
*   `--max-gap <int>`: Longest gap to fill in bp; larger gaps are skipped
    (likely novel sequence; default: 1000000).
*   `--query-depth <int>`: Query depth threshold for LASTZ (default: 50).
*   `--lastz-args <string>`: Additional arguments passed directly to LASTZ
    (overrides preset settings).
*   `-o, --outfile <file>`: Output PSL filename.
*   `-p, --parallel <int>`: Number of parallel LASTZ jobs (1..=1024; default:
    8).

## Examples

1.  **Default hybrid alignment**:
    ```bash
    pgr align hybrid ref.fa query.fa -o out.psl
    ```

2.  **Reuse an existing PSL** (from pgi, or any other aligner):
    ```bash
    pgr align hybrid ref.fa query.fa --avail-psl anchors.psl -o out.psl
    ```

3.  **A larger seeding buffer and a close-species preset**:
    ```bash
    pgr align hybrid ref.fa query.fa --preset set01 --overlap 2000 -o out.psl
    ```

4.  **Feed the combined PSL into the chain/net pipeline**:
    ```bash
    pgr align hybrid ref.fa query.fa -o out.psl
    pgr pl chainnet --syn ref.fa query.fa out.psl -o chain_out
    ```