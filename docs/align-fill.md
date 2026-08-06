# pgr align fill

`pgr align fill` combines the speed of `pgr align pgi` with the sensitivity of
**LASTZ** following the FastGA-gapfill idea: pgi produces coarse anchors, and
LASTZ fills the colinear gaps between consecutive same-strand anchors. The two
PSL sets are emitted together.

```bash
pgr align fill [OPTIONS] <target> <query>
```

Only **syntenic (colinear, same-strand)** searches are supported. Feed the
combined PSL to `pgr pl chainnet --syn` for the final chain/net/axt/maf. For
filling everything the anchors do not cover (contig ends, anchor-free
contigs), use `pgr align rest`.

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
    lifted back to genomic coordinates (the `pgr psl lift` logic).
5.  **Combine**: the anchors and the LASTZ records are written together. No
    dedup happens here — `pgr pl chainnet` handles the overlap/merge
    downstream.

`--overlap` defaults to 1 kb (the paper's box overlap), `--min-gap` to 100 bp
(smaller gaps are left to pgi) and `--max-gap` is unlimited by default.

## Options

*   `--avail-psl <file>`: Precomputed PSL anchor file (skips the internal
    `align pgi` run). Any aligner's PSL works (pgi, FastGA, minimap2...).
*   `--preset <set01..set07>`: Predefined LASTZ parameter set (see
    `pgr align lastz` for the list).
*   `--overlap <int>`: Box expansion beyond the gap in bp (LASTZ seeding
    buffer; default: 1000).
*   `--min-gap <int>`: Shortest gap to fill in bp; smaller gaps are left to
    pgi (default: 100).
*   `--max-gap <int>`: Longest gap to fill in bp; larger gaps are skipped
    (likely novel sequence; default: no limit).
*   `--query-depth <int>`: Query depth threshold for LASTZ (default: 50).
*   `--lastz-args <string>`: Additional arguments passed directly to LASTZ
    (overrides preset settings).
*   `-o, --outfile <file>`: Output PSL filename.
*   `-p, --parallel <int>`: Number of parallel LASTZ jobs (1..=1024; default:
    8).

## Examples

1.  **Default gap fill**:
    ```bash
    pgr align fill ref.fa query.fa -o out.psl
    pgr pl chainnet --syn ref.fa query.fa out.psl -o chain_out
    ```

2.  **Reuse an existing PSL** (from pgi, or any other aligner):
    ```bash
    pgr align fill ref.fa query.fa --avail-psl anchors.psl -o out.psl
    ```

3.  **A larger seeding buffer and a close-species preset**:
    ```bash
    pgr align fill ref.fa query.fa --preset set01 --overlap 2000 -o out.psl
    ```
