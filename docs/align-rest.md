# pgr align rest

The name `rest` means "the rest of the genome". Given a PSL of anchors (from
`pgr align pgi`, or any aligner's PSL via `--avail-psl`), `pgr align rest`
looks **beyond those alignments**: it computes the target-side regions the
anchors do not cover (trim -> excise small anchors -> whole-genome holes) and
tries again to find homology for them with **LASTZ**. The query side gets the
same treatment, and a k-mer prefilter pairs target/query holes so only likely
pairs are aligned. `pgr align fill` only fills the gaps between anchors;
`rest` fills everything else.

```bash
pgr align rest [OPTIONS] <target> <query>
```

Only **syntenic (colinear, same-strand)** searches are supported. Feed the
combined PSL to `pgr pl chainnet --syn` for the final chain/net/axt/maf.

**Note**: `lastz` must be installed and available in your `PATH`.
**Note**: the inputs are converted to 2bit in a tempdir and each hole is
extracted with `pgr 2bit range` (random access — no whole-genome sequences are
kept in memory).

## How it works

1.  **Anchors**: `pgr align pgi` runs first, or an existing PSL is reused with
    `--avail-psl` (producible by pgi, FastGA, minimap2...).
2.  **Per-side holes** (independent 1D runlist operations, no 2D
    coordinates): the anchor spans on each side are trimmed (`--trim`),
    anchors shorter than `--min-anchor` are excised, and the whole-genome
    complement is computed for every contig.
3.  **Query side**: the query holes are extracted with `pgr 2bit range` and
    concatenated into one multi-sequence FASTA (LASTZ supports multiple query
    sequences).
4.  **Prefilter**: every hole is sampled into a k-mer hash set (closed
    syncmers by default) and target/query holes sharing at least
    `--min-shared` sampled k-mers are paired — only those pairs are aligned,
    avoiding full-query scans per hole.
5.  **Fill**: each paired target hole is extracted as a single-sequence FASTA
    and aligned with LASTZ against its paired query hole (LAV output).
    Unpaired target holes are skipped by default (`--unmatched full` aligns
    them against the merged query-holes set instead; `--sampler none` skips
    the prefilter entirely).
6.  **Lift**: LASTZ LAV is converted to PSL and the sub-range coordinates are
    lifted back to genomic coordinates (the `pgr psl lift` logic).
7.  **Combine**: the anchors and the LASTZ records are written together. No
    dedup happens here — `pgr pl chainnet` handles the overlap/merge
    downstream.

`--trim` defaults to 500 bp (also covers the 1-11 bp pgi boundary error),
`--min-anchor` to 500 bp and `--max-gap` is unlimited by default. The
prefilter defaults to syncmer s-mer 17 / window 5 / min-shared 1: fast with a
small coverage cost; `--smer 15` recovers most of it, `--sampler none`
disables the prefilter for the full-coverage path.

## Options

*   `--avail-psl <file>`: Precomputed PSL anchor file (skips the internal
    `align pgi` run). Any aligner's PSL works (pgi, FastGA, minimap2...).
*   `--preset <set01..set07>`: Predefined LASTZ parameter set (see
    `pgr align lastz` for the list).
*   `--trim <int>`: Shrink each anchor span by this many bp on both ends
    before the complement (default: 500).
*   `--min-anchor <int>`: Excise anchors shorter than this in bp before the
    complement (default: 500).
*   `--max-gap <int>`: Skip holes longer than this in bp (likely novel
    sequence; default: no limit).
*   `--sampler <syncmer|minimizer|none>`: Fragment prefilter sampler
    (default: syncmer; none = every target hole against the full query-holes
    set).
*   `--smer <int>`: Syncmer s-mer length (default: 17; 15 for higher
    coverage).
*   `--kmer <int>`: Minimizer k-mer length (default: 17; used with
    `--sampler minimizer`).
*   `--window <int>`: Prefilter sampler window size (default: 5).
*   `--min-shared <int>`: Minimum shared sampled k-mers to pair a
    target/query hole (default: 1).
*   `--top-k <int>`: Keep only the top-K query holes per target hole by
    shared k-mers (default: no limit).
*   `--unmatched <skip|full>`: Unpaired target holes: skip (default) or
    align against the full query-holes set.
*   `--query-depth <int>`: Query depth threshold for LASTZ (default: 50).
*   `--lastz-args <string>`: Additional arguments passed directly to LASTZ
    (overrides preset settings).
*   `-o, --outfile <file>`: Output PSL filename.
*   `-p, --parallel <int>`: Number of parallel LASTZ jobs (1..=1024; default:
    8).

## Examples

1.  **Default complement fill**:
    ```bash
    pgr align rest ref.fa query.fa -o out.psl
    pgr pl chainnet --syn ref.fa query.fa out.psl -o chain_out
    ```

2.  **Reuse an existing PSL** (from pgi, or any other aligner):
    ```bash
    pgr align rest ref.fa query.fa --avail-psl anchors.psl -o out.psl
    ```

3.  **Drop tiny anchors and skip novel-sequence holes**:
    ```bash
    pgr align rest ref.fa query.fa --min-anchor 1000 --max-gap 100000 -o out.psl
    ```
