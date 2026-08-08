# pgr rept

Repeat detection and masking. `pgr rept` provides commands that detect
repetitive regions of a genome and emit runlist JSON ready for
`pgr fa mask`.

Command names use two dimensions:

*   prefix `e` = external repeat library required, `s` = self (genome only);
*   suffix `kmer` / `align` = detection mechanism (k-mer counting or
    alignment).

## Repeat libraries

Repeat-masking workflows compare a genome against a library of known repeat
consensus sequences. Download every library into a dedicated directory and
run the commands below from there:

```bash
mkdir -p "$HOME/data/repeats"
cd "$HOME/data/repeats"
```

Files under `tests/pgr/` (e.g. `tncentral.fa.gz`) are small test fixtures for
the test suite, not a working repeat database.

Libraries are stored as gzipped FASTA; `pgr` reads any `.gz` directly, so no
conversion is needed. BGZF is only relevant for random-access workflows
(`pgr fa range`, paf graph TSVs), which repeat libraries are not part of.

Before use, every library is normalized with `pgr fa filter` (uppercase,
IUPAC codes -> N, dashes stripped, duplicate IDs removed). Cleaning is
required for `pgr rept e-align`: RepBase consensus sequences contain
ambiguous (IUPAC) codes — a native feature of repeat consensi — that
otherwise make the alignment pass output nothing.

### TnCentral

TnCentral is a database of prokaryotic insertion sequences. Download the
complete IS library and prepare it for `pgr rept e-kmer`:

The site may reject command-line downloaders (`curl` returns HTTP 403 even
with a browser User-Agent). If that happens, download from a browser instead:
open <https://tncentral.ncc.unesp.br/data_download>, click "download Fasta
format" under "TnCentral Dataset", and save the archive into
`$HOME/data/repeats/`. Unpack it, rename the FASTA inside to `tncentral.fa`,
then clean and gzip it as below.

```bash
# Download and unpack
curl -LO https://tncentral.ncc.unesp.br/api/download_blast/nc/tn_in_is
unzip -j tn_in_is 'tncentral_integrall_isfinder.fa'

# A few records in this file lost their header newline (embedded '>' in
# sequence lines); split them back out first. The rule relies on NCBI
# accessions ending in a digit, which holds for this file.
# Then clean: uppercase, IUPAC -> N, strip dashes, drop duplicate IDs
perl -ne 'if (/^>/) { print; next } if (/>/) { s/>([A-Za-z0-9_.]*\d)/\n>$1\n/g } print' \
    'tncentral_integrall_isfinder.fa' |
    pgr fa filter stdin --upper --iupac --dash --uniq |
    gzip -9 -c > tncentral.fa.gz

# Sanity check and quality filter
pgr fa size tncentral.fa.gz
pgr dist mini tncentral.fa.gz -k 17 -w 5 -p 8 |
    tva filter stdin --ge 5:0.9
```

### RepBase

RepBase is the classic repeat database used by RepeatMasker. The distribution
tarball contains an EMBL-format library; convert it to FASTA with readseq:

```bash
# Download and unpack
curl -LO https://github.com/wang-q/ubuntu/releases/download/20190906/repeatmaskerlibraries-20140131.tar.gz
tar xvfz repeatmaskerlibraries-20140131.tar.gz Libraries/RepeatMaskerLib.embl

# https://sourceforge.net/projects/readseq/
java -jar ~/bin/readseq.jar -f fa Libraries/RepeatMaskerLib.embl

# Clean: uppercase, IUPAC -> N, strip dashes, drop duplicate IDs
pgr fa filter Libraries/RepeatMaskerLib.embl.fasta --upper --iupac --dash --uniq |
    gzip -9 -c > repbase.fa.gz

# Sanity check and quality filter
pgr fa size repbase.fa.gz
pgr dist mini repbase.fa.gz -k 17 -w 5 -p 8 |
    tva filter stdin --ge 5:0.9
```

### Dfam

Dfam is a curated database of transposable element families, the repeat
database RepeatMasker ships with. Each release provides consensus sequences
in `families/`. Two products fit pgr's consensus-based masking workflows
(`pgr rept e-kmer`, or a `pgr align pgi`-based masker):

*   `Dfam-RepeatMasker.lib.gz` — curated family consensus in FASTA, the same
    library RepeatMasker consumes via `-lib`. Download and use directly.
*   `Dfam-curated_only-1.embl.gz` — the same curated families as EMBL records
    with per-family metadata; convert to FASTA with readseq (same as RepBase).

Download and prepare (Dfam 4.0; `current/` in the URL also points to the
latest release):

```bash
# FASTA variant: ready to use, no conversion
curl -LO https://dfam.org/releases/Dfam_4.0/families/Dfam-RepeatMasker.lib.gz

# Clean: uppercase, IUPAC -> N, strip dashes, drop duplicate IDs
gzip -dc Dfam-RepeatMasker.lib.gz |
    pgr fa filter stdin --upper --iupac --dash --uniq |
    gzip -9 -c > dfam.fa.gz

# Sanity check and quality filter
pgr fa size dfam.fa.gz
pgr dist mini dfam.fa.gz -k 17 -w 5 -p 8 |
    tva filter stdin --ge 5:0.9
```

## Example run: E. coli MG1655

With the three libraries prepared above, an end-to-end pass on E. coli
MG1655 (`tests/genome/mg1655.fa.gz`, NC_000913, 4,641,652 bp) looks like
this:

```bash
# Run each library serially: FastK-based commands are not safe to run in
# parallel (three concurrent e-kmer runs crashed FastK with SIGSEGV)
for lib in tncentral repbase dfam; do
    pgr rept e-kmer "$HOME/data/repeats/$lib.fa.gz" tests/genome/mg1655.fa.gz \
        --keep-index -o "$lib.json"
done

# Tandem repeats complement the interspersed-repeats pass
pgr rept trf tests/genome/mg1655.fa.gz -o trf.json

# Coverage per library, and pairwise overlap
pgr fa size tests/genome/mg1655.fa.gz > mg1655.sizes
pgr runlist stat mg1655.sizes tncentral.json
pgr runlist statop mg1655.sizes repbase.json dfam.json
```

Measured on MG1655 (RepeatMasker reference run with the Dfam library:
49,379 bp, 1.06%; the same run with the TnCentral library masks 163,249 bp
raw, 89,743 bp after the ≥50 bp / ≥70% identity filter used by `e-align` —
see notes/ecoli-repeats.md §2.7):

| Library | Intervals | Covered (bp / %) | RM overlap |
| :--- | :--- | :--- | :--- |
| TnCentral | 48 | 56,969 / 1.23% | 90.7% of RM |
| RepBase | 38 | 42,763 / 0.92% | 86.0% of RM |
| Dfam | 39 | 42,386 / 0.91% | 85.5% of RM |
| `trf` | 84 | 18,768 / 0.40% | 0.8% of RM |

Notes:

*   The three libraries agree on a common core (~42.5 kb): 99.8% of Dfam
    hits fall inside RepBase and 99.6% of RepBase hits inside TnCentral.
    TnCentral adds ~14 kb of extra coverage (E. coli IS elements), which is
    fine for a "mask more, miss less" masking pass.
*   TnCentral covers 90.7% of the RepeatMasker intervals but masks ~15% more
    than it; RepBase/Dfam cover 85–86% and miss ~7 kb. For a TE-poor
    prokaryote all numbers are small; the real differences show on
    TE-rich eukaryotic genomes.
*   `trf` does not overlap `e-kmer` at all (tandem vs interspersed
    complement). `e-kmer` (TnCentral) + `trf` covers ~75.7 kb (1.63%) and
    ~91.5% of the RepeatMasker intervals.
*   Run libraries serially: concurrent `e-kmer` runs crashed FastK with
    SIGSEGV on RepBase. `--keep-index` caches each library table for reuse
    (needs a writable directory next to the library).

## RepeatMasker (reference)

RepeatMasker remains the reference annotation tool. Example run through a
native installation (4.2.4; TRF + RMBlast 2.14.1 — the CBP build, which also
runs on old glibc — configured via `perl ./configure`), using TnCentral as a
custom library, then converting its `.out` to a GFF runlist for comparison
(`<rm_dir>` is the RepeatMasker installation directory):

```bash
<rm_dir>/RepeatMasker genome.fa -lib tncentral.fa -pa 8 -e rmblast -dir rm_out
perl <rm_dir>/util/rmOutToGFF3.pl rm_out/genome.fa.out > genome.rm.gff
pgr gff runlist genome.rm.gff -o genome.rm.json
```

Notes:

*   `-lib` takes any FASTA library (TnCentral, RepBase, Dfam, ...) directly;
    it does **not** accept gzipped files — decompress first. `-species` needs a
    FamDB/Dfam installation, which is not configured here.
*   RepeatMasker's `configure` validates that RMBLAST_DIR contains six
    executables (rmblastn, dustmasker, makeblastdb, blastdbcmd,
    blastdb_aliastool, blastn). The install used for the reference runs
    satisfies this with a merged directory where rmblastn/makeblastdb are the
    CBP build and the other four are symlinks to the official package purely
    to pass the check; the `-lib` flow only ever calls rmblastn and
    makeblastdb.
*   On the 10-strain E. coli cohort (notes/ecoli-repeats.md §2.7), all `e-kmer`
    / `e-align` hits fall inside the RepeatMasker intervals (99%+), while
    RepeatMasker's raw output is more permissive: TnCentral hits total
    ~2.15 Mb (4.1% of the cohort) before filtering. After applying the same
    thresholds as `e-align` (span ≥ 50 bp, identity ≥ 70%), the totals
    converge (RepeatMasker 1.37 Mb vs LASTZ 1.36 Mb) and LASTZ covers 93.6%
    of the RepeatMasker intervals.

## Subcommands

| Subcommand | Description |
| :--- | :--- |
| `e-kmer` | Identify repeats against an external library (k-mer) |
| `e-align` | Identify repeats against an external library (alignment) |
| `masker` | Simulate RepeatMasker (TRF + `rmblastn` library search + TRF) |
| `s-kmer` | Identify repetitive regions by self k-mer depth (no library) |
| `s-align` | Identify repetitive regions by self alignment |
| `trf` | Identify tandem repeats via `trf` |

All six emit runlist JSON ready for `pgr fa mask`:

```bash
pgr rept e-kmer tests/pgr/tncentral.fa.gz tests/genome/mg1655.fa.gz \
    > tests/pgr/mg1655.ir.json

pgr rept masker tests/pgr/tncentral.fa.gz tests/genome/mg1655.fa.gz \
    > tests/pgr/mg1655.rmask.json

pgr runlist stat tests/genome/mg1655.chr.sizes tests/pgr/mg1655.ir.json

pgr rept s-kmer tests/genome/mg1655.fa.gz \
    > tests/pgr/mg1655.rept.json

pgr rept trf tests/genome/mg1655.fa.gz \
    > tests/pgr/mg1655.trf.json

pgr rept s-align tests/genome/mg1655.fa.gz \
    > tests/pgr/mg1655.salign.json

pgr runlist stat tests/genome/mg1655.chr.sizes tests/pgr/mg1655.rm.json
pgr runlist statop tests/genome/mg1655.chr.sizes tests/pgr/mg1655.ir.json tests/pgr/mg1655.rm.json
```

---

## e-kmer

Identify repeats in a genome against an external repeat library
(Dfam, RepBase, TnCentral) by k-mer analysis, mimicking `RepeatMasker`.

### Usage

```bash
pgr rept e-kmer [OPTIONS] <repeat> <infile>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `repeat` | | | File | Repeat database FASTA (Dfam, RepBase, etc.) |
| `infile` | | | File | Input genome FASTA (`.fa.gz` supported) |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `kmer` | `-k` | `--kmer` | Int | K-mer size (default: 17) |
| `fk` | | `--fill-kmer` | Int | Fill holes between repetitive k-mers (default: 2) |
| `min` | | `--min-len` | Int | Minimum length of repetitive fragments (default: 300) |
| `ff` | | `--fill-fragment` | Int | Fill holes between repetitive fragments (default: 10) |
| `keep_index` | | `--keep-index` | Flag | Keep the built repeat table next to the library for reuse |

### Dependencies

*   `FastK`, `Profex` (from FastK suite)

### Differences from RepeatMasker

`e-kmer` is a fast k-mer-based approximation, not a full `RepeatMasker` replacement.

*   **No repeat annotation**: `RepeatMasker` classifies each hit into a repeat family/class (`repeatmasker.out`). `e-kmer` only reports genomic intervals — it never labels a region with its repeat family.
*   **k-mer sensitivity**: detection relies on k-mers shared with the repeat database. Highly diverged copies sharing few exact k-mers are missed or split into fragments; the `--fill-kmer` / `--fill-fragment` steps bridge small gaps but cannot recover long diverged copies.
*   **Intervals only**: the output is a runlist of `chr:start-end` intervals (JSON), not a masked sequence. Feed it to `pgr fa mask --runlist` to soft-mask the genome.
*   **External tools**: requires `FastK` / `Profex` in `$PATH`, plus a repeat database (Dfam, RepBase, etc.). Interval merging is done internally by `pgr runlist`.
*   **Use case**: suitable for a quick, cheap repeat-masking pass on large genomes. For annotation-grade results (family/class labels, consensus coverage), use `RepeatMasker`.

### Caching the repeat table

`e-kmer` builds a FastK table from the repeat library on every run (in a
temporary directory). Pass `--keep-index` to save that table next to the
library (`<library>.repeat.k<k>.ktab` plus hidden part files and a
`.complete` marker); later runs reuse it directly instead of rebuilding.
The cache is invalidated automatically when the library file changes
(mtime). Same convention as `pgr align pgi --keep-index`.

### Combining multiple libraries

`e-kmer` takes a single repeat library as its positional `repeat` argument.
To mask against several libraries (e.g. RepBase and Dfam) at once, run once
per library and merge the resulting runlists with `pgr runlist merge`:

```bash
pgr rept e-kmer repbase.fa.gz genome.fa -o repbase.json
pgr rept e-kmer dfam.fa.gz   genome.fa -o dfam.json

pgr runlist merge repbase.json dfam.json -o libs.json
pgr fa mask genome.fa --runlist libs.json -o masked.fa
```

Running per library keeps the `--keep-index` cache valid for each library and
lets you tune parameters (and diagnose hits) per library. The genome is
scanned once per library, which is negligible for typical repeat libraries.

## e-align

Identify repeats in a genome against an external repeat library (Dfam,
RepBase, TnCentral) by alignment, mimicking the masking behavior of
`RepeatMasker` without its annotation post-processing. The library is
aligned to the genome with `pgr align pgi` (reference = genome, query =
library); alignment blocks are filtered by identity and length, merged into
intervals, and written as a runlist JSON ready for `pgr fa mask`.

Compared with `e-kmer`, `e-align` is slower (full alignment instead of
k-mer counting) but more accurate: it reports only blocks with enough
identity, so it is the preferred choice when masking quality matters.

### Usage

```bash
pgr rept e-align [OPTIONS] <repeat> <infile>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `repeat` | | | File | Repeat database FASTA (Dfam, RepBase, etc.) |
| `infile` | | | File | Input genome FASTA (`.fa.gz` supported) |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `kmer` | `-k` | `--kmer` | Int | k-mer size for indexing (default: 31) |
| `smer` | | `--smer` | Int | Syncmer s-mer length (default: 8) |
| `window` | | `--window` | Int | Syncmer window (default: 5) |
| `freq` | `-f` | `--freq` | Int | Max k-mer frequency to keep as seed (default: 50) |
| `min_shared` | | `--min-shared` | Int | Min shared seed length (default: 12) |
| `min_identity` | | `--min-identity` | Float | Min alignment identity (default: 0.70) |
| `min_len` | | `--min-len` | Int | Min length of repetitive fragments (default: 50) |
| `fill_fragment` | | `--fill-fragment` | Int | Fill holes between fragments (default: 10) |
| `parallel` | `-p` | `--parallel` | Int | Number of threads (default: 8) |
| `keep_index` | | `--keep-index` | Flag | Keep built pgi indexes next to the inputs |

`--min-identity` uses the gap-compressed identity
(`(matches + repeat_matches) / (matches + mismatches + repeat_matches)`,
insert bases excluded) — unlike `pgr sd`, whose block identity includes
insert bases in the denominator.

### Notes

*   The input genome must not be soft-masked: lowercase (soft-masked) repeat
    regions fragment the alignment and drastically underestimate coverage.
    `e-align` warns when it detects lowercase; uppercase the genome first
    (`tr a-z A-Z`) if warned. (`e-kmer` is case-insensitive and unaffected.)
*   The defaults (k=31, freq=50, min-shared=12) are tuned for sensitivity,
    matching the empirically chosen `pgr sd` pgi parameters (TnCentral × 10
    E. coli: ~47% closer to a LASTZ `set01` baseline than the previous
    k=40/freq=100/min-shared=16 defaults). Raise k/freq/min-shared for
    speed or specificity.
*   `--keep-index` caches the pgi indexes next to the inputs for reuse, same
    convention as `pgr align pgi`.

---

## masker

Simulate RepeatMasker 4.2.4's `-lib` pipeline: per 60 kb / 2 kb batch, run
TRF PERFECT (young simple repeats, excised), the rmblastn library search
(`general_search_parameters`, parameter-for-parameter), then TRF DIVERGED
(old simple repeats) on the masked query, and merge everything into a runlist
JSON ready for `pgr fa mask`. This is the closest pgr equivalent to a real
RepeatMasker run: it covers both the library hits (IS elements etc.) and the
simple repeats, with RepeatMasker's scoring matrix, word size, gap penalties,
cutoff, and TRF parameters (notes/ecoli-repeats.md §2.8).

### Usage

```bash
pgr rept masker [OPTIONS] <repeat> <infile>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `repeat` | | | File | Repeat database FASTA (`.fa.gz` supported) |
| `infile` | | | File | Input genome FASTA (`.fa.gz` supported) |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `cutoff` | | `--cutoff` | Int | RepeatMasker cutoff score (default: 225) |
| `speed` | | `--speed` | slow/default/quick/rush | Search speed tier; sets `-word_size` 8/9/11/13 (default: default → 9) |
| `matrix_gc` | | `--matrix-gc` | Int | Fixed GC % for scoring matrix selection (default: per fragment) |
| `min_len` | | `--min-len` | Int | Min length of repetitive fragments (default: 0 = RepeatMasker raw hits) |
| `fill_fragment` | | `--fill-fragment` | Int | Fill holes between fragments (default: 0 = RepeatMasker raw hits) |
| `parallel` | `-p` | `--parallel` | Int | Total threads across rmblastn processes (default: 8; 4 per process, like RepeatMasker) |
| `frag` | | `--frag` | Int | Max fragment length before splitting (default: 60000, RepeatMasker `-frag`; 0 = whole chromosome) |
| `rmblast_dir` | | `--rmblast-dir` | Dir | Directory with makeblastdb/rmblastn (optional; falls back to `$PATH`) |

### Dependencies

*   `makeblastdb`, `rmblastn` (RMBlast ≥ 2.13; validated with the CBP build
    of 2.14.1)
*   `trf`

### Notes

*   Replicated RepeatMasker 4.2.4 `general_search_parameters`: `-gapopen 24
    -gapextend 6`, `-mask_level 101`, `-complexity_adjust`, `-dust no`,
    xdrops 450/225/112, `-num_alignments 9999999`, and the GC-keyed
    `20p##g.matrix` scoring matrix (RepeatMasker `chooseMatrices`, selected
    per fragment, RepeatMasker's 60 kb / 2 kb batching).
*   TRF stages use RepeatMasker's own parameters: PERFECT (2/7/7/80/10/50/10,
    copy > 4) then DIVERGED (2/3/5/75/20/33/7, copy > 5), with PERFECT simple
    repeats and library hits X-masked between stages (RepeatMasker excises
    them; X-masking is hit-set equivalent and keeps coordinates simple).
*   RepeatMasker's annotation post-processing (family/class, fragment
    re-joining, boundary refinement) is not replicated; intervals are the
    stage outputs merged per chromosome. `--min-len` / `--fill-fragment` can
    be raised to match `e-align`'s filtering.
*   The official NCBI prebuilt RMBlast 2.14.1 requires glibc ≥ 2.29 and will
    not run on CentOS 7 (glibc 2.17); on old systems use the CBP build
    (glibc ≤ 2.16) or bioconda's package, and point `--rmblast-dir` at it.
*   On the 10-strain E. coli cohort, RepeatMasker's full output (IS + simple
    repeats) is fully covered (100.0%); our extra ~0.6% is element-end
    flanks that RepeatMasker's boundary refinement trims.
*   Soft-masked (lowercase) genomes are warned about: rmblastn skips
    lowercase regions, so uppercase the genome first (`tr a-z A-Z`).

---

## s-kmer

Identify repetitive regions in a genome using self k-mer depth analysis;
no repeat library is needed.

### Usage

```bash
pgr rept s-kmer [OPTIONS] <infile>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infile` | | | File | Input genome FASTA (`.fa.gz` supported) |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `kmer` | `-k` | `--kmer` | Int | K-mer size (default: 17) |
| `fk` | | `--fill-kmer` | Int | Fill holes between repetitive k-mers (default: 2) |
| `min` | | `--min-len` | Int | Min length of repetitive fragments (default: 100) |
| `ff` | | `--fill-fragment` | Int | Fill holes between repetitive fragments (default: 10) |

### Dependencies

*   `FastK`, `Profex`

---

## trf

Identify tandem repeats in a genome via `trf`.

### Usage

```bash
pgr rept trf [OPTIONS] <infile>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infile` | | | File | Input genome FASTA |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `trf_match` | | `--trf-match` | Int | TRF matching weight (default: 2) |
| `trf_mismatch` | | `--trf-mismatch` | Int | TRF mismatching penalty (default: 7) |
| `delta` | | `--delta` | Int | Indel penalty (default: 7) |
| `pm` | | `--pm` | Int | Match probability (default: 80) |
| `pi` | | `--pi` | Int | Indel probability (default: 10) |
| `min_score` | | `--min-score` | Int | Min alignment score (default: 50) |
| `max_period` | | `--max-period` | Int | Max period size (default: 2000) |

### Dependencies

*   `trf`

## s-align

Identify repetitive regions of a genome by self-alignment, without any
repeat library. This is the pgr-native port of the Cactus-style pipeline of
`scripts/pgr-repeat.sh`: the genome is split into overlapping windows, the
windows are aligned back to the genome with `lastz`, lifted to genomic
coordinates, and regions whose alignment depth exceeds a threshold are kept
(with 50%-overlap windows the baseline depth is 2; `--min-depth 4` means
at least two copies).

> Soft-masked (lowercase) genomes are detected and warned about: lowercase
> repeat regions are skipped by `lastz` and underestimate coverage, so
> uppercase the genome first (`tr a-z A-Z`) if warned.

> **What it detects**: self-alignment captures **every** region that appears
> more than once in the genome — transposable elements, segmental
> duplications (SD), tandem repeats, multi-copy gene families, etc. It does
> not distinguish repeat types and does not restrict hits to transposable
> elements (unlike the library-driven `e-kmer` / `e-align`, which only report
> regions matching the repeat library). Use `s-align` when the goal is
> broad "mask everything repetitive"; use `e-kmer` / `e-align` when only
> known transposable elements should be masked.
>
> **Before an SD search**: if the masked genome is later fed to a segmental
> duplication detector (the BISER-style workflow, which expects repeats —
> but not SDs — to be masked first), do **not** use `s-align` (or `s-kmer`):
> self-comparison would mask the SDs themselves and leave nothing to find.
> Mask with interspersed-repeat detection against a library (`e-kmer` /
> `e-align`) plus tandem repeats (`trf`) instead, following the reference
> recipe of TRF + RepeatMasker preprocessing.

### Usage

```bash
pgr rept s-align [OPTIONS] <infile>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infile` | | | File | Input genome FASTA (`.fa.gz` supported) |
| `outfile` | `-o` | `--outfile` | File | Output filename (default: stdout) |
| `window` | `-w` | `--window` | Int | Overlapping window length (default: 200) |
| `step` | | `--step` | Int | Window step size (default: 100) |
| `chunk_records` | | `--chunk-records` | Int | Split window output into chunks of N records (default: 10000) |
| `preset` | | `--preset` | Str | lastz parameter set (default: set01) |
| `parallel` | `-p` | `--parallel` | Int | Number of threads (default: 4) |
| `min_depth` | `-m` | `--min-depth` | Int | Minimum alignment depth to keep a region (default: 4) |

### Dependencies

*   `lastz`

## Notes

*   `pgr rept e-kmer` accepts any repeat FASTA (Dfam, RepBase, TnCentral).
*   RepeatMasker normally uses a species-specific library via FamDB, or a
    custom library via `-lib`.
*   A native pgr masking plan (Dfam full library + `pgr align pgi` alignment
    -> runlist -> `fa mask`) is tracked in the design notes
    [notes/design/repeat-masking.md](../notes/design/repeat-masking.md).
