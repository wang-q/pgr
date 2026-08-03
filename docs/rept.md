# pgr rept

Repeat detection and masking. `pgr rept` provides commands that detect
repetitive regions of a genome and emit runlist JSON ready for
`pgr fa mask`.

Command names use two dimensions:

*   prefix `e` = external repeat library required, `s` = self (genome only);
*   suffix `kmer` / `align` = detection mechanism (`align` variants are planned).

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
pgr dist seq tncentral.fa.gz -k 17 -w 5 -p 8 |
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
pgr dist seq repbase.fa.gz -k 17 -w 5 -p 8 |
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
pgr dist seq dfam.fa.gz -k 17 -w 5 -p 8 |
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
spanr stat mg1655.sizes tncentral.json
spanr statop mg1655.sizes repbase.json dfam.json
```

Measured on MG1655 (RepeatMasker reference output: 49,379 bp, 1.06%):

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
singularity image, converting its `.out` to a GFF runlist for comparison:

```bash
singularity run ~/bin/repeatmasker_master.sif /app/RepeatMasker/RepeatMasker \
    ./genome.fa -xsmall -species "bacteria"

singularity run ~/bin/repeatmasker_master.sif /app/RepeatMasker/util/rmOutToGFF3.pl \
    ./genome.fa.out > mg1655.rm.gff

spanr gff tests/pgr/mg1655.rm.gff -o tests/pgr/mg1655.rm.json
```

## Subcommands

| Subcommand | Description |
| :--- | :--- |
| `e-kmer` | Identify repeats against an external library (k-mer) |
| `s-kmer` | Identify repetitive regions by self k-mer depth (no library) |
| `trf` | Identify tandem repeats via `trf` |

All three emit runlist JSON ready for `pgr fa mask`:

```bash
pgr rept e-kmer tests/pgr/tncentral.fa.gz tests/pgr/mg1655.fa.gz \
    > tests/pgr/mg1655.ir.json

spanr stat tests/pgr/mg1655.chr.sizes tests/pgr/mg1655.ir.json

pgr rept s-kmer tests/pgr/mg1655.fa.gz \
    > tests/pgr/mg1655.rept.json

pgr rept trf tests/pgr/mg1655.fa.gz \
    > tests/pgr/mg1655.trf.json

spanr stat tests/pgr/mg1655.chr.sizes tests/pgr/mg1655.rm.json
spanr statop tests/pgr/mg1655.chr.sizes tests/pgr/mg1655.ir.json tests/pgr/mg1655.rm.json
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
*   `spanr`

### Differences from RepeatMasker

`e-kmer` is a fast k-mer-based approximation, not a full `RepeatMasker` replacement.

*   **No repeat annotation**: `RepeatMasker` classifies each hit into a repeat family/class (`repeatmasker.out`). `e-kmer` only reports genomic intervals — it never labels a region with its repeat family.
*   **k-mer sensitivity**: detection relies on k-mers shared with the repeat database. Highly diverged copies sharing few exact k-mers are missed or split into fragments; the `--fill-kmer` / `--fill-fragment` steps bridge small gaps but cannot recover long diverged copies.
*   **Intervals only**: the output is a runlist of `chr:start-end` intervals (JSON), not a masked sequence. Feed it to `pgr fa mask --runlist` to soft-mask the genome.
*   **External tools**: requires `FastK` / `Profex` and `spanr` in `$PATH`, plus a repeat database (Dfam, RepBase, etc.).
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
per library and merge the resulting runlists with `spanr merge`:

```bash
pgr rept e-kmer repbase.fa.gz genome.fa -o repbase.json
pgr rept e-kmer dfam.fa.gz   genome.fa -o dfam.json

spanr merge repbase.json dfam.json -o libs.json
pgr fa mask genome.fa --runlist libs.json -o masked.fa
```

Running per library keeps the `--keep-index` cache valid for each library and
lets you tune parameters (and diagnose hits) per library. The genome is
scanned once per library, which is negligible for typical repeat libraries.

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
*   `spanr`

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
*   `spanr`

## Notes

*   `pgr rept e-kmer` accepts any repeat FASTA (Dfam, RepBase, TnCentral).
*   RepeatMasker normally uses a species-specific library via FamDB, or a
    custom library via `-lib`.
*   A native pgr masking plan (Dfam full library + `pgr align pgi` alignment
    -> runlist -> `fa mask`) is tracked in the design notes
    [notes/design/repeat-masking.md](../notes/design/repeat-masking.md).
