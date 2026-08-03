# pgr rept

Repeat detection and masking. `pgr rept` provides commands that detect
repetitive regions of a genome and emit runlist JSON ready for
`pgr fa mask`.

Command names use two dimensions:

*   prefix `e` = external repeat library required, `s` = self (genome only);
*   suffix `kmer` / `align` = detection mechanism (`align` variants are planned).

## Repeat libraries

Repeat-masking workflows compare a genome against a library of known repeat
consensus sequences.

### TnCentral

TnCentral is a database of prokaryotic insertion sequences. Download the
complete IS library and prepare it for `pgr rept e-kmer`:

```bash
# Download and unpack
curl -LO https://tncentral.ncc.unesp.br/api/download_blast/nc/tn_in_is
unzip -j tn_in_is 'tncentral_integrall_isfinder.fa'
gzip -9 -c 'tncentral_integrall_isfinder.fa' > tncentral.fa.gz

# Sanity check and quality filter
pgr fa size tests/pgr/tncentral.fa.gz
pgr dist seq tests/pgr/tncentral.fa.gz -k 17 -w 5 -p 8 |
    spanr filter stdin --ge 5:0.9
```

### RepBase

RepBase is the classic repeat database used by RepeatMasker. The distribution
tarball contains an EMBL-format library; convert it to FASTA with readseq:

```bash
curl -LO https://github.com/wang-q/ubuntu/releases/download/20190906/repeatmaskerlibraries-20140131.tar.gz
tar xvfz repeatmaskerlibraries-20140131.tar.gz Libraries/RepeatMaskerLib.embl

# https://sourceforge.net/projects/readseq/
java -jar ~/bin/readseq.jar -f fa Libraries/RepeatMaskerLib.embl
mv Libraries/RepeatMaskerLib.embl.fasta repbase.fa
gzip -9 -k repbase.fa
```

### Dfam

Dfam is a curated database of transposable element families. Its full
consensus FASTA can be used directly as the candidate set for masking
workflows (`pgr rept e-kmer`, or a `pgr align pgi`-based masker):

*   Dfam website: <https://www.dfam.org/>
*   Dfam releases / downloads: <https://www.dfam.org/releases/>
*   The full consensus FASTA from the releases page works with `pgr rept e-kmer`
    and with RepeatMasker's `-lib` option.

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
    [notes/design/repeatmasker.md](../notes/design/repeatmasker.md).
