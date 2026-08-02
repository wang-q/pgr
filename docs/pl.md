# pgr pl

Integrated pipelines for genomic analysis.

`pgr pl` provides high-level workflows that combine multiple tools (both internal `pgr` commands and external binaries) to perform complex tasks like repeat masking, multiple sequence alignment construction, and UCSC-style chain/net processing.

## Subcommands

| Subcommand | Description |
| :--- | :--- |
| `chainnet` | Native chain/net pipeline (psl -> chain -> net -> axt -> maf, no kent-tools) |
| `ir` | Identify interspersed repeats (RepeatMasker-like) |
| `p2m` | Pairwise to Multiple alignment pipeline |
| `prefilter` | Prefilter genome/metagenome by amino acid minimizers |
| `rept` | Identify repetitive regions using k-mer analysis |
| `trf` | Identify tandem repeats via `trf` |
| `ucsc` | UCSC chain/net pipeline (psl -> chain -> net -> maf) |

---

## ir

Identify interspersed repeats in a genome. This command mimics the functionality of `RepeatMasker` by using k-mer analysis against a repeat database.

### Usage

```bash
pgr pl ir [OPTIONS] <repeat> <infile>
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

### Dependencies

*   `FastK`, `Profex` (from FastK suite)
*   `spanr`

### Differences from RepeatMasker

`ir` is a fast k-mer-based approximation, not a full `RepeatMasker` replacement.

*   **No repeat annotation**: `RepeatMasker` classifies each hit into a repeat family/class (`repeatmasker.out`). `ir` only reports genomic intervals — it never labels a region with its repeat family.
*   **k-mer sensitivity**: detection relies on k-mers shared with the repeat database. Highly diverged copies sharing few exact k-mers are missed or split into fragments; the `--fill-kmer` / `--fill-fragment` steps bridge small gaps but cannot recover long diverged copies.
*   **Intervals only**: the output is a runlist of `chr:start-end` intervals (JSON), not a masked sequence. Feed it to `pgr fa mask --runlist` to soft-mask the genome.
*   **External tools**: requires `FastK` / `Profex` and `spanr` in `$PATH`, plus a repeat database (Dfam, RepBase, etc.) for `ir`.
*   **Use case**: suitable for a quick, cheap repeat-masking pass on large genomes. For annotation-grade results (family/class labels, consensus coverage), use `RepeatMasker`.

---

## p2m

Pairwise to Multiple (p2m) Pipeline. Constructs a "core" Multiple Sequence Alignment (MSA) from multiple pairwise alignment files (Block FA).

### Usage

```bash
pgr pl p2m [OPTIONS] <infiles>...
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infiles` | | | Files | Input Block FA files (2 or more) |
| `outdir` | `-o` | `--outdir` | Dir | Output directory (default: "PL-p2m") |

### Logic

1.  **Reference-Based**: The first species of the first input file is treated as the reference target.
2.  **Intersection**: Only genomic regions covered by *all* input files are retained.
3.  **Stitching**: Aligned sequences are sliced and joined to form a gap-free core alignment.

### Dependencies

*   `spanr`

---

## prefilter

Prefilter genome/metagenome assembly by amino acid minimizers. Filters sequences by comparing them against protein references.

### Usage

```bash
pgr pl prefilter [OPTIONS] <infile> <reference>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infile` | | | File | Input genome/metagenome FASTA |
| `reference` | | | File | Reference protein FASTA |
| `chunk` | `-c` | `--chunk-size` | Int | Chunk size in bytes (default: 100000) |
| `len` | | `--len` | Int | Min amino acid length (default: 15) |
| `kmer` | `-k` | `--kmer` | Int | K-mer size (default: 7) |
| `window` | `-w` | `--window` | Int | Window size (default: 1) |
| `parallel` | `-p` | `--parallel` | Int | Number of threads (default: 1) |

---

## rept

Identify repetitive regions in a genome using k-mer analysis (self-comparison).

### Usage

```bash
pgr pl rept [OPTIONS] <infile>
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
pgr pl trf [OPTIONS] <infile>
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

---

## chainnet

Native chain/net pipeline. Runs the pairwise genome alignment workflow
(psl -> chain -> net -> axt -> maf) entirely with `pgr` commands — no external
kent-tools required. Output has been verified byte-for-byte identical to the
UCSC kent-tools pipeline (`pgr pl ucsc`) for all intermediate files.

### Usage

```bash
pgr pl chainnet [OPTIONS] <target> <query> <psl>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `target` | | | File | Target (reference) genome FASTA |
| `query` | | | File | Query genome FASTA |
| `psl` | | | Path | PSL file or directory containing PSL files |
| `outdir` | `-o` | `--outdir` | Dir | Output directory (default: stdout) |
| `gap_model` | | `--gap-model` | Str | Linear gap cost: "loose" (default) or "medium" |
| `min_score` | | `--min-score` | Int | Min alignment score (default: 1000) |
| `tname` | | `--t-name` | Str | Custom target name prefix |
| `qname` | | `--q-name` | Str | Custom query name prefix |
| `syn` | | `--syn` | Flag | Generate syntenic alignments only |

### Pipeline Steps

1.  **Prep**: `pgr fa size` + `pgr fa to-2bit` for target and query.
2.  **Chain**: `pgr psl chain` + `pgr chain anti-repeat`.
3.  **Merge**: `pgr chain sort`.
4.  **PreNet**: `pgr chain pre-net`.
5.  **Net**: `pgr chain net` + `pgr net syntenic` + `pgr net subset` +
    `pgr chain stitch` + `pgr net split`.
6.  **Axt**: `pgr net to-axt` | `pgr axt sort`.
7.  **Maf**: `pgr axt to-maf`.

### Dependencies

None external — uses only the `pgr` binary itself.

---

## ucsc

UCSC chain/net pipeline. Converts PSL alignments to MAF format via Chain/Net processing.

### Usage

```bash
pgr pl ucsc [OPTIONS] <target> <query> <psl>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `target` | | | File | Target (reference) genome FASTA |
| `query` | | | File | Query genome FASTA |
| `psl` | | | Path | PSL file or directory containing PSL files |
| `outdir` | `-o` | `--outdir` | Dir | Output directory (default: stdout) |
| `gap_model` | | `--gap-model` | Str | Gap cost: "loose" (default) or "medium" |
| `min_score` | | `--min-score` | Int | Min alignment score (default: 1000) |
| `tname` | | `--t-name` | Str | Custom target name prefix |
| `qname` | | `--q-name` | Str | Custom query name prefix |
| `syn` | | `--syn` | Flag | Generate syntenic alignments only |

### Pipeline Steps

1.  **axtChain**: Chain together alignments.
2.  **chainMergeSort**: Merge sorted chains.
3.  **chainPreNet**: Remove chains unlikely to be netted.
4.  **chainNet**: Create alignment nets.
5.  **netSyntenic**: Add synteny information.
6.  **netToAxt**: Convert net to AXT.
7.  **axtToMaf**: Convert AXT to MAF.

### Dependencies

Requires UCSC Kent tools in `$PATH`:
`axtChain`, `chainAntiRepeat`, `chainMergeSort`, `chainPreNet`, `chainNet`, `netSyntenic`, `netChainSubset`, `chainStitchId`, `netSplit`, `netToAxt`, `axtSort`, `axtToMaf`, `netFilter`, `chainSplit`.
