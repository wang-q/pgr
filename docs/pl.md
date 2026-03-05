# pgr pl

Integrated pipelines for genomic analysis.

`pgr pl` provides high-level workflows that combine multiple tools (both internal `pgr` commands and external binaries) to perform complex tasks like repeat masking, multiple sequence alignment construction, and UCSC-style chain/net processing.

## Subcommands

| Subcommand | Description |
| :--- | :--- |
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
| `fk` | | `--fk` | Int | Fill holes between repetitive k-mers (default: 2) |
| `min` | | `--min` | Int | Minimum length of repetitive fragments (default: 300) |
| `ff` | | `--ff` | Int | Fill holes between repetitive fragments (default: 10) |

### Dependencies

*   `FastK`, `Profex` (from FastK suite)
*   `spanr`

---

## p2m

Pairwise to Multiple (p2m) Pipeline. Constructs a "core" Multiple Sequence Alignment (MSA) from multiple pairwise alignment files (Block FASTA).

### Usage

```bash
pgr pl p2m [OPTIONS] <infiles>...
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infiles` | | | Files | Input Block FASTA files (2 or more) |
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
pgr pl prefilter [OPTIONS] <infile> <match>
```

### Arguments

| Argument | Short | Long | Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| `infile` | | | File | Input genome/metagenome FASTA |
| `match` | | | File | Reference protein FASTA |
| `chunk` | `-c` | `--chunk` | Int | Chunk size in bytes (default: 100000) |
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
| `fk` | | `--fk` | Int | Fill holes between repetitive k-mers (default: 2) |
| `min` | | `--min` | Int | Min length of repetitive fragments (default: 100) |
| `ff` | | `--ff` | Int | Fill holes between repetitive fragments (default: 10) |

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
| `match` | | `--match` | Int | Matching weight (default: 2) |
| `mismatch` | | `--mismatch` | Int | Mismatching penalty (default: 7) |
| `delta` | | `--delta` | Int | Indel penalty (default: 7) |
| `pm` | | `--pm` | Int | Match probability (default: 80) |
| `pi` | | `--pi` | Int | Indel probability (default: 10) |
| `minscore` | | `--minscore` | Int | Min alignment score (default: 50) |
| `maxperiod` | | `--maxperiod` | Int | Max period size (default: 2000) |

### Dependencies

*   `trf`
*   `spanr`

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
| `lineargap` | | `--lineargap` | Str | Gap cost: "loose" (default) or "medium" |
| `minscore` | | `--minscore` | Int | Min alignment score (default: 1000) |
| `tname` | `-t` | `--tname` | Str | Custom target name prefix |
| `qname` | `-q` | `--qname` | Str | Custom query name prefix |
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
