# pgr - Practical Genome Refiner

[![Build](https://github.com/wang-q/pgr/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/pgr/actions)
[![codecov](https://codecov.io/gh/wang-q/pgr/branch/master/graph/badge.svg?token=8toyNHCsVU)](https://codecov.io/gh/wang-q/pgr)
[![license](https://img.shields.io/github/license/wang-q/pgr)](https://github.com//wang-q/pgr)

`pgr` is a command-line toolkit for working with genomes and genome-derived
data: sequences, alignments, variation, phylogenies, and related formats.

<!-- TOC -->
* [pgr - Practical Genome Refiner](#pgr---practical-genome-refiner)
  * [Install](#install)
  * [Usage](#usage)
  * [Synopsis](#synopsis)
    * [`pgr help`](#pgr-help)
  * [Examples](#examples)
<!-- TOC -->

## Install

Current release: 0.1.0

```bash
cargo install --path . --force #--offline

# test
cargo test -- --test-threads=1
```

## Usage

After installation, the `pgr` binary should be available in your `PATH`:

```bash
pgr help
pgr fa --help
pgr fas --help
```

## Synopsis

### `pgr help`

```text
`pgr` - Practical Genome Refiner

Usage: pgr [COMMAND]

Commands:
  ms        Hudson's ms simulator tools
  axt       Manipulate AXT alignment files
  chain     Manipulate Chain alignment files
  chaining  Chaining alignment blocks
  clust     Clustering operations
  dist      Distance/Similarity metrics
  lav       Manipulate LAV alignment files
  maf       Manipulate MAF alignment files
  mat       Matrix operations
  net       Manipulate Net alignment files
  psl       Manipulate PSL alignment files
  pl        Run integrated pipelines
  2bit      Manage 2bit files
  fa        Manipulate FASTA files
  fas       Manipulate block FA files
  fq        Manipulate FASTQ files
  help      Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

Subcommand groups:

* Simulation:
    * ms    - Hudson's ms simulator tools: to-dna

* Sequences:
    * 2bit - 2bit query and extraction
    * fa   - FASTA operations: info, records, transform, indexing
    * fas  - Block FA operations: info, subset, transform, file, variation
    * fq   - FASTQ interleaving and conversion

* Genome alignments:
    * chaining - Chaining alignments: psl
    * chain - Chain operations: sort, filter, transform, to-net
    * net   - Net operations: info, subset, transform, convert
    * axt   - AXT sorting and conversion
    * lav   - Convert to PSL
    * maf   - Convert to Block FA
    * psl   - PSL statistics, manipulation, and conversion

* Clustering:
    * clust - Algorithms: cc, dbscan, k-medoids, mcl

* Distance:
    * dist  - Metrics: hv

* Matrix:
    * mat   - Processing: compare, format, subset, to-pair, to-phylip

* Pipelines:
    * pl - Workflows: p2m, trf, ir, rept, ucsc

```

## Examples

This repository contains many subcommands and end-to-end workflows. Extended
and curated examples are collected in:

- docs/usage_examples.md

Below are a few quick examples to get started:

```bash
# Basic FASTA statistics
pgr fa size tests/fasta/ufasta.fa

# Block FA summary
pgr fas stat tests/fas/example.fas --outgroup

# 2bit range extraction
pgr 2bit range tests/genome/mg1655.2bit NC_000913:1-100
```

## External dependencies

Some subcommands depend on external executables:

- `pgr pl ucsc` requires the UCSC kent-tools suite, including programs such as
  `faToTwoBit`, `axtChain`, `chainAntiRepeat`, `chainMergeSort`, `chainPreNet`,
  `chainNet`, `netSyntenic`, `netChainSubset`, `chainStitchId`, `netSplit`,
  `netToAxt`, `axtSort`, `axtToMaf`, `netFilter`, `netClass`, and `chainSplit`.
- `pgr pl trf` depends on `trf` and `spanr`.
- `pgr pl rept` and `pgr pl ir` depend on `FastK`, `Profex`, and `spanr`.
- `pgr pl p2m` depends on `spanr`.
- `pgr fas refine` depends on an external multiple sequence alignment tool such as
  `clustalw` (default), `muscle`, or `mafft`.

## Help text style guide

* **`about`**: Third-person singular (e.g., "Counts...", "Calculates...").
* **`after_help`**: Uses raw string `r###"..."###`.
    * **Description**: Detailed explanation.
    * **Notes**: Bullet points starting with `*`.
        * Standard note for `fa`/`fas`: `* Supports both plain text and gzipped (.gz) files`
        * Standard note for `fa`/`fas`: `* Reads from stdin if input file is 'stdin'`
        * Standard note for `twobit`: `* 2bit files are binary and require random access (seeking)`
        * Standard note for `twobit`: `* Does not support stdin or gzipped inputs`
    * **Examples**: Numbered list (`1.`, `2.`) with code blocks indented by 3 spaces.
* **Arguments**:
    * **Input**: `infiles` (multiple) or `infile` (single).
        * Help: `Input [FASTA|block FA|2bit] file(s) to process`.
    * **Output**: `outfile` (`-o`, `--outfile`).
        * Help: `Output filename. [stdout] for screen`.
* **Terminology**:
    * `pgr fa` -> "FASTA"
    * `pgr fas` -> "block FA"
    * `pgr twobit` -> "2bit"

## Author

Qiang Wang <wang-q@outlook.com>

## License

MIT.

Copyright by Qiang Wang.

Written by Qiang Wang <wang-q@outlook.com>, 2024-
