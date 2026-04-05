# pgr - Practical Genome Refiner

[![Build](https://github.com/wang-q/pgr/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/pgr/actions)
[![codecov](https://codecov.io/gh/wang-q/pgr/branch/master/graph/badge.svg)](https://codecov.io/gh/wang-q/pgr)
[![license](https://img.shields.io/github/license/wang-q/pgr)](https://github.com//wang-q/pgr)

`pgr` is a command-line toolkit for working with genomes and genome-derived
data: sequences, alignments, variation, phylogenies, and related formats.

It is designed as a practical “Swiss Army knife” for day-to-day bioinformatics
workflows, with a focus on:

- Format-aware utilities for common genomics file types (FASTA/FASTQ/2bit, AXT/PSL/Chain/Net/MAF, GFF, Newick)
- Interoperable outputs (tabular `cluster` / `pair` conventions, Newick for trees)
- Pipeline-friendly behavior (stdin/stdout where possible, predictable output, composable subcommands)
- Performance and robustness (Rust implementation, zero-panic policy for malformed inputs)

High-level capabilities include:

- Sequences: FASTA/FASTQ inspection, filtering, slicing, conversion, and 2bit querying
- Alignments: sorting, filtering, conversion, and coordinate/range utilities across UCSC formats
- Clustering & trees: distance/matrix processing, multiple clustering algorithms, tree cutting and visualization
- Pipelines & plots: integrated workflows (optionally using external tools) and LaTeX/TikZ figure generation

## Install

Current release: 0.2.0

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

## Author

Qiang Wang <wang-q@outlook.com>

## License

MIT.

Copyright by Qiang Wang.

Written by Qiang Wang <wang-q@outlook.com>, 2024-
