# pgr - Practical Genome Refiner

[![Build](https://github.com/wang-q/pgr/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/pgr/actions)
[![codecov](https://codecov.io/gh/wang-q/pgr/branch/master/graph/badge.svg)](https://codecov.io/gh/wang-q/pgr)
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
