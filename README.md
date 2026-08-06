# pgr - Practical Genome Refiner

[![Build](https://github.com/wang-q/pgr/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/pgr/actions)
[![codecov](https://codecov.io/gh/wang-q/pgr/branch/master/graph/badge.svg)](https://codecov.io/gh/wang-q/pgr)
[![license](https://img.shields.io/github/license/wang-q/pgr)](https://github.com//wang-q/pgr)

`pgr` is a command-line toolkit for working with genomes and genome-derived data: sequences,
alignments, variation, and related formats.

It is designed as a practical “Swiss Army knife” for day-to-day bioinformatics workflows, with a
focus on:

- Format-aware utilities for common genomics file types (FASTA/FASTQ/2bit, AXT/PSL/Chain/Net/MAF,
  GFF)
- Interoperable outputs (tabular conventions, FASTA/MAF for alignments)
- Pipeline-friendly behavior (stdin/stdout where possible, predictable output, composable
  subcommands)
- Performance and robustness (Rust implementation, zero-panic policy for malformed inputs)

High-level capabilities include:

- Sequences: FASTA/FASTQ inspection, filtering, slicing, conversion, 2bit querying, and pbit
  population archive compression
- Alignments: sorting, filtering, conversion, and coordinate/range utilities across UCSC formats
- Pangenome: PAF implicit graph indexing, querying, and conversion (BED/MAF/GFA/VCF)
- Pipelines & plots: integrated workflows (optionally using external tools) and LaTeX/TikZ figure
  generation

## Install

Current release: 0.4.0

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

## Command naming conventions

`pgr` organizes commands in two levels. The naming rules make the command line
predictable:

**First-level commands are named after the input format or the task domain:**

* Input formats: `fa`, `fas`, `fq`, `2bit`, `gff`, `rg`, `axt`, `chain`,
  `net`, `maf`, `paf`, `psl`, `lav`, `ms`
* Task domains: `dist`, `sd`, `rept`, `runlist`, `pl`, `plot`, `align`,
  `pgi`, `pbit`

**Second-level commands follow one of three naming patterns:**

1. **Operations within one format** (the majority, 70+ commands):
   `fa mask/sort/dedup/filter/rc/size`, `psl lift/stats/swap`,
   `paf query/graph`, `sd align/cluster/cross`,
   `runlist span/compare/merge`. Because the input and output share the same
   format, the operation name is what distinguishes one command from another.
2. **Format conversions are named after the output**, with a uniform `to-`
   prefix (about 25 commands across 12 families):
   `to-psl`, `to-maf`, `to-fas`, `to-paf`, `to-vcf`, `to-gfa`, `to-bed`,
   `to-chain`, `to-axt`, `to-hv`, `to-fa`, `to-2bit`, `to-dna`, `to-xlsx`,
   `to-rg`. This is the project-wide rule that answers "input or output":
   conversion commands are named after the output format.
3. **A few commands are named after the artifact or the argument**:
   `gff rg`, `gff runlist` (output format, without the `to-` prefix),
   `chain net`, `psl chain` (output format), `fa masked`, `2bit masked`
   (output property), `fa range`, `2bit range`, `pbit range`,
   `runlist genome` (input argument concept), `paf graph/index` (artifact),
   `plot dot/hh/nrps/venn` (output chart type).

**Rule of thumb**: a second-level command is named by its operation when the
format does not change, by the output (`to-X`) when it crosses formats, and
by the artifact or argument when neither applies.

## Examples

This repository contains many subcommands and end-to-end workflows. Extended and curated examples
are collected in:

- docs/usage_examples.md
- docs/rept.md (repeat masking: libraries, RepeatMasker, pgr rept e-kmer/s-kmer/trf)

Below are a few quick examples to get started:

```bash
# Basic FASTA statistics
pgr fa size tests/fasta/ufasta.fa

# Block FA summary
pgr fas stat tests/fas/example.fas --outgroup

# 2bit range extraction
pgr 2bit range tests/genome/mg1655.2bit NC_000913:1-100

# Create a pbit population archive from a reference and sample assemblies
pgr pbit create -r tests/pgr/pseudocat.fa -i tests/pgr/pseudopig.fa -o tmp.pbit

# Extract a region from all samples in the archive
pgr pbit range tmp.pbit scaffold_1:1-1000 -o tmp.fa
```

## External dependencies

Some subcommands depend on external executables:

- `pgr pl ucsc` requires the UCSC kent-tools suite, including programs such as `faToTwoBit`,
  `axtChain`, `chainAntiRepeat`, `chainMergeSort`, `chainPreNet`,`chainNet`, `netSyntenic`,
  `netChainSubset`, `chainStitchId`, `netSplit`,`netToAxt`, `axtSort`, `axtToMaf`, `netFilter`,
  `netClass`, and `chainSplit`.
- `pgr rept trf` depends on `trf`.
- `pgr rept s-kmer` and `pgr rept e-kmer` depend on `FastK` and `Profex`.
- `pgr fas refine` depends on an external multiple sequence alignment tool such as
  `clustalw` (default), `muscle`, or `mafft`.

## Author

Qiang Wang [wang-q@outlook.com](mailto:wang-q@outlook.com)

## License

MIT.

Copyright by Qiang Wang.

Written by Qiang Wang [wang-q@outlook.com](mailto:wang-q@outlook.com), 2024-
