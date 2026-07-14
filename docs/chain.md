# pgr chain

The `pgr chain` module provides tools for manipulating **UCSC Chain** format alignment files. These tools are core components of whole-genome alignment pipelines such as `pgr pl ucsc`.

## Overview

- **Role**: High-level processing and filtering of Chain alignments.
- **Input**: Chain format files (text).
- **Output**: Processed Chain files or Net files.
- **Complements**:
  - Upstream: `lastz` -> `axtChain` (or `pgr psl chain`) to generate chains.
  - Downstream: `pgr net` (to generate Net), `pgr maf` (to generate MAF).

## Subcommands

### 1. `pgr chain sort`: Sort chains

Sorts a Chain file by score in descending order.

- **Purpose**: Tools such as `chainPreNet` and `chainNet` usually require score-sorted input.
- **Arguments**:
  - `infiles`: Input Chain file(s).
  - `--input-list`: File containing a list of input paths (one per line).
  - `--save-id`: Preserve original Chain IDs (default: renumber from 1).

### 2. `pgr chain split`: Split chains

Splits a Chain file into separate files by target or query sequence name.

- **Purpose**: Parallel processing or splitting a large file into chromosome-organized files.
- **Arguments**:
  - `infiles`: Input Chain file(s).
  - `-o, --outdir <dir>`: Output directory (required).
  - `--by-query`: Split by query sequence name (default: target).
  - `--lump <N>`: Group results into at most N files by hashing sequence names, avoiding many small files.

### 3. `pgr chain stitch`: Stitch chain fragments

Joins chain fragments sharing the same Chain ID into a single chain per ID.

- **Purpose**: Repair cases where the same chain was broken by parallel processing or file splitting.
- **Behavior**: Fragments are grouped by ID, checked for consistent target/query name and query strand, converted to blocks, sorted by target start, and rebuilt into one chain. Scores are summed.
- **Arguments**:
  - `infile`: Input Chain file.
  - `-o, --outfile <file>`: Output Chain file.

### 4. `pgr chain anti-repeat`: Repeat and degeneracy filter

Filters out chains composed mainly of repetitive or low-complexity sequence.

- **Purpose**: Improve alignment quality by removing biologically meaningless false-positive alignments.
- **Mechanism**:
  - **Degeneracy filter**: Checks whether the alignment is mostly low-complexity sequence (for example, `ATATAT...`).
  - **Repeat filter**: Checks whether the alignment falls in soft-masked (lowercase) regions.
- **Arguments**:
  - `target_2bit`: Target genome 2bit file.
  - `query_2bit`: Query genome 2bit file.
  - `infile`: Input Chain file.
  - `--min-score`: Minimum score threshold (default: 5000).
  - `--no-check-score`: Chains above this score skip checks (default: 200000).
  - `-o, --outfile <file>`: Output Chain file.

### 5. `pgr chain pre-net`: Pre-net filtering

Removes chains that are fully covered by higher-scoring chains and therefore cannot contribute to a net.

- **Purpose**: Significantly reduce Chain file size and speed up the subsequent `chainNet` step.
- **Mechanism**: Uses a bitmap to track coverage on target and query sequences. Higher-scoring chains are processed first; lower-scoring chains whose blocks are already fully covered are dropped.
- **Arguments**:
  - `infile`: Input Chain file.
  - `t_sizes`: Target chromosome sizes file.
  - `q_sizes`: Query chromosome sizes file.
  - `--pad`: Padding around blocks (default: 1).
  - `--dots <N>`: Print a progress dot every N chains.
  - `--incl-hap`: Include haplotype sequences (`_hap`, `_alt`).
  - `-o, --outfile <file>`: Output Chain file.

### 6. `pgr chain net`: Build nets

Converts a Chain file into Net format (syntenic nets).

- **Purpose**: Net format represents high-level correspondences between genomes, distinguishing orthologs and paralogs and handling inversions and translocations.
- **Output**: Two Net files, one in target orientation and one in query orientation.
- **Arguments**:
  - `infile`: Input Chain file (should be sorted by score descending).
  - `t_sizes`: Target chromosome sizes file.
  - `q_sizes`: Query chromosome sizes file.
  - `out_target_net`: Output target Net file.
  - `out_query_net`: Output query Net file.
  - `--min-space`: Minimum gap size to fill (default: 25).
  - `--min-fill`: Minimum fill to record (default: `min-space / 2`).
  - `--min-score`: Minimum Chain score threshold (default: 2000).
  - `--incl-hap`: Include haplotype sequences (`_hap`, `_alt`).

## Typical workflow (UCSC pipeline)

```bash
# 1. Sort
pgr chain sort raw.chain -o sorted.chain

# 2. Pre-net filtering - remove redundant chains covered by higher-scoring alignments
pgr chain pre-net sorted.chain t.sizes q.sizes -o pre.chain

# 3. Build nets
pgr chain net pre.chain t.sizes q.sizes t.net q.net

# 4. Add synteny information (optional, usually with pgr net syntenic)
# ...
```
