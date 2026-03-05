# pgr fas

`pgr fas` provides a comprehensive suite of tools for manipulating **Block FASTA** files. Block FASTA is a format used to represent multiple sequence alignments (MSA), where each "block" consists of aligned sequences from different species or genomic regions.

## Subcommands

The subcommands are organized into the following categories:

*   **Info**: Extract information or statistics from block FA files.
    *   `check`: Check genome locations against a `chrom.sizes` file.
    *   `cover`: Calculate covered regions on chromosomes.
    *   `link`: Extract bi-lateral or multi-lateral range links.
    *   `name`: List species names present in the files.
    *   `stat`: Calculate alignment statistics (length, differences, etc.).
*   **Subset**: Filter and extract specific parts of the data.
    *   `filter`: Filter blocks by species presence or sequence length.
    *   `slice`: Extract specific alignment slices using a runlist.
    *   `subset`: Extract a subset of species from blocks.
*   **Transform**: Modify or combine block FA files.
    *   `concat`: Concatenate sequence pieces of the same species.
    *   `consensus`: Generate consensus sequences using POA (Partial Order Alignment).
    *   `join`: Join multiple files based on a common target sequence.
    *   `multiz`: Merge block FA files using a multiz-like banded DP algorithm.
    *   `refine`: Realign sequences within blocks using built-in or external tools.
    *   `replace`: Replace sequence headers using a mapping file.
*   **File**: Create or split block FA files.
    *   `create`: Create block FA files from range links.
    *   `separate`: Separate blocks into individual files per species.
    *   `split`: Split blocks into per-alignment or per-chromosome files.
*   **Variation**: Call variants from alignments.
    *   `to-vcf`: Export substitutions to VCF format.
    *   `to-xlsx`: Export substitutions and indels to an Excel file.
    *   `variation`: List variations (substitutions/indels) in TSV format.

---

## Info Commands

### check

Checks if the genome locations specified in the block headers are valid against a reference genome.

```bash
pgr fas check [OPTIONS] --genome <genome.fa> <infiles>...
```

*   `-r, --genome <path>`: Path to the reference genome FA file (required).
*   `--name <name>`: Check sequences for a specific species only.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### cover

Outputs the regions on chromosomes covered by the alignments in JSON format.

```bash
pgr fas cover [OPTIONS] <infiles>...
```

*   `--name <name>`: Only output regions for this species.
*   `--trim <int>`: Trim alignment borders by N bases to avoid overlaps (useful for lastz results).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### link

Outputs links between ranges (genomic coordinates) found in the alignment blocks.

```bash
pgr fas link [OPTIONS] <infiles>...
```

*   `--pair`: Output bilateral (pairwise) links.
*   `--best`: Output best-to-best bilateral links based on sequence similarity.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### name

Extracts all species names found in the block FA files.

```bash
pgr fas name [OPTIONS] <infiles>...
```

*   `-c, --count`: Also output the number of occurrences of each species name.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### stat

Calculates basic statistics for each alignment block (length, comparable bases, differences, gaps, etc.).

```bash
pgr fas stat [OPTIONS] <infiles>...
```

*   `--outgroup`: Treat the last sequence in each block as an outgroup (excludes it from some stats).
*   `-o, --outfile <file>`: Output filename (default: stdout).

---

## Subset Commands

### filter

Filters blocks based on species presence and sequence length, and optionally formats sequences.

```bash
pgr fas filter [OPTIONS] <infiles>...
```

*   `--name <name>`: Keep blocks containing this species.
*   `--ge <int>`: Keep sequences with length >= this value.
*   `--le <int>`: Keep sequences with length <= this value.
*   `--upper`: Convert sequences to uppercase.
*   `--dash`: Remove dashes (gaps) from sequences.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### slice

Extracts specific slices of alignments based on a provided runlist (JSON).

```bash
pgr fas slice [OPTIONS] --required <runlist.json> <infiles>...
```

*   `-r, --required <file>`: JSON file describing ranges to extract (required).
*   `--name <name>`: Reference species name (default: first species).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### subset

Extracts a subset of species from the alignment blocks.

```bash
pgr fas subset [OPTIONS] --required <name.lst> <infiles>...
```

*   `-r, --required <file>`: File with a list of species names to keep, one per line (required).
*   `--strict`: Skip blocks that do not contain *all* the required names.
*   `-o, --outfile <file>`: Output filename (default: stdout).

---

## Transform Commands

### concat

Concatenates sequence pieces of the same species from multiple blocks.

```bash
pgr fas concat [OPTIONS] --required <name.lst> <infiles>...
```

*   `-r, --required <file>`: File with a list of species names to keep/order (required).
*   `--phylip`: Output in relaxed PHYLIP format instead of FASTA.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### consensus

Generates a consensus sequence for each block using Partial Order Alignment (POA).

```bash
pgr fas consensus [OPTIONS] <infiles>...
```

*   `--engine <builtin|spoa>`: POA engine to use (default: builtin).
*   `--match <int>`: Score for matching bases (default: 5).
*   `--mismatch <int>`: Score for mismatching bases (default: -4).
*   `--gap-open <int>`: Gap opening penalty (default: -8).
*   `--gap-extend <int>`: Gap extension penalty (default: -6).
*   `--algorithm <local|global|semi_global>`: Alignment mode (default: global).
*   `--cname <name>`: Name for the consensus sequence (default: consensus).
*   `--outgroup`: Indicates the last sequence is an outgroup.
*   `-p, --parallel <int>`: Number of threads (default: 1).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### join

Joins multiple block FA files by a common target sequence.

```bash
pgr fas join [OPTIONS] <infiles>...
```

*   `--name <name>`: Target species name (default: first species).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### multiz

Merges multiple block FA files in a shared reference coordinate system using a multiz-like banded DP algorithm.

```bash
pgr fas multiz [OPTIONS] --ref <name> <infiles>...
```

*   `-r, --ref <name>`: Reference sequence name present in all inputs (required).
*   `--radius <int>`: Banded DP radius around reference diagonal (default: 30).
*   `--min-width <int>`: Minimum window width to merge (default: 1).
*   `--mode <core|union>`: Merge mode (default: core).
*   `--score-matrix <file>`: Score matrix file (LASTZ format).
*   `--gap-model <constant|medium|loose>`: Gap model (default: medium).
*   `--gap-open <int>`: Gap open cost.
*   `--gap-extend <int>`: Gap extension cost.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### refine

Realigns sequences within blocks using built-in or external MSA tools.

```bash
pgr fas refine [OPTIONS] <infiles>...
```

*   `--msa <program>`: Aligning program: `builtin` (default), `clustalw`, `mafft`, `muscle`, `spoa`, `none`.
*   `--outgroup`: Indicates presence of outgroups.
*   `--chop <int>`: Chop head and tail indels.
*   `--quick`: Quick mode, only aligns indel-adjacent regions.
*   `--pad <int>`: In quick mode, enlarge indel regions (default: 50).
*   `--fill <int>`: In quick mode, fill holes between indels (default: 50).
*   `-p, --parallel <int>`: Number of threads (default: 1).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### replace

Replaces headers in block FA files using a mapping file.

```bash
pgr fas replace [OPTIONS] --required <replace.tsv> <infiles>...
```

*   `-r, --required <file>`: TSV file containing replacement rules (required).
*   `-o, --outfile <file>`: Output filename (default: stdout).

---

## File Commands

### create

Creates block FA files from links of ranges (e.g., from `pgr fas link`).

```bash
pgr fas create [OPTIONS] --genome <genome.fa> <infiles>...
```

*   `-r, --genome <file>`: Path to the reference genome FA file (required).
*   `--name <name>`: Set a species name for ranges (if not multi-genome).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### separate

Separates block FA files into individual files per species.

```bash
pgr fas separate [OPTIONS] <infiles>...
```

*   `-s, --suffix <string>`: File extension for output files (default: .fasta).
*   `--rc`: Reverse-complement sequences if the strand is '-'.
*   `-o, --outdir <dir>`: Output directory (default: stdout).

### split

Splits block FA files into per-alignment or per-chromosome files.

```bash
pgr fas split [OPTIONS] <infiles>...
```

*   `--chr`: Split files by chromosomes.
*   `--simple`: Simplify headers by keeping only species names.
*   `-s, --suffix <string>`: File extension for output files (default: .fas).
*   `-o, --outdir <dir>`: Output directory (default: stdout).

---

## Variation Commands

### to-vcf

Exports substitutions (SNPs) to VCF format.

```bash
pgr fas to-vcf [OPTIONS] <infiles>...
```

*   `--sizes <file>`: Chrom sizes file to emit `##contig` headers.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### to-xlsx

Exports variations (substitutions and indels) to an Excel file with formatting.

```bash
pgr fas to-xlsx [OPTIONS] <infiles>...
```

*   `--indel`: Include indels.
*   `--outgroup`: Indicates presence of outgroups.
*   `--nosingle`: Omit singleton variations.
*   `--nocomplex`: Omit complex variations.
*   `--min <float>`: Minimal frequency.
*   `--max <float>`: Maximal frequency.
*   `--wrap <int>`: Wrap length for visualization (default: 50).
*   `-o, --outfile <file>`: Output filename (default: variations.xlsx).

### variation

Lists variations (substitutions and indels) in TSV format.

```bash
pgr fas variation [OPTIONS] <infiles>...
```

*   `--indel`: List indels.
*   `--outgroup`: Indicates presence of outgroups.
*   `-o, --outfile <file>`: Output filename (default: stdout).
