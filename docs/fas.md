# pgr fas

`pgr fas` provides a comprehensive suite of tools for manipulating **Block FA** files. Block FA is a format used to represent multiple sequence alignments (MSA), where each "block" consists of aligned sequences from different species or genomic regions.

## Subcommands

The subcommands are organized into the following categories:

*   **Info**: Extract information or statistics from block FA files.
    *   `check`: Check genome locations against a reference genome FA file.
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
    *   `variation`: List variations (substitutions) in TSV format.

---

## Info Commands

### check

Checks if the genome locations specified in the block headers are valid against a reference genome FA file.

```bash
pgr fas check [OPTIONS] --genome <genome> <infiles>...
```

*   `-g, --genome <path>`: Path to the reference genome FA file (required).
*   `--name <name>`: Check sequences for a specific species only.
*   `-o, --outfile <file>`: Output filename (default: stdout).

Output format (tab-separated): each line contains the entry range followed by its status (`OK` or `FAILED`).

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
*   `--best`: Output nearest-neighbor bilateral links based on sequence distance (deduplicated).
*   `-o, --outfile <file>`: Output filename (default: stdout).

Output format: each line is tab-separated. By default, all ranges in a block are printed on one line. With `--pair` or `--best`, each line contains two ranges.

### name

Extracts all species names found in the block FA files.

```bash
pgr fas name [OPTIONS] <infiles>...
```

*   `-C, --count`: Also output the number of occurrences of each species name.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### stat

Calculates basic statistics for each alignment block (length, comparable bases, differences, gaps, etc.).

```bash
pgr fas stat [OPTIONS] <infiles>...
```

*   `--outgroup`: Treat the last sequence in each block as an outgroup (excludes it from some stats).
*   `-o, --outfile <file>`: Output filename (default: stdout).

Output columns (tab-separated):

*   `target`: Target range of the block.
*   `length`: Alignment length including gaps.
*   `comparable`: Number of positions with unambiguous bases in all sequences.
*   `difference`: Number of polymorphic positions among comparable bases.
*   `gap`: Number of positions where every sequence contains a gap.
*   `ambiguous`: Number of positions with at least one ambiguous base and no gap.
*   `D`: Mean pairwise divergence over all sequence pairs.
*   `indel`: Total span size of all indel regions.

---

## Subset Commands

### filter

Filters blocks based on species presence and sequence length, and optionally formats sequences.

```bash
pgr fas filter [OPTIONS] <infiles>...
```

*   `--name <name>`: Species whose sequence is used for length filtering. Blocks not containing this species are skipped. Defaults to the first species in each block.
*   `--min-len <int>`: Keep blocks where the selected species' alignment length (including gaps) is >= this value.
*   `--max-len <int>`: Keep blocks where the selected species' alignment length (including gaps) is <= this value.
*   `--upper`: Convert sequences to uppercase.
*   `--dash`: Remove dashes (gaps) from sequences.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### slice

Extracts specific slices of alignments based on a provided runlist (JSON).

```bash
pgr fas slice [OPTIONS] --runlist <runlist.json> <infiles>...
```

*   `--runlist <file>`: JSON file describing ranges to extract (required).
*   `--name <name>`: Reference species name (default: first species of the first block).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### subset

Extracts a subset of species from the alignment blocks.

```bash
pgr fas subset [OPTIONS] --required <name.lst> <infiles>...
```

*   `-R, --required <file>`: File with a list of species names to keep, one per line (required).
*   `--strict`: Skip blocks that do not contain *all* the required names.
*   `-o, --outfile <file>`: Output filename (default: stdout).

---

## Transform Commands

### concat

Concatenates sequence pieces of the same species from multiple blocks.

```bash
pgr fas concat [OPTIONS] --required <name.lst> <infiles>...
```

*   `-R, --required <file>`: File with a list of species names to keep/order (required).
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
*   `--align-mode <local|global|semi_global>`: Alignment mode (default: global).
*   `--consensus-name <name>`: Name for the consensus sequence (default: consensus).
*   `--outgroup`: Indicates the last sequence is an outgroup.
*   `-p, --parallel <int>`: Number of threads (default: 1).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### join

Joins multiple block FA files by a common target sequence.

```bash
pgr fas join [OPTIONS] <infiles>...
```

*   `--name <name>`: Target species name. Defaults to the first species of the first block and is used as the common target for all blocks.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### multiz

Merges multiple block FA files in a shared reference coordinate system using a multiz-like banded DP algorithm.

```bash
pgr fas multiz [OPTIONS] --ref-name <name> <infiles>...
```

*   `-r, --ref-name <name>`: Reference sequence name present in all inputs (required).
*   `--radius <int>`: Banded DP radius around reference diagonal (default: 30).
*   `--min-width <int>`: Minimum window width to merge (default: 1).
*   `--mode <core|union>`: Merge mode (default: core).
*   `--score-scheme <file>`: Score scheme file (LASTZ format) or preset name (e.g. hoxd55).
*   `--gap-model <constant|medium|loose>`: Gap model (default: medium).
*   `--align-gap-open <int>`: Alignment gap open cost.
*   `--align-gap-extend <int>`: Alignment gap extension cost.
*   `--match-score <int>`: Match score (default: 2).
*   `--mismatch-score <int>`: Mismatch penalty (default: -1).
*   `--gap-score <int>`: Gap penalty (default: -2).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### refine

Realigns sequences within blocks using built-in or external MSA tools.

```bash
pgr fas refine [OPTIONS] <infiles>...
```

*   `--engine <program>`: Aligning program: `builtin` (default), `clustalw`, `mafft`, `muscle`, `spoa`, `none`.
*   `--outgroup`: Indicates presence of outgroups.
*   `--chop <int>`: Chop head and tail indels (default: 0, disabled).
*   `--quick`: Quick mode, only aligns indel-adjacent regions.
*   `--indel-pad <int>`: In quick mode, enlarge indel regions (default: 50).
*   `--fill <int>`: In quick mode, fill holes between indels (default: 50).
*   `-p, --parallel <int>`: Number of threads (default: 1).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### replace

Replaces headers in block FA files using a mapping file.

```bash
pgr fas replace [OPTIONS] --replace-tsv <replace.tsv> <infiles>...
```

*   `--replace-tsv <file>`: TSV file containing replacement rules (required). Each line is a tab-separated list:
    *   One field: if the name uniquely matches one header in a block, the whole block is dropped.
    *   Two fields: `original_name<TAB>new_name` replaces the matching header.
    *   Three or more fields: duplicates the block once for every replacement name after the first.
    *   If a block contains multiple matching headers, the block is kept unchanged and a warning is emitted.
*   `-o, --outfile <file>`: Output filename (default: stdout).

A header that appears more than once within the same block is also treated as multiple matching headers, and the block will be kept unchanged.

---

## File Commands

### create

Creates block FA files from links of ranges (e.g., from `pgr fas link`).

```bash
pgr fas create [OPTIONS] --genome <genome> <infiles>...
```

*   `-g, --genome <file>`: Path to the reference genome FA file (required).
*   `--name <name>`: Set a species name for ranges (default: inferred from header).
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
*   `--no-single`: Omit singleton variations.
*   `--no-complex`: Omit complex variations.
*   `--min-freq <float>`: Minimal frequency.
*   `--max-freq <float>`: Maximal frequency.
*   `--wrap <int>`: Wrap length for visualization (default: 50).
*   `-o, --outfile <file>`: Output filename (default: variations.xlsx).

### variation

Lists variations (substitutions) in TSV format.

```bash
pgr fas variation [OPTIONS] <infiles>...
```

*   `--outgroup`: Indicates presence of outgroups.
*   `-o, --outfile <file>`: Output filename (default: stdout).

Output columns (tab-separated):

*   `#target`: Target range of the block.
*   `chr`: Chromosome name.
*   `chr_pos`: Position on the chromosome.
*   `range`: Chromosome position formatted as `chr:pos`.
*   `pos`: Position within the alignment (1-based).
*   `tbase`, `qbase`, `bases`, `mutant_to`, `freq`, `pattern`, `obase`: Fields from the substitution record (see `pgr::libs::alignment::Substitution`).
