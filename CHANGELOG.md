# Change Log

## Unreleased - ReleaseDate

### New Features

#### Pangenome (PAF)

* **`pgr paf`** - PAF (Pairwise mApping Format) file manipulation (9 subcommands).
  Treats pairwise alignments as an implicit pangenome graph; query commands
  traverse the graph on demand without materializing it.
  * `index`: Build a persistent `.paf.idx` index (interval tree + BFS transitive
    closure) over one or more PAF files, with BGZF/gzip support and lazy CIGAR
    loading.
  * `query`: Project target intervals through the alignment network; supports
    `--transitive` (A→B→C chaining), `-b/--batch` batch queries, and `-o bed`
    BED output.
  * `to-bed`: Convert PAF records (or query results) to BED intervals.
  * `to-fas`: Extract pairwise alignments as block FA via on-demand FASTA slicing.
  * `to-maf`: Reconstruct MAF from PAF CIGAR; supports multi-way POA MSA and
    minus-strand output.
  * `to-vcf`: Emit VCF from POA multiple sequence alignment; supports
    substitutions and indels with left-alignment, plus rGFA tags.
  * `to-gfa`: Produce local GFA with unchopping, GFA header, and crush mode.
  * `graph`: Materialize a coarse whole-genome graph (splits nodes only at
    variants ≥ `--min-var-len`, default 100 bp).
  * `stat`: Report alignment statistics; supports topology-only mode.

* **`pgr maf to-paf`** - Convert MAF alignments to PAF format (new `maf`
  subcommand; `maf` now has 2 subcommands).

#### Population Genome Compression (pbit)

* **`pgr pbit`** - Native population genome compression format (6 subcommands).
  Stores a reference genome as standard 2bit and each sample as LZ-diff or
  PAF-driven CIGAR delta encoding, preserving random access while compressing
  large cohorts of homologous samples.
  * `create`: Build a new `.pbit` archive from a reference FASTA and sample
    FASTAs (with optional PAF for delta encoding).
  * `append`: Add samples to an existing `.pbit` archive.
  * `range`: Extract coordinate intervals from all samples.
  * `some`: Extract complete contigs by name from all samples.
  * `stat`: Display archive statistics, sample list, and contig list.
  * `to-fa`: Export sample sequences to per-sample FASTA files.

### Enhancements

* **`pgr fas`**: Added outgroup support for `stat`; improved `concat` order
  handling and duplicate header detection; renamed `--aligner` to `--engine`;
  consolidated shared logic into `libs/fas_multiz`.
* **CLI standardization**: Centralized repeated `clap` argument definitions into
  a shared `cmd_pgr::args` module; unified flag names and short aliases across
  commands.
* **Error handling**: Replaced `unwrap`/`expect`/`unreachable!` with
  `anyhow::Result` propagation and contextual errors across the codebase.
* **Code quality**: Extracted shared logic into `libs/` modules (POA args, FASTA
  utilities, range reversal, etc.); modernized imports and deduplicated code.
* **Logging**: Added `log`/`env_logger` based logging in pipeline commands.
* **Tests**: Switched all CLI tests to temporary directories; added
  comprehensive fixtures for `paf`, `pbit`, `2bit`, and `dist`; added
  regression tests for minus-strand CIGAR and edge cases.
* **Docs**: Migrated scattered docs to dedicated per-command files under
  `docs/`; added `docs/formats/` for format references; standardized
  terminology to "Block FA".

### Core Libraries

* **`src/libs/paf/`** - PAF implicit graph core: record I/O, interval tree
  index (`coitrees`), bidirectional CIGAR handling, BFS transitive closure,
  and graph materialization.
* **`src/libs/pbit/`** - pbit archive format: LZ-diff and PAF-driven CIGAR
  delta encoding, with random-access decompression.
* **`src/libs/poa/`** - Partial Order Alignment engine (Spoa C++ → Rust port),
  used by `paf to-maf`/`to-vcf` and `fas refine`.
* **CIGAR support** in `libs/paf/cigar.rs`: splits `M` into `=`/`X`,
  case-insensitive matching, and bounds-checked operations.
* **`src/libs/loc.rs`** - FASTA random-access index (`.loc`) with BGZF support,
  shared by `fa range`, `paf`, and `pbit`.

### Changed

* Migrated `clust`, `cut`, `eval`, `mat`, `nwk` command modules and associated
  libraries (`libs/phylo/`, `libs/clust/`, `libs/cut/`, `libs/eval/`,
  `libs/pairmat/`) to the `necom` project. `pgr` now focuses on genome data
  processing (sequences, alignments, pangenome, pipelines, plotting).
* `pgr dist` reduced from 3 to 2 subcommands (`vector` migrated to `necom`).
* `pgr pl` reduced from 7 to 6 subcommands (`condense` migrated to `necom`).
* `pgr maf` expanded from 1 to 2 subcommands (added `to-paf`).

### Removed

* `nom` dependency (Newick parsing moved to `necom`).
* `xxhash-rust` dependency (unused; never appeared in `Cargo.toml`).

### Dependencies

* **Added**: `coitrees` (interval tree for PAF/loc), `intspan` (interval set),
  `bincode` + `serde` + `serde_json` (PAF index persistence), `lru` (caching),
  `minimizer-iter` (pbit/paf), `cmd_lib` (pipelines), `which` (executable
  lookup), `flate2` (gzip), `log` + `env_logger` (logging), `itertools`,
  `rand`.
* **Removed**: `nom`.

## 0.2.0 - 2026-04-05

### New Features

#### Sequence Manipulation

* **`pgr fa`** - FASTA file manipulation toolkit (18 subcommands)
  * `size`: Display sequence lengths
  * `count`: Count sequences and total bases
  * `masked`: Report masked (lowercase/N) regions
  * `n50`: Calculate N50 statistics
  * `one`: Extract a single sequence by name or index
  * `some`: Extract multiple sequences by name list
  * `order`: Reorder sequences
  * `split`: Split sequences into separate files
  * `window`: Generate sliding windows
  * `replace`: Replace sequence contents
  * `rc`: Reverse complement sequences
  * `filter`: Filter sequences by length or pattern
  * `dedup`: Remove duplicate sequences
  * `mask`: Mask regions of sequences
  * `six-frame`: Generate six-frame translations
  * `to-2bit`: Convert FASTA to 2bit format
  * `gz`: Compress/decompress FASTA files
  * `range`: Extract specific ranges from indexed FASTA

* **`pgr fq`** - FASTQ file manipulation (2 subcommands)
  * `to-fa`: Convert FASTQ to FASTA
  * `interleave`/`il`: Interleave paired-end FASTQ files

* **`pgr fas`** - Block FA (multiple alignment) manipulation (19 subcommands)
  * `check`: Validate block FA format
  * `cover`: Calculate coverage statistics
  * `link`: Generate link information
  * `name`: Manipulate sequence names
  * `stat`: Display alignment statistics
  * `filter`: Filter alignments
  * `slice`: Extract slices from alignments
  * `subset`: Create subset of alignments
  * `concat`: Concatenate alignments
  * `consensus`: Generate consensus sequences
  * `join`: Join alignments
  * `refine`: Refine alignments using external tools
  * `replace`: Replace sequences in alignments
  * `create`: Create block FA from other formats
  * `separate`: Separate alignments into files
  * `split`: Split alignments
  * `to-vcf`: Convert to VCF format
  * `to-xlsx`: Export to Excel format
  * `variation`: Analyze variations
  * `multiz`: Handle Multiz alignment format

* **`pgr 2bit`** - 2bit file manipulation (5 subcommands)
  * `masked`: Report masked regions
  * `size`: Display sequence sizes
  * `range`: Extract sequence ranges
  * `some`: Extract specific sequences
  * `to-fa`: Convert to FASTA format

* **`pgr gff`** - GFF file operations (1 subcommand)
  * `rg`: GFF range query operations

#### Alignment Formats

* **`pgr axt`** - AXT pairwise alignment manipulation (4 subcommands)
  * `sort`: Sort AXT alignments
  * `to-fas`: Convert to block FA format
  * `to-maf`: Convert to MAF format
  * `to-psl`: Convert to PSL format

* **`pgr chain`** - Chain alignment manipulation (6 subcommands)
  * `anti-repeat`: Remove repetitive alignments
  * `net`: Generate net from chain
  * `pre-net`: Pre-net processing
  * `sort`: Sort chain files
  * `split`: Split chain files
  * `stitch`: Stitch chain alignments

* **`pgr net`** - Net alignment manipulation (6 subcommands)
  * `class`: Classify net alignments
  * `filter`: Filter net alignments
  * `split`: Split net files
  * `subset`: Create subset of net
  * `syntenic`: Extract syntenic regions
  * `to-axt`: Convert to AXT format

* **`pgr maf`** - MAF multiple alignment manipulation (1 subcommand)
  * `to-fas`: Convert to block FA format

* **`pgr psl`** - PSL alignment manipulation (8 subcommands)
  * `chain`: Convert to chain format
  * `histo`: Generate histogram statistics
  * `lift`: Lift coordinates
  * `rc`: Reverse complement
  * `stats`: Display alignment statistics
  * `swap`: Swap target/query
  * `to-chain`: Convert to chain format
  * `to-range`: Convert to range format

* **`pgr lav`** - LAV alignment manipulation (2 subcommands)
  * `lastz`: LASTZ-specific operations
  * `to-psl`: Convert to PSL format

#### Distance

* **`pgr dist`** - Distance and similarity metrics (2 subcommands)
  * `hv`: Hypervector distance calculations
  * `seq`: Sequence distance calculations

#### Other Tools

* **`pgr ms`** - Hudson's ms simulator tools (1 subcommand)
  * `to-dna`: Convert ms output to DNA sequences

* **`pgr pl`** - Integrated pipelines (6 subcommands)
  * `p2m`: Pair-to-multi alignment pipeline
  * `prefilter`: Prefiltering pipeline
  * `trf`: Tandem repeat finder pipeline
  * `ir`: Inverted repeat detection pipeline
  * `rept`: Repeat analysis pipeline
  * `ucsc`: UCSC Kent tools wrapper pipeline

* **`pgr plot`** - Plotting tools (3 subcommands)
  * `hh`: Histogram plotting
  * `nrps`: NRPS visualization
  * `venn`: Venn diagram generation

### Core Libraries

* **`src/libs/poa/`** - Partial Order Alignment (POA) implementation

* **`src/libs/chain/`** - Chain/Net alignment processing logic

* **`src/libs/io.rs`** - I/O utilities

### Technical Features

* **Rust Implementation**: High-performance, memory-safe implementation
* **Standard Bioinformatics Formats**: Full support for FASTA, FASTQ, 2bit, GFF, AXT, Chain, Net, MAF, PSL, LAV
* **Parallel Computing**: Rayon-based parallelism for performance-critical operations
* **Zero Panic Policy**: Robust error handling for malformed inputs
* **Pipeline-Friendly**: stdin/stdout support where possible, predictable outputs
* **Comprehensive CLI**: Built with clap for excellent command-line experience
* **Testing**: Extensive integration tests using assert_cmd

### Dependencies

* **CLI**: clap 4.5.28
* **Error Handling**: anyhow 1.0.93
* **Bioinformatics**: noodles 0.104.0, bio 0.30.1
* **Parallelism**: rayon 1.10.0, crossbeam 0.8.4
* **Parsing**: regex 1.11.1
* **Data Structures**: petgraph 0.7.1, indexmap 2.13.0
* **Hashing**: rapidhash, fxhash, murmurhash3
* **Output**: rust_xlsxwriter 0.83.0, csv 1.4.0, tera 1.20.1

### Notes

* Some subcommands depend on external executables (UCSC Kent tools, clustalw/muscle/mafft, trf, FastK, etc.)
* Supports both plain text and gzipped (.gz) files for most text-based formats
* 2bit files require random access and do not support stdin or gzipped inputs

## 0.1.0 - 2025-02-08

* New binary `pgr`
