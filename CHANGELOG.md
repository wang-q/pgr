# Change Log

## Unreleased - ReleaseDate

### Changed

* Migrated `clust`, `cut`, `eval`, `mat`, `nwk` command modules and associated
  libraries (`libs/phylo/`, `libs/clust/`, `libs/cut/`, `libs/eval/`,
  `libs/pairmat/`) to the `necom` project. `pgr` now focuses on genome data
  processing (sequences, alignments, pangenome, pipelines, plotting).
* Removed `nom` dependency (Newick parsing moved to `necom`).

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
