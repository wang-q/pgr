# Change Log

## Unreleased - ReleaseDate

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

#### Clustering and Phylogenetic Analysis

* **`pgr clust`** - Clustering operations (9 subcommands)
  * `hier`: Hierarchical clustering (NN-chain algorithm)
  * `nj`: Neighbor-Joining tree construction
  * `upgma`: UPGMA tree construction
  * `cc`: Connected components clustering
  * `cut`: Tree cutting for cluster extraction
  * `dbscan`: DBSCAN clustering
  * `k-medoids`: K-medoids clustering
  * `mcl`: Markov Clustering (MCL)
  * `eval`: Cluster evaluation metrics

* **`pgr dist`** - Distance and similarity metrics (3 subcommands)
  * `hv`: Hypervector distance calculations
  * `seq`: Sequence distance calculations
  * `vector`: Vector distance calculations

* **`pgr mat`** - Matrix operations (6 subcommands)
  * `compare`: Compare distance matrices
  * `format`: Format matrix files
  * `subset`: Create matrix subset
  * `to-pair`: Convert to pair format
  * `to-phylip`: Convert to PHYLIP format
  * `transform`: Matrix transformations

* **`pgr nwk`** - Newick tree manipulation and visualization (17 subcommands)
  * `stat`: Tree statistics
  * `label`: Label tree nodes
  * `distance`: Calculate pairwise distances
  * `support`: Branch support operations
  * `order`: Order tree nodes
  * `prune`: Prune tree branches
  * `rename`: Rename tree nodes
  * `replace`: Replace node information
  * `reroot`: Reroot the tree
  * `subtree`: Extract subtrees
  * `topo`: Topological operations
  * `comment`: Add comments to trees
  * `indent`: Indent/format tree files
  * `to-dot`: Convert to Graphviz DOT format
  * `to-forest`: Convert to forest representation
  * `to-tex`: Convert to LaTeX/TikZ format
  * `cmp`: Compare trees

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

* **`src/libs/phylo/`** - Phylogenetic analysis core library
  * Tree structure definitions and traversals
  * Tree I/O operations
  * Tree statistics (`stat.rs`)
  * Tree cutting algorithms (`cut.rs`)
  * Tree manipulation algorithms (sorting, rerooting)

* **`src/libs/poa/`** - Partial Order Alignment (POA) implementation

* **`src/libs/chain/`** - Chain/Net alignment processing logic

* **`src/libs/clust/`** - Clustering algorithm implementations
  * `hier.rs`: Hierarchical clustering with NN-chain algorithm
  * `dbscan.rs`: DBSCAN implementation
  * `mcl.rs`: Markov Clustering implementation
  * `k_medoids.rs`: K-medoids implementation

* **`src/libs/io.rs`** - I/O utilities

### Technical Features

* **Rust Implementation**: High-performance, memory-safe implementation
* **Standard Bioinformatics Formats**: Full support for FASTA, FASTQ, 2bit, GFF, AXT, Chain, Net, MAF, PSL, LAV, Newick
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
* **Parsing**: nom 8.0.0, regex 1.11.1
* **Data Structures**: petgraph 0.7.1, indexmap 2.13.0
* **Hashing**: rapidhash, fxhash, murmurhash3, xxhash-rust
* **Output**: rust_xlsxwriter 0.83.0, csv 1.4.0, tera 1.20.1

### Notes

* Some subcommands depend on external executables (UCSC Kent tools, clustalw/muscle/mafft, trf, FastK, etc.)
* Supports both plain text and gzipped (.gz) files for most text-based formats
* 2bit files require random access and do not support stdin or gzipped inputs

## 0.1.0 - 2025-02-08

* New binary `pgr`
