use clap::*;
use pgr::libs::paf::fasta::{load_fasta_tsv, FastaStore};
use pgr::libs::paf::graph::PafGraph;
use std::collections::HashMap;

pub fn make_subcommand() -> Command {
    Command::new("graph")
        .about("Induces a coarse GFA graph from PAF alignments")
        .after_help(
            r###"
Builds a coarse pangenome graph (GFA v1.0) from pairwise PAF alignments
and a set of FASTA sequences, using a seqwish-style segment-level DSU
algorithm.

Algorithm:
* Splits each alignment at indels >= --min-var-len into match segments.
* Unions aligned segments via disjoint-set union (transitive closure).
* Derives graph nodes (DSU classes) + edges (path adjacencies) + novel
  segments (unaligned gaps), then emits S/L/P lines.

This builds a coarse global GFA: only large structural variations
(>= --min-var-len) split nodes; small indels stay within a node as
variations. For base-level regional graphs, see `pgr paf to-gfa`.

-f/--fasta-tsv (optional):
  TSV with two columns: genome_name <tab> bgzf_fasta_path
  Each genome_name must match a query/target name in the PAF index.
  Omit for topology-only mode (S lines emit '*' with LN:i: tags).

Notes:
* Input PAF files should contain cg:Z: tags for accurate splitting
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads PAF from stdin if input file is 'stdin'
* GFA node ids are 1-based; node 1 is the earliest segment by (seq, start)
* S lines carry rGFA tags: SN:Z (source seq), SO:i (0-based start), SR:i:0

Examples:
1. Build a coarse graph with default SV threshold (100bp):
   pgr paf graph alignments.paf -f genomes.tsv -o graph.gfa

2. Stricter threshold (only >= 500bp SVs split nodes):
   pgr paf graph alignments.paf -f genomes.tsv --min-var-len 500 -o graph.gfa

3. Read PAF from stdin:
   cat alignments.paf | pgr paf graph stdin -f genomes.tsv > graph.gfa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input PAF file (or 'stdin' for piped input)"),
        )
        .arg(
            Arg::new("fasta_tsv")
                .long("fasta-tsv")
                .short('f')
                .num_args(1)
                .help(
                    "TSV file: genome_name <tab> bgzf_fasta_path (optional for topology-only mode)",
                ),
        )
        .arg(
            Arg::new("min_var_len")
                .long("min-var-len")
                .num_args(1)
                .default_value("100")
                .value_parser(clap::value_parser!(i32))
                .help("Minimum indel length to split at (default: 100)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let infile = matches.get_one::<String>("infile").unwrap();
    let tsv_path = matches.get_one::<String>("fasta_tsv");
    let min_var_len = matches
        .get_one::<i32>("min_var_len")
        .copied()
        .unwrap_or(100);
    let outfile = crate::cmd_pgr::args::get_outfile(matches);

    // Load FASTA sequences via TSV + FastaStore (optional for topology-only mode).
    let seqs: HashMap<String, Vec<u8>> = if let Some(tsv) = tsv_path {
        let seq_to_file = load_fasta_tsv(tsv)?;
        let mut store = FastaStore::new(&seq_to_file)?;
        let mut map = HashMap::new();
        for name in seq_to_file.keys() {
            map.insert(name.clone(), store.fetch_full(name)?);
        }
        map
    } else {
        HashMap::new()
    };

    // Read PAF.
    let paf_reader = pgr::reader(infile)?;

    // Build the graph.
    let seqs_ref = if seqs.is_empty() { None } else { Some(&seqs) };
    let graph = PafGraph::build(paf_reader, seqs_ref, min_var_len)?;

    // Write GFA.
    let writer = pgr::writer(outfile)?;
    graph.write_gfa(writer)?;

    Ok(())
}
