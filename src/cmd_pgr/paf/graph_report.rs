use clap::*;
use pgr::libs::paf::graph::PafGraph;
use std::collections::HashMap;

pub fn make_subcommand() -> Command {
    Command::new("graph-report")
        .about("Reports coarse GFA topology metrics from PAF alignments")
        .after_help(
            r###"
Computes a topology report (TSV: key<TAB>value) over the coarse pangenome
graph induced from PAF alignments, using the same build path as `pgr paf
graph` (seqwish-style segment-level DSU).

Metrics (~25 dimensions):
* Basic topology: segments, links, paths, path_steps, total_segment_bp
* Segment length distribution: min/mean/median/max
* Node coverage: mean/median, singleton_nodes, reused_nodes,
  reused_nodes_cross_path
* Structure: components, largest_component_nodes, tips, isolated_nodes,
  self_loop_edges
* Path length distribution (steps and bp): min/median/max

Use this to assess graph quality without materializing GFA, or to compare
graphs built with different --min-var-len thresholds.

Notes:
* Input PAF files should contain cg:Z: tags for accurate splitting
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads PAF from stdin if input file is 'stdin'
* FASTA files (-f) are required to populate node sequences and lengths

Examples:
1. Report with default SV threshold (100bp):
   pgr paf graph-report alignments.paf -f refs.fa -o report.tsv

2. Compare thresholds:
   pgr paf graph-report aln.paf -f refs.fa --min-var-len 50  -o r50.tsv
   pgr paf graph-report aln.paf -f refs.fa --min-var-len 500 -o r500.tsv

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input PAF file (or 'stdin' for piped input)"),
        )
        .arg(
            Arg::new("fasta")
                .long("fasta")
                .short('f')
                .required(true)
                .num_args(1..)
                .help("FASTA file(s) providing sequence content and lengths"),
        )
        .arg(
            Arg::new("min_var_len")
                .long("min-var-len")
                .num_args(1)
                .default_value("100")
                .value_parser(clap::value_parser!(i32))
                .help("Minimum indel length to split at (default: 100)"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output TSV report filename. [stdout] for screen"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let infile = matches.get_one::<String>("infile").unwrap();
    let fasta_files: Vec<String> = matches
        .get_many::<String>("fasta")
        .unwrap()
        .cloned()
        .collect();
    let min_var_len = matches
        .get_one::<i32>("min_var_len")
        .copied()
        .unwrap_or(100);
    let outfile = matches.get_one::<String>("outfile").unwrap();

    // Load FASTA sequences into a name -> bytes map.
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    for fa_path in &fasta_files {
        let reader = pgr::reader(fa_path);
        let mut fa_in = noodles_fasta::io::Reader::new(reader);
        for result in fa_in.records() {
            let record = result?;
            let name = String::from_utf8(record.name().into())?;
            let seq_bytes: Vec<u8> = record.sequence()[..].to_vec();
            seqs.insert(name, seq_bytes);
        }
    }

    // Read PAF and build the graph.
    let paf_reader = pgr::reader(infile);
    let graph = PafGraph::build(paf_reader, &seqs, min_var_len)?;

    // Compute and write the report.
    let report = graph.report();
    let writer = pgr::writer(outfile);
    report.write_tsv(writer)?;

    Ok(())
}
