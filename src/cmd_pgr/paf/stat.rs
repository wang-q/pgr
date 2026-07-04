use clap::{ArgMatches, Command};
use pgr::libs::paf::fasta::load_all_seqs;
use pgr::libs::paf::graph::PafGraph;
/// Build the clap subcommand for stat.
pub fn make_subcommand() -> Command {
    let cmd = Command::new("stat")
        .about("Reports coarse pangenome graph topology metrics from PAF alignments")
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

-f/--fasta-tsv (optional):
  TSV with two columns: genome_name <tab> bgzf_fasta_path
  Each genome_name must match a query/target name in the PAF index.
  Omit for topology-only mode (node lengths inferred from segment coords).

Notes:
* Input PAF files should contain cg:Z: tags for accurate splitting
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads PAF from stdin if input file is 'stdin'

Examples:
1. Report with default SV threshold (100bp):
   pgr paf stat alignments.paf -f genomes.tsv -o report.tsv

2. Compare thresholds:
   pgr paf stat aln.paf -f genomes.tsv --min-var-len 50  -o r50.tsv
   pgr paf stat aln.paf -f genomes.tsv --min-var-len 500 -o r500.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PAF file (or 'stdin' for piped input)",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg());
    let cmd = crate::cmd_pgr::args::add_optional_fasta_tsv_arg(cmd);
    crate::cmd_pgr::args::add_min_var_len_arg(cmd)
}
/// Execute the stat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let tsv_path = args.get_one::<String>("fasta_tsv").map(|s| s.as_str());
    let min_var_len = args.get_one::<i32>("min_var_len").copied().unwrap_or(100);
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    // Load FASTA sequences via TSV + FastaStore (optional for topology-only mode).
    let seqs = load_all_seqs(tsv_path)?;

    // Read PAF and build the graph.
    let paf_reader = pgr::reader(infile)?;
    let seqs_ref = if seqs.is_empty() { None } else { Some(&seqs) };
    let graph = PafGraph::build(paf_reader, seqs_ref, min_var_len)?;

    // Compute and write the report.
    let report = graph.report();
    let writer = pgr::writer(outfile)?;
    report.write_tsv(writer)?;

    Ok(())
}
