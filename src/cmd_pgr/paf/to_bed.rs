use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for to-bed.
pub fn make_subcommand() -> Command {
    crate::cmd_pgr::args::add_query_args(crate::cmd_pgr::args::add_optional_fasta_tsv_arg(
        Command::new("to-bed"),
    ))
    .about("Queries PAF index and outputs BED3 coordinates")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs query coordinates as BED3 (name start end), one line per result.

This is the pipe-friendly coordinate-only view of `pgr paf query`.
All query options (region, --transitive, filters) are supported.

Notes:
* Input PAF files should contain cg:Z: tags for accurate projection
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'
* --merge-distance requires -f/--fasta-tsv (optional; for CIGAR recomputation)

Examples:
1. Single region to BED:
   pgr paf to-bed alignments.paf chr1:1000-5000

2. Batch query from BED regions:
   pgr paf to-bed alignments.paf.idx -b regions.bed

3. With transitive BFS and identity filter:
   pgr paf to-bed alignments.paf chr1:1000-5000 -t --min-identity 0.8

"###,
    )
    .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the to-bed command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opts = crate::cmd_pgr::args::query_options_from_args(args);
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut inputs: Vec<&str> = vec![opts.infile.as_str()];
    if let Some(tsv) = opts.fasta_tsv.as_deref() {
        inputs.push(tsv);
    }
    if let Some(s) = opts.subset_list.as_deref() {
        inputs.push(s);
    }
    if let Some(s) = opts.syntenic_filter.as_deref() {
        inputs.push(s);
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, inputs)?;
    let (idx, all_results, _fasta_store) = pgr::libs::paf::query::run_query(&opts)?;
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    for (_, results) in &all_results {
        pgr::libs::paf::to_bed::write_bed3(&idx, results, &mut writer)?;
    }
    writer.flush()?;
    Ok(())
}
