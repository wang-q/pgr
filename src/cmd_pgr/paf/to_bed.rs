use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for to-bed.
pub fn make_subcommand() -> Command {
    crate::cmd_pgr::args::add_query_args(Command::new("to-bed"))
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

Examples:
1. Single region to BED:
   pgr paf to-bed alignments.paf chr1:1000-5000

2. Batch query from BED regions:
   pgr paf to-bed alignments.paf.idx -b regions.bed

3. With transitive BFS and identity filter:
   pgr paf to-bed alignments.paf chr1:1000-5000 -t --min-identity 0.8

"###,
        )
}
/// Execute the to-bed command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opts = crate::cmd_pgr::args::query_options_from_args(args);
    let (idx, all_results) = pgr::libs::paf::query::run_query(&opts)?;
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    for (_, results) in &all_results {
        pgr::libs::paf::to_bed::write_bed3(&idx, results, &mut out)?;
    }
    out.flush()?;
    Ok(())
}
