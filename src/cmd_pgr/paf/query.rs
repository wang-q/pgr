use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;
/// Build the clap subcommand for query.
pub fn make_subcommand() -> Command {
    let cmd = Command::new("query")
        .about("Queries PAF index for coordinate projection")
        .after_help(
            r###"
Queries a PAF file or saved index for intervals overlapping a target
region and projects them to query coordinates via CIGAR.

Accepts either a PAF file (built on-the-fly) or a .paf.idx index
(loaded from disk, instant startup).

Two modes:
* Default: single-hop projection — finds all PAF records whose target
  interval overlaps the query region and lifts coordinates to the
  corresponding query sequence.
* --transitive: multi-hop BFS traversal — iteratively projects through
  intermediate sequences up to --max-depth hops.

Region input (one of):
* Positional <region>: single region (e.g. chr1:1000-5000)
* -b/--bed-regions <file>: BED file with multiple regions (one per line,
  tab-separated `name start end`), enabling batch query

Output: PAF (12 columns + gi/bi/cg tags). For BED/MAF/FAS/GFA/VCF
output, use the corresponding `pgr paf to-*` subcommand.

Notes:
* Input PAF files should contain cg:Z: tags for accurate projection
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'

Examples:
1. Single-hop projection from a PAF file:
   pgr paf query alignments.paf chr1:1000-5000

2. Query from a saved index (faster):
   pgr paf query alignments.paf.idx chr1:1000-5000

3. Transitive BFS with filters:
   pgr paf query alignments.paf chr1:1000-5000 --transitive --min-identity 0.8

4. Batch query:
   pgr paf query alignments.paf.idx -b regions.bed

"###,
        );
    crate::cmd_pgr::args::add_query_args(crate::cmd_pgr::args::add_optional_fasta_tsv_arg(cmd))
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the query command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opts = crate::cmd_pgr::args::query_options_from_args(args);
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let (idx, all_results, _fasta_store) = pgr::libs::paf::query::run_query(&opts)?;
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    for (_, results) in &all_results {
        pgr::libs::paf::query::output_paf(&mut writer, &idx, results)?;
    }
    writer.flush()?;
    Ok(())
}
