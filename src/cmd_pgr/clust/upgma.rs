use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::clust::upgma;
use std::io::Write;

/// Build the clap subcommand for upgma.
pub fn make_subcommand() -> Command {
    Command::new("upgma")
        .about("Constructs a phylogenetic tree using UPGMA")
        .after_help(
            r###"
Constructs a phylogenetic tree from a distance matrix using the UPGMA algorithm.

Notes:
* Input: PHYLIP distance matrix (relaxed or strict).
* Output: Newick tree.

Examples:
1. Build tree from matrix:
   pgr clust upgma matrix.phy -o tree.nwk

2. Pipe matrix to tree:
   cat matrix.phy | pgr clust upgma stdin > tree.nwk
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PHYLIP matrix file. [stdin] for standard input",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the upgma command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .ok_or_else(|| anyhow::anyhow!("missing required argument: infile"))?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    // Load matrix
    let matrix = pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(infile)?;

    // Build tree
    let tree = upgma::upgma(&matrix)?;

    // Output tree
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    writer.write_all((tree.to_newick() + "\n").as_ref())?;

    writer.flush()?;
    Ok(())
}
