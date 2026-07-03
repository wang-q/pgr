use clap::{ArgMatches, Command};
use pgr::libs::clust::nj;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("nj")
        .about("Constructs a phylogenetic tree using Neighbor-Joining")
        .after_help(
            r###"
Constructs a phylogenetic tree from a distance matrix using the Neighbor-Joining (NJ) algorithm.

Notes:
* Input: PHYLIP distance matrix (relaxed or strict).
* Output: Newick tree (midpoint rooted).
* NJ is a bottom-up clustering method suitable for variable evolutionary rates.

Examples:
1. Build tree from matrix:
   pgr clust nj matrix.phy -o tree.nwk

2. Pipe matrix to tree:
   cat matrix.phy | pgr clust nj stdin > tree.nwk
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PHYLIP matrix file. [stdin] for standard input",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    // Load matrix
    let matrix = pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(infile)?;

    // Build tree
    let tree = nj::nj(&matrix)?;

    // Output tree
    let mut writer = pgr::writer(outfile)?;
    writer.write_all((tree.to_newick() + "\n").as_ref())?;

    Ok(())
}
