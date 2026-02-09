use clap::*;
use pgr::libs::phylo::build;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("nj")
        .about("Construct a phylogenetic tree using Neighbor-Joining")
        .after_help(
            r###"
Constructs a phylogenetic tree from a distance matrix using the Neighbor-Joining (NJ) algorithm.

Notes:
* Input: PHYLIP distance matrix (relaxed or strict).
* Output: Newick tree (midpoint rooted).
* NJ is a bottom-up clustering method suitable for variable evolutionary rates.

Examples:
1. Build tree from matrix:
   pgr mat nj matrix.phy -o tree.nwk

2. Pipe matrix to tree:
   cat matrix.phy | pgr mat nj - > tree.nwk
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input PHYLIP matrix file. [stdin] for standard input"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();

    // Load matrix
    let matrix = intspan::NamedMatrix::from_relaxed_phylip(infile);

    // Build tree
    let tree = build::nj(&matrix)?;

    // Output tree
    let mut writer = intspan::writer(outfile);
    writer.write_all((tree.to_newick() + "\n").as_ref())?;

    Ok(())
}
