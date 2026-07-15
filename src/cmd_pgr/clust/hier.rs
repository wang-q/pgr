use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;
use std::str::FromStr;

use pgr::libs::clust::hier::{linkage_inplace, to_tree, Method};
use pgr::libs::pairmat::NamedMatrix;
use pgr::libs::phylo::tree::io::to_newick;
/// Build the clap subcommand for hier.
pub fn make_subcommand() -> Command {
    Command::new("hier")
        .about("Clusters entries via hierarchical clustering")
        .visible_alias("hclust")
        .after_help(
            r###"
This command performs hierarchical clustering (agglomerative) on a distance matrix.
It takes a PHYLIP format distance matrix as input and produces a Newick tree string.

Notes:
* Input matrix must be in PHYLIP format (strict or relaxed).
* If you have a pairwise list (name1 name2 dist), use `pgr mat to-phylip` first.
* The output Newick tree uses the linkage distance (merge height) as node height.
* For Ward's method, the input is assumed to be Euclidean distances (or similar).

Examples:
1. Basic usage with Ward's method (default):
   pgr clust hier matrix.phy > tree.nwk

2. Use Average Linkage (UPGMA):
   pgr clust hier matrix.phy --method average > tree.nwk

3. Use Single Linkage (Nearest Neighbor):
   pgr clust hier matrix.phy --method single > tree.nwk
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PHYLIP distance matrix file. [stdin] for standard input",
        ))
        .arg(crate::cmd_pgr::args::clust_method_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the hier command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .ok_or_else(|| anyhow::anyhow!("missing required argument: infile"))?;
    // clust_method has a clap default value, so unwrap is safe.
    let method_str = args.get_one::<String>("clust_method").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    // Parse method
    let method = Method::from_str(method_str)
        .map_err(|e| anyhow::anyhow!("invalid --clust-method '{}': {}", method_str, e))?;

    // Read matrix
    let matrix = NamedMatrix::from_relaxed_phylip(infile)?;

    // Perform clustering
    let (names, condensed) = matrix.into_parts();
    let steps = linkage_inplace(condensed, method);

    // Convert to tree
    let tree = to_tree(&steps, &names)?;

    // Format output
    let newick = to_newick(&tree);

    // Write output
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    writer.write_all((newick + "\n").as_ref())?;

    writer.flush()?;
    Ok(())
}
