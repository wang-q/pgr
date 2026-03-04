use clap::*;
use std::io::Write;
use std::str::FromStr;

use pgr::libs::clust::hier::{Method, linkage_inplace, to_tree};
use pgr::libs::pairmat::NamedMatrix;
use pgr::libs::phylo::tree::io::to_newick;

pub fn make_subcommand() -> Command {
    Command::new("hier")
        .about("Hierarchical clustering (dendrogram)")
        .visible_alias("hclust")
        .after_help(
            r#"
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
"#,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input PHYLIP distance matrix file. [stdin] for standard input"),
        )
        .arg(
            Arg::new("method")
                .long("method")
                .short('m')
                .default_value("ward")
                .help("Clustering method (single, complete, average, weighted, centroid, median, ward)"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .default_value("stdout")
                .help("Output file path. [stdout] for screen"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let infile = matches.get_one::<String>("infile").unwrap();
    let method_str = matches.get_one::<String>("method").unwrap();
    let outfile = matches.get_one::<String>("outfile").unwrap();

    // Parse method
    let method = Method::from_str(method_str)
        .map_err(|e: String| anyhow::anyhow!(e))?;

    // Read matrix
    let matrix = NamedMatrix::from_relaxed_phylip(infile);

    // Perform clustering
    let (names, condensed) = matrix.into_parts();
    let steps = linkage_inplace(condensed, method);

    // Convert to tree
    let tree = to_tree(&steps, &names);

    // Format output
    let newick = to_newick(&tree);

    // Write output
    let mut writer = pgr::writer(outfile);
    writer.write_all((newick + "\n").as_ref())?;

    Ok(())
}
