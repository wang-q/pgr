use anyhow::{anyhow, Context};
use clap::{ArgMatches, Command};
use pgr::libs::phylo::tree::{distance, Tree};
use std::collections::BTreeMap;
use std::io::{Read, Write};

/// Build the clap subcommand for distance.
pub fn make_subcommand() -> Command {
    Command::new("distance")
        .about("Calculates distances between nodes")
        .after_help(
            r###"
Calculates distances between nodes or generates distance matrices.

Notes:
* Modes:
    * `root`: Distance from each node to the root.
      Output: Node \t Distance
    * `parent`: Distance from each node to its parent.
      Output: Node \t Distance
    * `pairwise`: Distance between every pair of nodes.
      Output: Node1 \t Node2 \t Distance
    * `lca`: Distance from each node in a pair to their Lowest Common Ancestor (LCA).
      Output: Node1 \t Node2 \t Dist1 \t Dist2
    * `phylip`: A Phylip-formatted distance matrix.
      Note: `-I` and `-L` are ignored in this mode.

* The `-I` and `-L` options filter out internal or leaf nodes (except in 'phylip' mode).
* Input must be a valid Newick file.

Examples:
1. Distances to root (default):
   pgr nwk distance tree.nwk

2. Pairwise distances:
   pgr nwk distance tree.nwk --mode pairwise

3. Generate Phylip matrix:
   pgr nwk distance tree.nwk --mode phylip > matrix.phy

4. Distances to parent for leaves only:
   pgr nwk distance tree.nwk --mode parent -I
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::mode_arg(
            "root",
            &["root", "parent", "pairwise", "lca", "phylip"],
            "Set the mode for calculating distances",
        ))
        .arg(crate::cmd_pgr::args::internal_arg())
        .arg(crate::cmd_pgr::args::leaf_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the distance command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer: Box<dyn Write> = Box::new(
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?,
    );

    let infile = args.get_one::<String>("infile").unwrap();
    let mut input = String::new();
    let mut reader =
        pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
    reader
        .read_to_string(&mut input)
        .with_context(|| format!("Failed to read from {}", infile))?;

    // Attempt to parse Newick. If it fails, return error.
    let tree = Tree::from_newick(&input).with_context(|| "Failed to parse Newick")?;

    let mode = args.get_one::<String>("mode").unwrap();

    let skip_internal = args.get_flag("internal");
    let skip_leaf = args.get_flag("leaf");

    // ids with names
    let mut id_of = BTreeMap::new();
    let name_id_map = tree.get_name_id();

    for (name, id) in name_id_map {
        let node = tree
            .get_node(id)
            .ok_or_else(|| anyhow!("node {} not found in tree", id))?;
        let is_leaf = node.children.is_empty();

        if (is_leaf && !skip_leaf) || (!is_leaf && !skip_internal) {
            id_of.insert(name, id);
        }
    }

    match mode.as_str() {
        "root" => distance::dist_root(&tree, &id_of, &mut writer)?,
        "parent" => distance::dist_parent(&tree, &id_of, &mut writer)?,
        "pairwise" => distance::dist_pairwise(&tree, &id_of, &mut writer)?,
        "lca" => distance::dist_lca(&tree, &id_of, &mut writer)?,
        "phylip" => distance::dist_phylip(&tree, &id_of, &mut writer)?,
        _ => anyhow::bail!("unknown distance mode: {}", mode),
    }

    Ok(())
}
