use anyhow::{anyhow, Result};
use clap::*;
use pgr::libs::phylo::tree::{distance, Tree};
use std::collections::BTreeMap;
use std::io::{Read, Write};

// Create clap subcommand arguments
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
   pgr nwk distance tree.nwk -m phylip > matrix.phy

4. Distances to parent for leaves only:
   pgr nwk distance tree.nwk -m parent -I
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input filename. [stdin] for standard input"),
        )
        .arg(
            Arg::new("mode")
                .long("mode")
                .short('m')
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("root"),
                    builder::PossibleValue::new("parent"),
                    builder::PossibleValue::new("pairwise"),
                    builder::PossibleValue::new("lca"),
                    builder::PossibleValue::new("phylip"),
                ])
                .default_value("root")
                .help("Set the mode for calculating distances"),
        )
        .arg(
            Arg::new("Internal")
                .long("Internal")
                .short('I')
                .action(ArgAction::SetTrue)
                .help("Ignore internal nodes"),
        )
        .arg(
            Arg::new("Leaf")
                .long("Leaf")
                .short('L')
                .action(ArgAction::SetTrue)
                .help("Ignore leaf nodes"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> Result<()> {
    let mut writer: Box<dyn Write> = pgr::writer(args.get_one::<String>("outfile").unwrap())?;

    let infile = args.get_one::<String>("infile").unwrap();
    let input = match infile.as_str() {
        "stdin" => {
            let mut buffer = String::new();
            std::io::stdin().read_to_string(&mut buffer)?;
            buffer
        }
        _ => std::fs::read_to_string(infile)?,
    };

    // Attempt to parse Newick. If it fails, return error.
    let tree = Tree::from_newick(&input).map_err(|e| anyhow!("Failed to parse Newick: {:?}", e))?;

    let mode = args.get_one::<String>("mode").unwrap();

    let skip_internal = args.get_flag("Internal");
    let skip_leaf = args.get_flag("Leaf");

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
        _ => {}
    }

    Ok(())
}
