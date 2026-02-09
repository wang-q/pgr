use clap::*;
use pgr::libs::phylo::node::NodeId;
use pgr::libs::phylo::tree::Tree;
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
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());

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
    let tree = Tree::from_newick(&input)
        .map_err(|e| anyhow::anyhow!("Failed to parse Newick: {:?}", e))?;

    let mode = args.get_one::<String>("mode").unwrap();

    let skip_internal = args.get_flag("Internal");
    let skip_leaf = args.get_flag("Leaf");

    // ids with names
    let mut id_of = BTreeMap::new();
    let name_id_map = tree.get_name_id();

    for (name, id) in name_id_map {
        let node = tree.get_node(id).unwrap();
        let is_leaf = node.children.is_empty();

        if is_leaf && !skip_leaf {
            id_of.insert(name, id);
        } else if !is_leaf && !skip_internal {
            id_of.insert(name, id);
        }
    }

    match mode.as_str() {
        "root" => dist_root(&tree, &id_of, &mut writer),
        "parent" => dist_parent(&tree, &id_of, &mut writer),
        "pairwise" => dist_pairwise(&tree, &id_of, &mut writer),
        "lca" => dist_lca(&tree, &id_of, &mut writer),
        "phylip" => dist_phylip(&tree, &id_of, &mut writer),
        _ => unreachable!(),
    }

    Ok(())
}

fn dist_root(tree: &Tree, id_of: &BTreeMap<String, NodeId>, writer: &mut Box<dyn Write>) {
    let root = tree.get_root().unwrap();
    for (k, v) in id_of.iter() {
        let dist = {
            let (edge_sum, num_edges) = tree.get_distance(&root, v).unwrap();
            if edge_sum.abs() > 1e-9 {
                edge_sum
            } else {
                num_edges as f64
            }
        };
        writer
            .write_fmt(format_args!("{}\t{}\n", k, format_float(dist)))
            .unwrap();
    }
}

fn dist_parent(tree: &Tree, id_of: &BTreeMap<String, NodeId>, writer: &mut Box<dyn Write>) {
    for (k, v) in id_of.iter() {
        let parent = tree.get_node(*v).unwrap().parent;
        if parent.is_none() {
            writer.write_fmt(format_args!("{}\t0\n", k)).unwrap();
            continue;
        }
        let parent = parent.unwrap();

        let dist = {
            let (edge_sum, num_edges) = tree.get_distance(&parent, v).unwrap();
            if edge_sum.abs() > 1e-9 {
                edge_sum
            } else {
                num_edges as f64
            }
        };
        writer
            .write_fmt(format_args!("{}\t{}\n", k, format_float(dist)))
            .unwrap();
    }
}

fn dist_pairwise(tree: &Tree, id_of: &BTreeMap<String, NodeId>, writer: &mut Box<dyn Write>) {
    for (k1, v1) in id_of.iter() {
        for (k2, v2) in id_of.iter() {
            let dist = {
                let (edge_sum, num_edges) = tree.get_distance(v1, v2).unwrap();
                if edge_sum.abs() > 1e-9 {
                    edge_sum
                } else {
                    num_edges as f64
                }
            };
            writer
                .write_fmt(format_args!("{}\t{}\t{}\n", k1, k2, format_float(dist)))
                .unwrap();
        }
    }
}

fn dist_lca(tree: &Tree, id_of: &BTreeMap<String, NodeId>, writer: &mut Box<dyn Write>) {
    for (k1, v1) in id_of.iter() {
        for (k2, v2) in id_of.iter() {
            let lca = tree.get_common_ancestor(v1, v2).unwrap();

            let dist1 = {
                let (edge_sum, num_edges) = tree.get_distance(&lca, v1).unwrap();
                if edge_sum.abs() > 1e-9 {
                    edge_sum
                } else {
                    num_edges as f64
                }
            };

            let dist2 = {
                let (edge_sum, num_edges) = tree.get_distance(&lca, v2).unwrap();
                if edge_sum.abs() > 1e-9 {
                    edge_sum
                } else {
                    num_edges as f64
                }
            };
            writer
                .write_fmt(format_args!(
                    "{}\t{}\t{}\t{}\n",
                    k1,
                    k2,
                    format_float(dist1),
                    format_float(dist2)
                ))
                .unwrap();
        }
    }
}

fn dist_phylip(tree: &Tree, id_of: &BTreeMap<String, NodeId>, writer: &mut Box<dyn Write>) {
    let names: Vec<&String> = id_of.keys().collect();
    let n = names.len();

    // Phylip header
    writer.write_fmt(format_args!("    {}\n", n)).unwrap();

    for (i, name) in names.iter().enumerate() {
        let v1 = id_of.get(*name).unwrap();

        // Name padding to 10 chars usually, but let's just print name followed by tab/space
        // Phylip strict format requires 10 chars for name.
        // Relaxed format (which is common) allows longer names separated by whitespace.
        // Let's print name then spaces.
        writer.write_fmt(format_args!("{} ", name)).unwrap();

        for (j, other_name) in names.iter().enumerate() {
            let v2 = id_of.get(*other_name).unwrap();
            let dist = if i == j {
                0.0
            } else {
                let (edge_sum, num_edges) = tree.get_distance(v1, v2).unwrap();
                if edge_sum.abs() > 1e-9 {
                    edge_sum
                } else {
                    num_edges as f64
                }
            };

            writer.write_fmt(format_args!(" {:.6}", dist)).unwrap();
        }
        writer.write_all(b"\n").unwrap();
    }
}

fn format_float(val: f64) -> String {
    let rounded = (val * 1e6).round() / 1e6;
    format!("{}", rounded)
}
