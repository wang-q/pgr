use clap::*;
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("reroot")
        .about("Reroots a tree at a specified node or the longest branch")
        .after_help(
            r###"
Reroots a phylogenetic tree on a specific branch or node.

Notes:
* Target Nodes (`--node` / `-n`):
    * Can be terminal nodes (leaves) or internal nodes.
    * If multiple nodes are provided, the tree is rerooted on the edge leading to their Lowest Common Ancestor (LCA).
    * This effectively treats the specified nodes as the "ingroup" and the rest as the "outgroup".
    * The new root is placed at the midpoint of the parent edge of the LCA.
    * If the LCA is already the root, the original tree is returned unchanged.
* Default Behavior:
    * If NO nodes are specified, the tree is rerooted at the midpoint of the longest branch.
* Support Values (`--support-as-labels` / `-s`):
    * Standard Newick format stores support values as internal node labels.
    * Rerooting changes edge directions, which can misplace these labels.
    * Use `-s` to shift internal node labels along the path of rerooting to maintain their association with the correct split.
* Topology Changes:
    * The parent edge of the original root (if any) is ignored/merged.
    * Degree-2 nodes (nodes with 1 parent and 1 child) created during the process are automatically removed.

Examples:
1. Reroot at the longest branch (default):
   pgr nwk reroot input.nwk

2. Reroot at a specific node (ingroup):
   pgr nwk reroot input.nwk -n Homo

3. Reroot at the LCA of multiple nodes:
   pgr nwk reroot input.nwk -n Homo -n Pan

4. Reroot and preserve support values (internal node labels):
   pgr nwk reroot input.nwk -n Homo -s

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
            Arg::new("node")
                .long("node")
                .short('n')
                .num_args(1)
                .action(ArgAction::Append)
                .help("Node name"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
        .arg(
            Arg::new("support_as_labels")
                .long("support-as-labels")
                .short('s')
                .action(ArgAction::SetTrue)
                .help("Treat internal node labels as support values and shift them when rerooting"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let process_support = args.get_flag("support_as_labels");

    let infile = args.get_one::<String>("infile").unwrap();
    let mut trees = pgr::libs::phylo::reader::from_file(infile)?;
    
    // Process only the first tree for now, or loop if we want to support multi-tree
    // Since arguments are node names, it implies a single tree context or consistent naming.
    // We'll process the first tree to match previous behavior.
    if let Some(mut tree) = trees.pop() {
        // ids with names
        let id_of: BTreeMap<_, _> = tree.get_name_id();

        // All IDs matched
        let mut ids = BTreeSet::new();
        if let Some(nodes) = args.get_many::<String>("node") {
            for name in nodes {
                if let Some(&id) = id_of.get(name) {
                    ids.insert(id);
                }
            }
        }

        if !ids.is_empty() {
            let mut nodes: Vec<usize> = ids.iter().cloned().collect();
            let mut sub_root_id = nodes.pop().unwrap();

            for id in &nodes {
                sub_root_id = tree.get_common_ancestor(&sub_root_id, id).unwrap();
            }

            let old_root = tree.get_root().unwrap();
            if old_root == sub_root_id {
                let out_string = tree.to_newick();
                writer.write_all((out_string + "\n").as_ref())?;
                return Ok(());
            }

            let new_root = tree.insert_parent(sub_root_id).unwrap();

            // Reroot at the new node
            tree.reroot_at(new_root, process_support).unwrap();

            // Compress: remove degree-2 nodes (redundant internal nodes)
            // The old root likely became a degree-2 node.
            tree.remove_degree_two_nodes();

            let out_string = tree.to_newick();
            writer.write_all((out_string + "\n").as_ref())?;
        } else {
            // Default behavior: Root at the middle of the longest branch
            if let Some(longest_node) = tree.get_node_with_longest_edge() {
                let new_root = tree.insert_parent(longest_node).unwrap();
                tree.reroot_at(new_root, process_support).unwrap();
                tree.remove_degree_two_nodes();
                
                let out_string = tree.to_newick();
                writer.write_all((out_string + "\n").as_ref())?;
            } else {
                // No valid edge to split (e.g. single node tree), just print original
                let out_string = tree.to_newick();
                writer.write_all((out_string + "\n").as_ref())?;
            }
        }
    }

    Ok(())
}

