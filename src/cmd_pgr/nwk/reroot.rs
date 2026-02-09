use clap::*;
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use pgr::libs::phylo::tree::Tree;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("reroot")
        .about("Reroots a tree at a specified node or the longest branch")
        .after_help(
            r###"
Reroots a phylogenetic tree on a specific branch or node.

Notes:
* Target Selection:
    * Default: If NO nodes are specified, reroots at the midpoint of the longest branch.
    * Specified Nodes (`--node` / `-n`):
        * Reroots on the edge leading to the Lowest Common Ancestor (LCA) of specified nodes.
        * Treats specified nodes as the "ingroup" and the rest as the "outgroup".
        * The new root is placed at the midpoint of the parent edge of the LCA.
    * Lax Mode (`--lax` / `-l`):
        * If the LCA of specified nodes is already the root, use the *unspecified* nodes (complement) as the ingroup.
        * Useful for defining an outgroup by exclusion.

* Operations:
    * Reroot (Default): Creates a bifurcating root at the target edge.
    * Deroot (`--deroot` / `-d`): Splices out the ingroup to create a multifurcating root.

* Technical Details:
    * Support Values (`--support-as-labels` / `-s`):
        * Shifts internal node labels along the rerooting path to maintain split associations.
        * Necessary because rerooting changes edge directions.
    * Topology:
        * The original root's parent edge is merged.
        * Degree-2 nodes created during the process are removed.

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
        .arg(
            Arg::new("deroot")
                .long("deroot")
                .short('d')
                .action(ArgAction::SetTrue)
                .help("Deroot the tree (create a multifurcating root) (see Notes)"),
        )
        .arg(
            Arg::new("lax")
                .long("lax")
                .short('l')
                .action(ArgAction::SetTrue)
                .help("Lax mode: Use the complement if the specified nodes form the root (see Notes)"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let process_support = args.get_flag("support_as_labels");
    let deroot = args.get_flag("deroot");
    let lax = args.get_flag("lax");

    let infile = args.get_one::<String>("infile").unwrap();
    let mut trees = Tree::from_file(infile)?;

    // Process only the first tree for now, or loop if we want to support multi-tree
    // Since arguments are node names, it implies a single tree context or consistent naming.
    // We'll process the first tree to match previous behavior.
    if let Some(mut tree) = trees.pop() {
        if deroot {
            tree.deroot().map_err(|e| anyhow::anyhow!(e))?;
        } else {
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

                // Lax mode check
                if old_root == sub_root_id && lax {
                    // Identify leaves in the "ingroup" (specified nodes)
                    let mut specified_leaves = BTreeSet::new();
                    for &id in &ids {
                        // Find all descendant leaves of this node
                        if let Ok(subtree) = tree.get_subtree(&id) {
                            for sub_id in subtree {
                                if let Some(node) = tree.get_node(sub_id) {
                                    if node.children.is_empty() {
                                        specified_leaves.insert(sub_id);
                                    }
                                }
                            }
                        }
                    }

                    // Identify all leaves in the tree
                    let all_leaves: BTreeSet<usize> = tree.get_leaves().into_iter().collect();

                    // Complement = All Leaves - Specified Leaves
                    let complement_leaves: Vec<usize> =
                        all_leaves.difference(&specified_leaves).cloned().collect();

                    if !complement_leaves.is_empty() {
                        let mut comp_nodes = complement_leaves.clone();
                        let mut comp_lca = comp_nodes.pop().unwrap();
                        for id in &comp_nodes {
                            comp_lca = tree.get_common_ancestor(&comp_lca, id).unwrap();
                        }

                        if comp_lca != old_root {
                            sub_root_id = comp_lca;
                        }
                    }
                }

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
            } else {
                // Default behavior: Root at the middle of the longest branch
                if let Some(longest_node) = tree.get_node_with_longest_edge() {
                    let new_root = tree.insert_parent(longest_node).unwrap();
                    tree.reroot_at(new_root, process_support).unwrap();
                    tree.remove_degree_two_nodes();
                }
            }
        }

        let out_string = tree.to_newick();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}
