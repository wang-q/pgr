use super::utils as nwr;
use clap::*;
use pgr::libs::phylo::reader;
use std::collections::HashSet;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("prune")
        .about("Remove nodes from the Newick file")
        .after_help(
            r###"
This tool removes nodes from the Newick tree based on provided labels or patterns.

Notes:
* Target nodes can be specified by name (`--node`), file (`--file`), or regex (`--regex`).
* With `--invert`, specified nodes (along with their ancestors and descendants) are kept, and everything else is removed.
* Topology changes:
    * If a node removal leaves its parent with only one child, the parent is collapsed (spliced out).
    * Internal nodes that lose all children are also removed.

Examples:
1. Remove specific nodes by name:
   $ pgr nwk prune input.nwk -n Homo -n Pan

2. Remove nodes using a list in a file:
   $ pgr nwk prune input.nwk -f remove.txt

3. Keep a clade (e.g., Hominidae) and remove everything else (Invert mode):
   $ pgr nwk prune input.nwk -v -n Hominidae

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
            Arg::new("file")
                .long("file")
                .short('f')
                .num_args(1)
                .help("A file contains node names"),
        )
        .arg(
            Arg::new("regex")
                .long("regex")
                .short('r')
                .num_args(1)
                .action(ArgAction::Append)
                .help("Nodes match the regular expression"),
        )
        .arg(
            Arg::new("descendants")
                .long("descendants")
                .short('D')
                .action(ArgAction::SetTrue)
                .help("Include all descendants of internal nodes"),
        )
        .arg(
            Arg::new("invert")
                .long("invert")
                .short('v')
                .action(ArgAction::SetTrue)
                .help("Invert pruning: keep specified nodes, their ancestors and descendants"),
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
    //----------------------------
    // Args
    //----------------------------
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let infile = args.get_one::<String>("infile").unwrap();

    let trees = reader::from_file(infile)?;

    for mut tree in trees {
        //----------------------------
        // Operating
        //----------------------------

        // 1. Identify internals before pruning
        // We store them to check later if they become leaves
        let mut old_internals = vec![];
        if let Some(root) = tree.get_root() {
            let all_nodes = tree.levelorder(&root).unwrap_or_default();
            for id in all_nodes {
                if let Some(node) = tree.get_node(id) {
                    if !node.children.is_empty() {
                        old_internals.push(id);
                    }
                }
            }
        }

        // 2. Identify targets
        let target_ids = nwr::match_names(&tree, args);

        // 3. Determine nodes to remove
        let to_remove = if args.get_flag("invert") {
            let mut keep = HashSet::new();

            if let Some(root) = tree.get_root() {
                // Convert targets to HashSet for fast lookup
                let target_set: HashSet<usize> = target_ids.iter().cloned().collect();
                let mut is_in_clade = HashSet::new();

                // Pass 1: Downward propagation (Descendants)
                // Use levelorder which is topological (parents before children)
                let all_nodes = tree.levelorder(&root).unwrap_or_default();

                for &id in &all_nodes {
                    let mut kept_descendant = false;

                    // Check if self is target
                    if target_set.contains(&id) {
                        kept_descendant = true;
                    }
                    // Check if parent is in clade (propagates downwards)
                    else if let Some(node) = tree.get_node(id) {
                        if let Some(parent) = node.parent {
                            if is_in_clade.contains(&parent) {
                                kept_descendant = true;
                            }
                        }
                    }

                    if kept_descendant {
                        is_in_clade.insert(id);
                        keep.insert(id);
                    }
                }

                // Pass 2: Upward propagation (Ancestors)
                // Iterate in reverse (children before parents)
                for &id in all_nodes.iter().rev() {
                    if keep.contains(&id) {
                        if let Some(node) = tree.get_node(id) {
                            if let Some(parent) = node.parent {
                                keep.insert(parent);
                            }
                        }
                    }
                }

                // Collect nodes NOT in keep set
                all_nodes
                    .into_iter()
                    .filter(|id| !keep.contains(id))
                    .collect()
            } else {
                vec![]
            }
        } else {
            target_ids.into_iter().collect()
        };

        // 4. Remove nodes
        for id in to_remove {
            tree.remove_node(id, true);
        }

        // 5. Clean up: remove internals that became leaves
        for id in old_internals.into_iter().rev() {
            if let Some(node) = tree.get_node(id) {
                // If it's still there and has no children, it became a leaf
                // (Only remove if it wasn't originally a leaf - checked by old_internals logic)
                if node.children.is_empty() {
                    tree.remove_node(id, true);
                }
            }
        }

        // 6. Cleanup degree-2 nodes (Post-order)
        if let Some(root) = tree.get_root() {
            let nodes = tree.postorder(&root).unwrap_or_default();
            for id in nodes {
                if let Some(node) = tree.get_node(id) {
                    if node.children.len() == 1 {
                        if tree.get_root() == Some(id) {
                            // Root with 1 child -> promote child to root
                            let child_id = node.children[0];
                            tree.set_root(child_id);
                            tree.remove_node(id, false);
                        } else {
                            // Internal degree-2 -> collapse
                            tree.collapse_node(id).ok();
                        }
                    }
                }
            }
        }

        //----------------------------
        // Output
        //----------------------------
        let out_string = tree.to_newick();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}
