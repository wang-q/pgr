use super::utils as nwr;
use clap::*;
use pgr::libs::phylo::reader;
use pgr::libs::phylo::writer;
use std::collections::BTreeSet;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("subtree")
        .about("Extracts a subtree")
        .after_help(
            r###"
Extracts a subtree (clade) rooted at the lowest common ancestor (LCA) of all
provided nodes.

Notes:
* Node selection:
    * `-n`: Select nodes by exact name.
    * `-f`: Select nodes from a file.
    * `-r`: Select nodes by regular expression.
    * If no selection is provided, no output is generated.
* Monophyly check (`-M`):
    * Ensures the subtree defined by the LCA contains ONLY the selected named
      terminal nodes.
    * Useful to verify if a group is monophyletic.
* Condense (`--condense`):
    * Instead of extracting the subtree, it replaces the subtree with a single
      node in the original tree.
    * The new node inherits the edge length of the subtree root.
    * Added annotations: `member=<count>` and `tri=white`.

Examples:
1. Extract subtree for Human and Chimp:
   pgr nwk subtree tree.nwk -n Human -n Chimp

2. Extract subtree for all nodes matching "Homo":
   pgr nwk subtree tree.nwk -r "^Homo"

3. Condense the Hominini clade (LCA of Homo and Pan) into a single node "Hominini":
   pgr nwk subtree tree.nwk -n Homo -n Pan --condense Hominini

4. Check if a group is monophyletic (prints nothing if not):
   pgr nwk subtree tree.nwk -n Human -n Chimp -n Gorilla -M
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
                .help("Select nodes by exact name"),
        )
        .arg(
            Arg::new("file")
                .long("file")
                .short('f')
                .num_args(1)
                .help("Select nodes from a file"),
        )
        .arg(
            Arg::new("regex")
                .long("regex")
                .short('r')
                .num_args(1)
                .action(ArgAction::Append)
                .help("Select nodes by regular expression (case insensitive)"),
        )
        .arg(
            Arg::new("descendants")
                .long("descendants")
                .short('D')
                .action(ArgAction::SetTrue)
                .help("Include all descendants of selected internal nodes"),
        )
        .arg(
            Arg::new("monophyly")
                .long("monophyly")
                .short('M')
                .action(ArgAction::SetTrue)
                .help("Only print the subtree when it's monophyletic"),
        )
        .arg(
            Arg::new("condense")
                .long("condense")
                .short('c')
                .num_args(1)
                .help("Condense the subtree into a single node with this name"),
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

    let is_monophyly = args.get_flag("monophyly");
    let condense_name = args.get_one::<String>("condense");
    let is_condense = condense_name.is_some();

    let infile = args.get_one::<String>("infile").unwrap();
    let trees = reader::from_file(infile)?;

    if trees.is_empty() {
        return Ok(());
    }

    // We process the first tree only, as typical for this operation
    // or should we process all? The original code did `nwr::read_newick` which likely returned one tree.
    // Let's iterate but usually nwk tools handle multi-tree files by processing each.
    // But `condense` modifies the tree.
    // Let's assume we process all trees.

    // Note: Since we might modify the tree (condense), we need it to be mutable.
    // But we are iterating.
    // We can iterate and process.

    // Since `reader::from_file` returns Vec<Tree>, we can iterate mutably.
    let mut trees = trees;

    for tree in &mut trees {
        // IDs matching names
        let ids = nwr::match_names(tree, args);

        if ids.is_empty() {
            continue;
        }

        // Find LCA
        let mut nodes: Vec<usize> = ids.iter().cloned().collect();
        let mut sub_root_id = nodes.pop().unwrap();

        // If multiple nodes, find their common ancestor
        for id in &nodes {
            match tree.get_common_ancestor(&sub_root_id, id) {
                Ok(anc) => sub_root_id = anc,
                Err(_) => {
                    // Nodes might be in disjoint parts if tree structure is broken or multiple roots?
                    // Tree struct assumes single root for `get_common_ancestor` usually?
                    // Or if they are disconnected.
                    // If error, we can't find LCA.
                    continue;
                }
            }
        }

        // Monophyly check
        if is_monophyly {
            let mut descendants_named = BTreeSet::new();
            if let Ok(subtree_nodes) = tree.get_subtree(&sub_root_id) {
                for id in subtree_nodes {
                    let node = tree.get_node(id).unwrap();
                    // Check if it's a leaf/tip and has a name
                    if node.children.is_empty() && node.name.is_some() {
                        descendants_named.insert(id);
                    }
                }
            }

            // Check if the selected set matches the set of named descendants
            // The logic is: the subtree should contain EXACTLY the selected named tips.
            // If `ids` contains internal nodes, we should probably only compare tips?
            // Original code:
            // `if name_of.contains_key(id) && tree.get(id).unwrap().is_tip()`
            // `ids.ne(&descendants)`
            // So it compares the set of selected IDs (which presumably are all named, or descendants of named)
            // with the set of named tips in the subtree.

            // Wait, `match_names` might include internal nodes if `-D` is used.
            // If `-D` is used, `ids` includes descendants.
            // The original code used `id_of` (map name->id) to check if a node "has a name" in the original map context?
            // Actually `nwr::get_name_id` returns map of named nodes.
            // So `name_of.contains_key(id)` effectively checks `node.name.is_some()`.

            // Let's filter `ids` to tips only for comparison?
            // Or just compare sets as is.
            // If I selected an internal node explicitly, it is in `ids`.
            // But `descendants` collection logic in original code: `is_tip()`.
            // So original code compares `ids` (which might contain internal nodes?) with `descendants` (tips only).
            // This implies `ids` should only contain tips for this check to pass?
            // Unless `match_names` only returns tips? No, it returns whatever matches.
            // But `match_names` without `-D` returns the named nodes themselves.
            // If I select an internal node "A", `ids` has "A".
            // `descendants` has tips of "A".
            // "A" != tips of "A". So monophyly check fails for internal node selection?
            // That seems strict.
            // But if `-D` is used, `ids` includes descendants.
            // If I use `-n Human -n Chimp`, `ids` = {Human, Chimp}.
            // LCA is Ancestor. Subtree tips are {Human, Chimp}. Match!
            // If I use `-n Ancestor`, `ids` = {Ancestor}.
            // LCA is Ancestor. Subtree tips {Human, Chimp}.
            // {Ancestor} != {Human, Chimp}. Fails.
            // This seems to be the intended behavior of the original code.

            if ids != descendants_named {
                if is_condense {
                    // Even if monophyly fails, if condensing, we output the original tree?
                    // Original code:
                    // `if ids.ne(&descendants) { if is_condense { print tree } return Ok(()); }`
                    // So it aborts the operation and prints the UNMODIFIED tree.
                    let out_string = writer::write_newick(tree);
                    writer.write_fmt(format_args!("{}\n", out_string)).unwrap();
                }
                continue;
            }
        }

        // Output
        if is_condense {
            let name = condense_name.unwrap();
            
            // 1. Get info from sub_root
            let sub_root = tree.get_node(sub_root_id).unwrap();
            let parent_id_opt = sub_root.parent;
            let edge_len = sub_root.length;
            
            if let Some(parent_id) = parent_id_opt {
                // 2. Create new node
                let new_node_id = tree.add_node();
                if let Some(node) = tree.get_node_mut(new_node_id) {
                    node.set_name(name);
                    node.length = edge_len;
                    // Add properties
                    let mut props = std::collections::BTreeMap::new();
                    props.insert("member".to_string(), ids.len().to_string());
                    props.insert("tri".to_string(), "white".to_string());
                    node.properties = Some(props);
                }
                
                // 3. Remove old subtree
                tree.remove_node(sub_root_id, true);
                
                // 4. Link new node to parent
                // Note: remove_node disconnects sub_root from parent, so we can just add child.
                tree.add_child(parent_id, new_node_id).unwrap();
                
                // 5. Output full tree
                let out_string = writer::write_newick(tree);
                writer.write_fmt(format_args!("{}\n", out_string)).unwrap();
            } else {
                // Subtree root is tree root.
                // Replaces the entire tree with a single node?
                // Logic:
                // Clear tree? Or just make a new root.
                // Since we want to output the "condensed tree", which is just one node.
                // We can just construct a string or modify tree.
                // Let's modify tree for consistency.
                
                // Remove root (clears everything basically)
                tree.remove_node(sub_root_id, true);
                
                // Add new root
                let new_root = tree.add_node();
                tree.set_root(new_root);
                if let Some(node) = tree.get_node_mut(new_root) {
                    node.set_name(name);
                    // Root usually has no length
                    // Add properties
                     let mut props = std::collections::BTreeMap::new();
                    props.insert("member".to_string(), ids.len().to_string());
                    props.insert("tri".to_string(), "white".to_string());
                    node.properties = Some(props);
                }
                
                let out_string = writer::write_newick(tree);
                writer.write_fmt(format_args!("{}\n", out_string)).unwrap();
            }

        } else {
            // Extract subtree
            let out_string = writer::write_subtree(tree, sub_root_id, "");
            writer.write_fmt(format_args!("{}\n", out_string)).unwrap();
        }
    }

    Ok(())
}
