use clap::ArgMatches;
use pgr::libs::phylo::node::Node;
use pgr::libs::phylo::tree::Tree;
use regex::RegexBuilder;
use std::collections::{BTreeMap, BTreeSet};

// Named IDs that match the name rules
pub fn match_names(tree: &Tree, args: &ArgMatches) -> BTreeSet<usize> {
    // IDs with names
    let id_of: BTreeMap<_, _> = tree.get_name_id();

    // all matched IDs
    let mut ids = BTreeSet::new();

    // ids supplied by --node
    if args.contains_id("node") {
        for name in args.get_many::<String>("node").unwrap() {
            if id_of.contains_key(name) {
                let id = id_of.get(name).unwrap();
                ids.insert(*id);
            }
        }
    }

    // ids supplied by --file
    if args.contains_id("file") {
        let file = args.get_one::<String>("file").unwrap();
        for name in intspan::read_first_column(file).iter() {
            if id_of.contains_key(name) {
                let id = id_of.get(name).unwrap();
                ids.insert(*id);
            }
        }
    }

    // ids matched with --regex
    if args.contains_id("regex") {
        for regex in args.get_many::<String>("regex").unwrap() {
            let re = RegexBuilder::new(regex)
                .case_insensitive(true)
                .build()
                .unwrap();
            for name in id_of.keys() {
                if re.is_match(name) {
                    let id = id_of.get(name).unwrap();
                    ids.insert(*id);
                }
            }
        }
    }

    // Default is printing all named nodes
    let is_all =
        !(args.contains_id("node") || args.contains_id("file") || args.contains_id("regex"));

    if is_all {
        ids = id_of.values().cloned().collect();
    }

    // Include all descendants of internal nodes
    let is_descendants = if args.try_contains_id("descendants").is_ok() {
        args.get_flag("descendants")
    } else {
        false
    };

    if is_descendants {
        let nodes: Vec<Node> = ids
            .iter()
            .map(|e| tree.get_node(*e).unwrap().clone())
            .collect();
        for node in &nodes {
            if !node.is_leaf() {
                for id in &tree.get_subtree(&node.id).unwrap() {
                    if tree.get_node(*id).unwrap().name.is_some() {
                        ids.insert(*id);
                    }
                }
            }
        }
    }

    ids
}

// IDs that match the position rules
pub fn match_positions(tree: &Tree, args: &ArgMatches) -> BTreeSet<usize> {
    let mut skip_internal = if args.try_contains_id("Internal").is_ok() {
        args.get_flag("Internal")
    } else {
        false
    };
    let skip_leaf = if args.try_contains_id("Leaf").is_ok() {
        args.get_flag("Leaf")
    } else {
        false
    };

    let is_monophyly = if args.try_contains_id("monophyly").is_ok() {
        args.get_flag("monophyly")
    } else {
        false
    };

    if is_monophyly {
        skip_internal = true;
    }

    // all matched IDs
    let mut ids = BTreeSet::new();

    // inorder needs IsBinary
    tree.preorder(&tree.get_root().unwrap())
        .unwrap()
        .iter()
        .for_each(|id| {
            let node = tree.get_node(*id).unwrap();
            if node.is_leaf() && !skip_leaf {
                ids.insert(*id);
            }
            if !node.is_leaf() && !skip_internal {
                ids.insert(*id);
            }
        });

    ids
}
