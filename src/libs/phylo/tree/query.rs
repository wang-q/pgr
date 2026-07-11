use std::collections::{BTreeMap, BTreeSet};

use super::Tree;
use crate::libs::phylo::node::{Node, NodeId};
use clap::ArgMatches;
use regex::RegexBuilder;

/// Return the path from root to `id` (inclusive).
pub fn get_path_from_root(tree: &Tree, id: &NodeId) -> anyhow::Result<Vec<NodeId>> {
    let mut path = Vec::new();
    let mut current = *id;

    if tree.get_node(current).is_none() {
        anyhow::bail!("Node {} not found", current);
    }

    loop {
        path.push(current);
        match tree.nodes[current].parent {
            Some(p) => current = p,
            None => break,
        }
    }

    path.reverse();
    // Validate root
    if let Some(root) = tree.root {
        if path[0] != root {
            anyhow::bail!("Node is detached from root");
        }
    }

    Ok(path)
}

/// Find Lowest Common Ancestor (LCA) of two nodes.
pub fn get_common_ancestor(tree: &Tree, a: &NodeId, b: &NodeId) -> anyhow::Result<NodeId> {
    let path_a = get_path_from_root(tree, a)?;
    let path_b = get_path_from_root(tree, b)?;

    let mut lca = None;

    for (u, v) in path_a.iter().zip(path_b.iter()) {
        if u == v {
            lca = Some(*u);
        } else {
            break;
        }
    }

    lca.ok_or_else(|| anyhow::anyhow!("Nodes are not in the same tree (no common ancestor)"))
}

/// Find the Lowest Common Ancestor (LCA) of multiple nodes.
pub fn get_lca(tree: &Tree, nodes: &[NodeId]) -> anyhow::Result<NodeId> {
    if nodes.is_empty() {
        anyhow::bail!("Cannot find LCA of empty node set");
    }
    let mut lca = nodes[0];
    for &n in &nodes[1..] {
        lca = get_common_ancestor(tree, &lca, &n)?;
    }
    Ok(lca)
}

/// Calculate distance between two nodes.
/// Returns (weighted_distance, topological_distance).
pub fn get_distance(tree: &Tree, a: &NodeId, b: &NodeId) -> anyhow::Result<(f64, usize)> {
    let lca = get_common_ancestor(tree, a, b)?;

    let dist_to_lca = |start: &NodeId, end: &NodeId| -> (f64, usize) {
        let mut weighted = 0.0;
        let mut topo = 0;
        let mut curr = *start;

        while curr != *end {
            if let Some(node) = tree.get_node(curr) {
                weighted += node.length.unwrap_or(0.0);
                topo += 1;
                if let Some(p) = node.parent {
                    curr = p;
                } else {
                    break;
                }
            }
        }
        (weighted, topo)
    };

    let (w1, t1) = dist_to_lca(a, &lca);
    let (w2, t2) = dist_to_lca(b, &lca);

    Ok((w1 + w2, t1 + t2))
}

/// Distance between two nodes. Uses branch lengths if non-zero, else edge count.
pub fn node_distance(tree: &Tree, a: &NodeId, b: &NodeId) -> anyhow::Result<f64> {
    let (edge_sum, num_edges) = get_distance(tree, a, b)?;
    Ok(if edge_sum.abs() > 1e-9 {
        edge_sum
    } else {
        num_edges as f64
    })
}

/// Check if a set of nodes is monophyletic.
/// A set is monophyletic if the set of leaves in the subtree of their LCA
/// is exactly the same as the set of leaves reachable from the input nodes.
pub fn is_monophyletic(tree: &Tree, nodes: &[NodeId]) -> bool {
    if nodes.is_empty() {
        return false;
    }
    if nodes.len() == 1 {
        // A single node is monophyletic with respect to itself
        return true;
    }

    // 1. Find LCA
    let mut lca = nodes[0];
    for &n in &nodes[1..] {
        match get_common_ancestor(tree, &lca, &n) {
            Ok(anc) => lca = anc,
            Err(_) => return false, // Not in same tree
        }
    }

    // 2. Get all leaves under LCA
    let lca_leaves: BTreeSet<NodeId> = crate::libs::phylo::tree::stat::get_leaves(tree, lca)
        .into_iter()
        .collect();

    // 3. Get all leaves from input nodes
    let mut input_leaves = BTreeSet::new();
    for &n in nodes {
        let leaves = crate::libs::phylo::tree::stat::get_leaves(tree, n);
        for leaf in leaves {
            input_leaves.insert(leaf);
        }
    }

    // 4. Compare
    lca_leaves == input_leaves
}

/// Collect IDs of all named leaves (children.is_empty() && name.is_some()) in the subtree rooted at `id`.
pub fn get_named_leaves(tree: &Tree, id: NodeId) -> BTreeSet<NodeId> {
    let mut result = BTreeSet::new();
    let subtree_nodes = tree.get_subtree(&id);
    for nid in subtree_nodes {
        if let Some(node) = tree.get_node(nid) {
            if node.children.is_empty() && node.name.is_some() {
                result.insert(nid);
            }
        }
    }
    result
}

/// Get height of a node (max distance to any leaf in its subtree).
pub fn get_height(tree: &Tree, id: NodeId, weighted: bool) -> f64 {
    let node = match tree.get_node(id) {
        Some(n) => n,
        None => return 0.0,
    };

    if node.children.is_empty() {
        return 0.0;
    }

    node.children
        .iter()
        .map(|&child| {
            let dist = if weighted {
                tree.get_node(child).and_then(|n| n.length).unwrap_or(0.0)
            } else {
                1.0
            };
            dist + get_height(tree, child, weighted)
        })
        .fold(0.0, f64::max)
}

/// Count number of descendants (all nodes in subtree excluding self).
pub fn count_descendants(tree: &Tree, id: NodeId) -> usize {
    let mut count = 0;
    if let Some(node) = tree.get_node(id) {
        for &child in &node.children {
            count += 1 + count_descendants(tree, child);
        }
    }
    count
}

/// Find nodes matching a predicate.
pub fn find_nodes<F>(tree: &Tree, predicate: F) -> Vec<NodeId>
where
    F: Fn(&Node) -> bool,
{
    tree.nodes
        .iter()
        .filter(|n| !n.deleted && predicate(n))
        .map(|n| n.id)
        .collect()
}

/// Get node ID by name. Returns first match.
pub fn get_node_by_name(tree: &Tree, name: &str) -> Option<NodeId> {
    tree.nodes
        .iter()
        .find(|n| !n.deleted && n.name.as_deref() == Some(name))
        .map(|n| n.id)
}

/// Get root ID.
pub fn get_root(tree: &Tree) -> Option<NodeId> {
    tree.root
}

/// Find the medoid index among `ids` (min sum of pairwise branch-length distances).
pub fn tree_medoid(tree: &Tree, ids: &[NodeId]) -> Option<usize> {
    if ids.is_empty() {
        return None;
    }
    if ids.len() == 1 {
        return Some(0);
    }
    let mut min_sum_dist = f64::MAX;
    let mut best_idx = 0;
    for i in 0..ids.len() {
        let mut current_sum = 0.0;
        for j in 0..ids.len() {
            if i == j {
                continue;
            }
            let dist = get_distance(tree, &ids[i], &ids[j])
                .map(|(d, _)| d)
                .unwrap_or(f64::MAX);
            current_sum += dist;
        }
        if current_sum < min_sum_dist {
            min_sum_dist = current_sum;
            best_idx = i;
        }
    }
    Some(best_idx)
}

/// Lax-mode complement LCA: when the LCA of `specified_ids` equals `root_id`,
/// treat the unspecified leaves as the ingroup and return their LCA.
///
/// Returns `Some(comp_lca)` when the complement is non-empty and its LCA
/// differs from `root_id`; otherwise returns `None`.
pub fn lax_complement_lca(
    tree: &Tree,
    specified_ids: &BTreeSet<NodeId>,
    root_id: NodeId,
) -> Option<NodeId> {
    // Leaves under each specified node
    let mut specified_leaves = BTreeSet::new();
    for &id in specified_ids {
        let subtree = tree.get_subtree(&id);
        for sub_id in subtree {
            if let Some(node) = tree.get_node(sub_id) {
                if node.children.is_empty() {
                    specified_leaves.insert(sub_id);
                }
            }
        }
    }

    let all_leaves: BTreeSet<NodeId> = tree.get_leaves().into_iter().collect();
    let complement_leaves: Vec<NodeId> =
        all_leaves.difference(&specified_leaves).cloned().collect();

    if complement_leaves.is_empty() {
        return None;
    }

    let mut comp_nodes = complement_leaves.clone();
    let mut comp_lca = comp_nodes.pop()?;
    for id in &comp_nodes {
        comp_lca = tree.get_common_ancestor(&comp_lca, id).ok()?;
    }

    if comp_lca == root_id {
        None
    } else {
        Some(comp_lca)
    }
}

/// Format a node's label with extra columns (`dup`, `taxid`, `species`, `full`).
pub fn format_label_columns(node: &Node, name: &str, columns: &[String]) -> String {
    let mut out = String::from(name);
    if columns.is_empty() {
        return out;
    }
    for column in columns {
        match column.as_str() {
            "dup" => out.push_str(&format!("\t{}", name)),
            "taxid" => out.push_str(&format!(
                "\t{}",
                node.get_property("T").map(|s| s.as_str()).unwrap_or("")
            )),
            "species" => out.push_str(&format!(
                "\t{}",
                node.get_property("S").map(|s| s.as_str()).unwrap_or("")
            )),
            "full" => {
                let comment = node
                    .properties
                    .as_ref()
                    .map(|p| {
                        let pairs: Vec<String> =
                            p.iter().map(|(k, v)| format!(":{}={}", k, v)).collect();
                        if pairs.is_empty() {
                            String::new()
                        } else {
                            format!("[{}]", pairs.join(" "))
                        }
                    })
                    .unwrap_or_default();
                out.push_str(&format!("\t{}", comment));
            }
            _ => {}
        }
    }
    out
}

/// Returns IDs of named nodes matching the name selection rules from CLI args.
pub fn match_names(tree: &Tree, args: &ArgMatches) -> anyhow::Result<BTreeSet<usize>> {
    // IDs with names
    let id_of: BTreeMap<_, _> = tree.get_name_id();

    // all matched IDs
    let mut ids = BTreeSet::new();

    // ids supplied by --node
    if args.contains_id("node") {
        for name in args.get_many::<String>("node").unwrap() {
            if let Some(id) = id_of.get(name) {
                ids.insert(*id);
            }
        }
    }

    // ids supplied by --name-list
    if args.contains_id("name_list") {
        let file = args.get_one::<String>("name_list").unwrap();
        for name in crate::libs::io::read_names::<Vec<String>>(file)?.iter() {
            if let Some(id) = id_of.get(name) {
                ids.insert(*id);
            }
        }
    }

    // ids matched with --regex
    if args.contains_id("regex") {
        for regex in args.get_many::<String>("regex").unwrap() {
            let re = RegexBuilder::new(regex).case_insensitive(true).build()?;
            for name in id_of.keys() {
                if re.is_match(name) {
                    if let Some(id) = id_of.get(name) {
                        ids.insert(*id);
                    }
                }
            }
        }
    }

    // Default is printing all named nodes
    let is_all =
        !(args.contains_id("node") || args.contains_id("name_list") || args.contains_id("regex"));

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
            .filter_map(|e| tree.get_node(*e).cloned())
            .collect();
        for node in &nodes {
            if !node.is_leaf() {
                let subtree_ids = tree.get_subtree(&node.id);
                for id in &subtree_ids {
                    if let Some(n) = tree.get_node(*id) {
                        if n.name.is_some() {
                            ids.insert(*id);
                        }
                    }
                }
            }
        }
    }

    Ok(ids)
}

/// Returns IDs of nodes matching the position selection rules from CLI args.
pub fn match_positions(tree: &Tree, args: &ArgMatches) -> BTreeSet<usize> {
    let mut skip_internal = if args.try_contains_id("internal").is_ok() {
        args.get_flag("internal")
    } else {
        false
    };
    let skip_leaf = if args.try_contains_id("leaf").is_ok() {
        args.get_flag("leaf")
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

    let Some(root_id) = tree.get_root() else {
        return ids;
    };
    let preorder_ids = tree.preorder(&root_id);

    preorder_ids.iter().for_each(|id| {
        if let Some(node) = tree.get_node(*id) {
            if node.is_leaf() && !skip_leaf {
                ids.insert(*id);
            }
            if !node.is_leaf() && !skip_internal {
                ids.insert(*id);
            }
        }
    });

    ids
}
