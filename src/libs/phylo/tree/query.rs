use std::collections::BTreeSet;

use super::Tree;
use crate::libs::phylo::node::{Node, NodeId};
pub fn get_path_from_root(tree: &Tree, id: &NodeId) -> Result<Vec<NodeId>, String> {
    let mut path = Vec::new();
    let mut current = *id;

    if tree.get_node(current).is_none() {
        return Err(format!("Node {} not found", current));
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
            return Err("Node is detached from root".to_string());
        }
    } else if !path.is_empty() {
        // Tree has no root but node exists? Should not happen in valid tree.
    }

    Ok(path)
}

/// Find Lowest Common Ancestor (LCA) of two nodes.
pub fn get_common_ancestor(tree: &Tree, a: &NodeId, b: &NodeId) -> Result<NodeId, String> {
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

    lca.ok_or_else(|| "Nodes are not in the same tree (no common ancestor)".to_string())
}

/// Calculate distance between two nodes.
/// Returns (weighted_distance, topological_distance).
pub fn get_distance(tree: &Tree, a: &NodeId, b: &NodeId) -> Result<(f64, usize), String> {
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
