//! Shared tree-depth helpers for SVG and Forest layout.

use super::super::Tree;
use crate::libs::phylo::node::NodeId;

/// Depth of `id` from the root (root has depth 0).
pub(super) fn node_depth(tree: &Tree, id: NodeId) -> usize {
    let mut depth = 0usize;
    let mut curr = id;
    while let Some(node) = tree.get_node(curr) {
        match node.parent {
            Some(p) => {
                depth += 1;
                curr = p;
            }
            None => break,
        }
    }
    depth
}

// max depth of this node's children
pub(super) fn branch_depth(tree: &Tree, id: NodeId) -> usize {
    let self_depth = node_depth(tree, id);
    match tree.get_subtree(&id) {
        Ok(nodes) => nodes
            .iter()
            .map(|nid| node_depth(tree, *nid))
            .max()
            .unwrap_or(self_depth),
        Err(_) => self_depth,
    }
}
