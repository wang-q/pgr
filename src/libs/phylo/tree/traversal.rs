use super::Tree;
use crate::libs::phylo::node::NodeId;
use std::collections::VecDeque;

/// Get node IDs in preorder traversal (Root -> Children)
pub fn preorder(tree: &Tree, start_node: NodeId) -> Vec<NodeId> {
    let mut result = Vec::new();
    let mut stack = vec![start_node];

    while let Some(id) = stack.pop() {
        if let Some(node) = tree.get_node(id) {
            result.push(id);
            // Push children in reverse order so they are processed in order
            for &child in node.children.iter().rev() {
                stack.push(child);
            }
        }
    }

    result
}

/// Get node IDs in postorder traversal (Children -> Root)
pub fn postorder(tree: &Tree, start_node: NodeId) -> Vec<NodeId> {
    let mut result = Vec::new();
    fn helper(tree: &Tree, id: NodeId, result: &mut Vec<NodeId>) {
        if let Some(node) = tree.get_node(id) {
            for &child in &node.children {
                helper(tree, child, result);
            }
            result.push(id);
        }
    }

    helper(tree, start_node, &mut result);
    result
}

/// Get node IDs in levelorder traversal (BFS)
pub fn levelorder(tree: &Tree, start_node: NodeId) -> Vec<NodeId> {
    let mut result = Vec::new();
    let mut queue = VecDeque::new();
    queue.push_back(start_node);

    while let Some(id) = queue.pop_front() {
        if let Some(node) = tree.get_node(id) {
            result.push(id);
            for &child in &node.children {
                queue.push_back(child);
            }
        }
    }

    result
}

/// Extract a subtree rooted at `node_id`.
/// Returns a new Tree.
pub fn extract_subtree(tree: &Tree, node_id: NodeId) -> anyhow::Result<Tree> {
    if tree.get_node(node_id).is_none() {
        anyhow::bail!("Node {} not found", node_id);
    }

    let mut new_tree = Tree::new();
    let mut id_map = std::collections::HashMap::new();
    let mut stack = vec![(node_id, None::<NodeId>)];

    while let Some((old_id, new_parent_opt)) = stack.pop() {
        let Some(old_node) = tree.get_node(old_id) else {
            // Skip deleted nodes that are still referenced as children.
            continue;
        };

        let new_id = new_tree.add_node();
        id_map.insert(old_id, new_id);

        if let Some(new_node) = new_tree.get_node_mut(new_id) {
            new_node.name = old_node.name.clone();
            new_node.length = old_node.length;
            new_node.properties = old_node.properties.clone();
        }

        if let Some(new_parent) = new_parent_opt {
            super::ops::add_child(&mut new_tree, new_parent, new_id)?;
        } else {
            new_tree.set_root(new_id);
        }

        for &child in old_node.children.iter().rev() {
            stack.push((child, Some(new_id)));
        }
    }

    Ok(new_tree)
}
