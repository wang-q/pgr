use super::Tree;
use crate::libs::phylo::node::NodeId;

/// Add a child to a parent node.
/// Updates both parent's `children` list and child's `parent` field.
pub fn add_child(tree: &mut Tree, parent_id: NodeId, child_id: NodeId) -> Result<(), String> {
    // Validation
    if parent_id == child_id {
        return Err("Cannot add node as child of itself".to_string());
    }
    if tree.get_node(parent_id).is_none() {
        return Err(format!("Parent node {} not found or deleted", parent_id));
    }
    if tree.get_node(child_id).is_none() {
        return Err(format!("Child node {} not found or deleted", child_id));
    }

    // Check if child already has a parent
    let child_parent = tree.nodes[child_id].parent;
    if let Some(old_parent) = child_parent {
        return Err(format!(
            "Node {} already has parent {}",
            child_id, old_parent
        ));
    }

    // Link
    tree.nodes[child_id].parent = Some(parent_id);
    tree.nodes[parent_id].children.push(child_id);

    Ok(())
}

/// Soft remove a node and its descendants (optional recursive).
/// If recursive is false, children are orphaned (parent set to None).
pub fn remove_node(tree: &mut Tree, id: NodeId, recursive: bool) {
    if id >= tree.nodes.len() || tree.nodes[id].deleted {
        return;
    }

    // 1. Handle Parent Relation
    if let Some(parent_id) = tree.nodes[id].parent {
        // Remove self from parent's children list
        if let Some(parent) = tree.get_node_mut(parent_id) {
            parent.children.retain(|&child| child != id);
        }
    }

    // 2. Handle Children
    let children = tree.nodes[id].children.clone();
    for child_id in children {
        if recursive {
            remove_node(tree, child_id, true);
        } else {
            // Orphan the child
            if let Some(child) = tree.get_node_mut(child_id) {
                child.parent = None;
            }
        }
    }

    // 3. Mark as deleted
    if let Some(node) = tree.get_node_mut(id) {
        node.deleted = true;
        node.children.clear();
        node.parent = None;
    }

    // 4. Update root if needed
    if tree.root == Some(id) {
        tree.root = None;
    }
}

/// Collapse a node, removing it and connecting its children to its parent.
/// Edge lengths are summed (parent->node + node->child).
pub fn collapse_node(tree: &mut Tree, id: NodeId) -> Result<(), String> {
    if tree.get_node(id).is_none() {
        return Err(format!("Node {} not found", id));
    }
    if tree.root == Some(id) {
        return Err("Cannot collapse root node".to_string());
    }

    // 1. Get info
    let (parent_id, parent_edge) = {
        let node = tree.get_node(id).unwrap();
        // Safety: Checked root above, so parent must exist
        (node.parent.unwrap(), node.length)
    };

    let children_info: Vec<(NodeId, Option<f64>)> = {
        let node = tree.get_node(id).unwrap();
        node.children
            .iter()
            .map(|&c| {
                let child_node = tree.nodes.get(c).unwrap();
                (c, child_node.length)
            })
            .collect()
    };

    // 2. Re-parent children
    let mut new_children_ids = Vec::new();
    for (child_id, child_edge) in children_info {
        let new_edge = match (parent_edge, child_edge) {
            (Some(p), Some(c)) => Some(p + c),
            (Some(p), None) => Some(p),
            (None, Some(c)) => Some(c),
            (None, None) => None,
        };

        // Update child
        if let Some(child) = tree.get_node_mut(child_id) {
            child.parent = Some(parent_id);
            child.length = new_edge;
        }
        new_children_ids.push(child_id);
    }

    // 3. Update parent's children list
    if let Some(parent) = tree.get_node_mut(parent_id) {
        if let Some(pos) = parent.children.iter().position(|&x| x == id) {
            parent.children.splice(pos..pos + 1, new_children_ids);
        }
    }

    // 4. Mark deleted
    if let Some(node) = tree.get_node_mut(id) {
        node.deleted = true;
        node.children.clear();
        node.parent = None;
    }

    Ok(())
}

/// Compact the tree by removing soft-deleted nodes and remapping IDs.
/// This invalidates all existing NodeIds held outside!
pub fn compact(tree: &mut Tree) {
    let mut old_to_new = std::collections::HashMap::new();
    let mut new_nodes = Vec::with_capacity(tree.nodes.len());
    let mut new_idx = 0;

    // 1. Build mapping and new node list (without edges first)
    for old_node in &tree.nodes {
        if !old_node.deleted {
            old_to_new.insert(old_node.id, new_idx);
            // Create a shallow copy with updated ID but empty edges (will fill later)
            let mut new_node = old_node.clone();
            new_node.id = new_idx;
            new_node.parent = None; // Reset relations
            new_node.children.clear();
            new_nodes.push(new_node);
            new_idx += 1;
        }
    }

    // 2. Reconstruct edges using the mapping
    for (old_idx, node) in tree.nodes.iter().enumerate() {
        if node.deleted {
            continue;
        }

        let new_self_idx = *old_to_new.get(&old_idx).unwrap();

        // Remap parent
        if let Some(old_parent) = node.parent {
            if let Some(&new_parent) = old_to_new.get(&old_parent) {
                new_nodes[new_self_idx].parent = Some(new_parent);
            }
        }

        // Remap children
        for &old_child in &node.children {
            if let Some(&new_child) = old_to_new.get(&old_child) {
                new_nodes[new_self_idx].children.push(new_child);
            }
        }
    }

    // 3. Update root
    if let Some(old_root) = tree.root {
        tree.root = old_to_new.get(&old_root).copied();
    }

    // 4. Swap
    tree.nodes = new_nodes;
}

/// Insert a node in the middle of the desired node and its parent.
/// Returns the new parent node ID.
pub fn insert_parent(tree: &mut Tree, id: NodeId) -> Result<NodeId, String> {
    let node = tree.get_node(id).ok_or(format!("Node {} not found", id))?;
    let parent = node.parent.ok_or("Node has no parent")?;
    let length = node.length;
    let new_length = length.map(|l| l / 2.0);

    let new_node = tree.add_node();

    // Link parent -> new_node
    add_child(tree, parent, new_node)?;
    if let Some(n) = tree.get_node_mut(new_node) {
        n.length = new_length;
    }

    // Unlink parent -> id
    if let Some(p_node) = tree.get_node_mut(parent) {
        p_node.children.retain(|&c| c != id);
    }
    // Update id parent
    if let Some(node) = tree.get_node_mut(id) {
        node.parent = None;
    }

    // Link new_node -> id
    add_child(tree, new_node, id)?;
    if let Some(node) = tree.get_node_mut(id) {
        node.length = new_length;
    }

    Ok(new_node)
}

/// Swap parent-child link of a node.
/// This reverses the edge between the node and its parent.
pub fn swap_parent(
    tree: &mut Tree,
    id: NodeId,
    _prev_edge: Option<f64>,
) -> Result<Option<f64>, String> {
    let node = tree.get_node(id).ok_or(format!("Node {} not found", id))?;
    let parent = node.parent.ok_or("Node has no parent")?;

    // Swap lengths
    let child_len = node.length;
    let parent_len = tree.get_node(parent).ok_or("Parent not found")?.length;

    if let Some(p_node) = tree.get_node_mut(parent) {
        p_node.length = child_len;
    }
    if let Some(node) = tree.get_node_mut(id) {
        node.length = parent_len;
    }

    // Unlink parent -> id
    if let Some(p_node) = tree.get_node_mut(parent) {
        p_node.children.retain(|&c| c != id);
    }
    // Unlink id -> parent
    if let Some(node) = tree.get_node_mut(id) {
        node.parent = None;
    }

    // Link id -> parent (parent becomes child)
    // We must clear parent's parent pointer to satisfy add_child check (as parent is "moving down")
    if let Some(p_node) = tree.get_node_mut(parent) {
        p_node.parent = None;
    }

    add_child(tree, id, parent)?;

    Ok(None)
}

/// Insert a new parent node for a pair of nodes (LCA-based).
/// Returns the new parent node ID.
pub fn insert_parent_pair(tree: &mut Tree, id1: NodeId, id2: NodeId) -> Result<NodeId, String> {
    let old = tree.get_common_ancestor(&id1, &id2)?;

    // Get original edge lengths
    let edge1 = tree.get_node(id1).ok_or("Node 1 not found")?.length;
    let edge2 = tree.get_node(id2).ok_or("Node 2 not found")?.length;

    // New node with parent (old) has no edge length
    let new = tree.add_node();
    add_child(tree, old, new)?;

    // Move children to new node
    // 1. Unlink from their current parents
    let p1 = tree.get_node(id1).and_then(|n| n.parent);
    if let Some(p) = p1 {
        if let Some(p_node) = tree.get_node_mut(p) {
            p_node.children.retain(|&c| c != id1);
        }
    }

    let p2 = tree.get_node(id2).and_then(|n| n.parent);
    if let Some(p) = p2 {
        if let Some(p_node) = tree.get_node_mut(p) {
            p_node.children.retain(|&c| c != id2);
        }
    }

    if let Some(node) = tree.get_node_mut(id1) {
        node.parent = None;
    }
    if let Some(node) = tree.get_node_mut(id2) {
        node.parent = None;
    }

    // 2. Link to new
    add_child(tree, new, id1)?;
    if let Some(node) = tree.get_node_mut(id1) {
        node.length = edge1;
    }

    add_child(tree, new, id2)?;
    if let Some(node) = tree.get_node_mut(id2) {
        node.length = edge2;
    }

    Ok(new)
}

/// Remove nodes that have a parent and exactly one child (degree 2 nodes).
/// This is often used after rerooting to clean up the tree.
pub fn remove_degree_two_nodes(tree: &mut Tree) {
    loop {
        // Find a node that is:
        // 1. Not deleted
        // 2. Has a parent (not root)
        // 3. Has exactly 1 child
        let to_remove = if let Some(_root) = tree.get_root() {
            tree.find_nodes(|n| {
                n.parent.is_some() && // Not root
                n.children.len() == 1 // Degree 2 (1 parent, 1 child)
            })
            .first()
            .cloned()
        } else {
            None
        };

        if let Some(id) = to_remove {
            // Ignore result, just proceed
            let _ = collapse_node(tree, id);
        } else {
            break;
        }
    }
}

/// Deroot the tree by splicing out one of the root's children if the root is bifurcating.
/// This effectively merges the two edges connected to the root into a single edge,
/// removing the root node's structural role and making the tree multifurcating at the top level.
/// The "heavier" child (with more descendants) is the one collapsed into the root.
pub fn deroot(tree: &mut Tree) -> Result<(), String> {
    let root = tree.root.ok_or("Empty tree")?;
    let children = tree.get_node(root).unwrap().children.clone();

    if children.len() != 2 {
        return Err("Root is not bifurcating (degree != 2)".to_string());
    }

    let c1 = children[0];
    let c2 = children[1];

    // Weight = 1 (self) + descendants
    let weight1 = 1 + super::query::count_descendants(tree, c1);
    let weight2 = 1 + super::query::count_descendants(tree, c2);

    // Collapse the heavier one. If equal, pick first (c1).
    let target = if weight1 >= weight2 { c1 } else { c2 };

    collapse_node(tree, target)
}

/// Reroot the tree at the specified node.
/// This reverses the direction of edges along the path from the old root to the new root.
pub fn reroot_at(
    tree: &mut Tree,
    new_root_id: NodeId,
    process_support_values: bool,
) -> Result<(), String> {
    if tree.get_node(new_root_id).is_none() {
        return Err(format!("Node {} not found", new_root_id));
    }

    let old_root_id = tree.root.ok_or("Tree has no root")?;
    if old_root_id == new_root_id {
        return Ok(());
    }

    // 1. Get path from old root to new root
    let path = tree.get_path_from_root(&new_root_id)?;

    // 1.5 Process Support Values (Labels)
    // Shift internal node labels along the path to align with edge reversals
    if process_support_values {
        let new_root_is_leaf = tree
            .get_node(new_root_id)
            .map(|n| n.children.is_empty())
            .unwrap_or(false);

        // Capture original names
        let names: Vec<Option<String>> = path
            .iter()
            .map(|&id| tree.get_node(id).unwrap().name.clone())
            .collect();

        for i in 0..path.len() {
            let node_id = path[i];
            // Only modify internal nodes (leaves keep Taxon names)
            // Note: All nodes on path except possibly the last one are ancestors, thus internal.
            let is_leaf = (i == path.len() - 1) && new_root_is_leaf;

            if !is_leaf {
                let new_name = if i < path.len() - 1 {
                    // Take from next node, UNLESS next node is a leaf (Taxon)
                    let next_is_leaf = (i + 1 == path.len() - 1) && new_root_is_leaf;
                    if next_is_leaf {
                        None
                    } else {
                        names[i + 1].clone()
                    }
                } else {
                    // New root (internal): takes label from old root
                    names[0].clone()
                };

                if let Some(node) = tree.get_node_mut(node_id) {
                    node.name = new_name;
                }
            }
        }
    }

    // 2. Collect edge lengths along the path
    // path[i]'s length represents edge (path[i-1] -> path[i])
    let mut lengths = Vec::new();
    for &id in &path {
        lengths.push(tree.get_node(id).unwrap().length);
    }

    // 3. Reverse edges
    for i in (1..path.len()).rev() {
        let child_id = path[i];
        let parent_id = path[i - 1];
        let length = lengths[i];

        // a. Remove child from parent's children
        if let Some(parent) = tree.nodes.get_mut(parent_id) {
            parent.children.retain(|&x| x != child_id);
        }

        // b. Add parent to child's children
        if let Some(child) = tree.nodes.get_mut(child_id) {
            child.children.push(parent_id);
        }

        // c. Update parent's parent pointer and length
        if let Some(parent) = tree.nodes.get_mut(parent_id) {
            parent.parent = Some(child_id);
            parent.length = length;
        }
    }

    // 4. Finalize new root
    if let Some(new_root) = tree.nodes.get_mut(new_root_id) {
        new_root.parent = None;
        new_root.length = None;
    }

    tree.root = Some(new_root_id);

    Ok(())
}

/// Prune nodes that match a predicate.
/// Warning: This removes the matching nodes AND their descendants.
pub fn prune_where<F>(tree: &mut Tree, predicate: F)
where
    F: Fn(&crate::libs::phylo::node::Node) -> bool,
{
    // We need to collect IDs first to avoid borrowing issues
    let to_remove: Vec<NodeId> = tree
        .nodes
        .iter()
        .filter(|n| !n.deleted && predicate(n))
        .map(|n| n.id)
        .collect();

    for id in to_remove {
        remove_node(tree, id, true);
    }
}
