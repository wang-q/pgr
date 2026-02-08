use super::tree::Tree;
use super::NodeId;
use std::collections::HashMap;

/// Sort the children of each node by their name (label).
///
/// # Arguments
///
/// * `tree` - The tree to modify.
/// * `descending` - If true, sort in descending order (Z-A).
///
/// # Example
/// ```
/// use pgr::libs::phylo::tree::Tree;
/// use pgr::libs::phylo::algo;
///
/// let mut tree = Tree::from_newick("(B,A);").unwrap();
/// algo::sort_by_name(&mut tree, false);
/// assert_eq!(tree.to_newick(), "(A,B);");
/// ```
pub fn sort_by_name(tree: &mut Tree, descending: bool) {
    if tree.is_empty() {
        return;
    }

    let root = tree.get_root().unwrap();
    // We can use levelorder to visit all nodes, but we just need to iterate all nodes.
    // However, to propagate sorting up the tree (inheriting labels), we should use post-order.
    // Or, we can just sort all nodes. For a single pass, if we just sort children, the order of visiting nodes matters
    // IF the sort criteria depends on the sorted order of children (which it does for recursive labeling).
    // So we MUST use post-order (leaves first, then parents).
    let ids = match tree.postorder(&root) {
        Ok(v) => v,
        Err(_) => return,
    };

    // Pre-collect names to avoid borrowing issues during sort
    let mut name_map: HashMap<NodeId, String> = HashMap::new();
    for &id in &ids {
        if let Some(node) = tree.get_node(id) {
            name_map.insert(id, node.name.clone().unwrap_or_default());
        }
    }

    for id in ids {
        // We need to re-borrow name_map for each node if we update it?
        // Actually, if a node is unnamed, we want to USE its children's names.
        // But we are iterating post-order. Children are already visited and sorted.
        // So when we are at `id`, its children are sorted.
        // We can determine the "effective name" of `id` now.

        // Let's get the node
        let children = if let Some(node) = tree.get_node(id) {
            node.children.clone()
        } else {
            continue;
        };

        if children.is_empty() {
            continue;
        }

        // Sort children
        // To avoid double borrow (tree immutable in sort_by closure vs tree mutable in get_node_mut),
        // we can collect sort keys first.
        let mut child_keys: HashMap<NodeId, String> = HashMap::new();
        for &child_id in &children {
            child_keys.insert(child_id, get_sort_key(tree, &name_map, child_id));
        }

        if let Some(node) = tree.get_node_mut(id) {
            node.children.sort_by(|&a, &b| {
                // Get name or derive from children
                let name_a = child_keys.get(&a).map(|s| s.as_str()).unwrap_or("");
                let name_b = child_keys.get(&b).map(|s| s.as_str()).unwrap_or("");

                if descending {
                    name_b.cmp(name_a)
                } else {
                    name_a.cmp(name_b)
                }
            });
        }

        // After sorting children, update name_map for THIS node if it has no name
        // This ensures parents can use this node's derived name
        if let Some(node) = tree.get_node(id) {
            if node.name.as_deref().unwrap_or("").is_empty() {
                // Derived name is the name/key of the first child (after sort)
                if let Some(&first_child) = node.children.first() {
                    let child_key = get_sort_key(tree, &name_map, first_child);
                    name_map.insert(id, child_key);
                }
            }
        }
    }
}

fn get_sort_key(tree: &Tree, name_map: &HashMap<NodeId, String>, id: NodeId) -> String {
    if let Some(name) = name_map.get(&id) {
        if !name.is_empty() {
            return name.clone();
        }
    }
    // Fallback if not in map or empty (should be in map if visited post-order)
    // But if it was a leaf with empty name, it stays empty?
    // If it's a leaf, name_map has ""
    // If it's an internal node, we might have updated name_map in the loop.

    // If we are here, it means name_map has "" or missing.
    // If missing, try to get from tree
    if let Some(node) = tree.get_node(id) {
        return node.name.clone().unwrap_or_default();
    }

    String::new()
}

/// Sort the children of each node by the number of descendants (also known as ladderize).
///
/// # Arguments
///
/// * `tree` - The tree to modify.
/// * `descending` - If true, nodes with more descendants come first.
///
/// # Example
/// ```
/// use pgr::libs::phylo::tree::Tree;
/// use pgr::libs::phylo::algo;
///
/// // ((A,B),C)
/// // (A,B) has 2 descendants (leaves), C has 1 descendant.
/// let mut tree = Tree::from_newick("((A,B),C);").unwrap();
///
/// // Ascending: C (1) < (A,B) (2)
/// algo::sort_by_descendants(&mut tree, false);
/// assert_eq!(tree.to_newick(), "(C,(A,B));");
/// ```
pub fn sort_by_descendants(tree: &mut Tree, descending: bool) {
    if tree.is_empty() {
        return;
    }

    let root = tree.get_root().unwrap();
    let ids = match tree.levelorder(&root) {
        Ok(v) => v,
        Err(_) => return,
    };

    // Calculate descendant counts
    // Since `get_subtree` (preorder) returns the node itself + descendants,
    // the count is simply len() - 1 (if we strictly mean descendants) or len() (subtree size).
    // Subtree size is stable.
    let mut size_map: HashMap<NodeId, usize> = HashMap::new();

    // Optimization: Calculate sizes bottom-up (postorder) instead of calling get_subtree for each node.
    // But get_subtree is O(N) per call, making this O(N^2).
    // Let's use postorder to do it in O(N).
    if let Ok(post_ids) = tree.postorder(&root) {
        for &id in &post_ids {
            let mut count = 0;
            if let Some(node) = tree.get_node(id) {
                if node.is_leaf() {
                    count = 1; // Count self as 1 unit of "size"
                } else {
                    count = 1; // Self
                    for child in &node.children {
                        count += size_map.get(child).unwrap_or(&0);
                    }
                }
            }
            size_map.insert(id, count);
        }
    }

    // Now sort
    for id in ids {
        if let Some(node) = tree.get_node_mut(id) {
            if node.children.is_empty() {
                continue;
            }

            node.children.sort_by(|a, b| {
                let size_a = size_map.get(a).unwrap_or(&0);
                let size_b = size_map.get(b).unwrap_or(&0);
                if descending {
                    size_b.cmp(size_a)
                } else {
                    size_a.cmp(size_b)
                }
            });
        }
    }
}

/// Sort the children of each node based on a list of names.
///
/// Nodes are ordered by the position of their descendants in the provided list.
/// If a node has multiple descendants in the list, the minimum position is used.
/// Nodes with no descendants in the list are placed at the end.
///
/// # Arguments
///
/// * `tree` - The tree to modify.
/// * `order_list` - A list of names defining the desired order.
///
/// # Example
/// ```
/// use pgr::libs::phylo::tree::Tree;
/// use pgr::libs::phylo::algo;
///
/// let mut tree = Tree::from_newick("(A,B,C);").unwrap();
/// let order = vec!["C".to_string(), "B".to_string(), "A".to_string()];
/// algo::sort_by_list(&mut tree, &order);
/// assert_eq!(tree.to_newick(), "(C,B,A);");
/// ```
pub fn sort_by_list(tree: &mut Tree, order_list: &[String]) {
    if tree.is_empty() {
        return;
    }

    let root = tree.get_root().unwrap();

    // Map name -> position
    let mut pos_map: HashMap<String, usize> = HashMap::new();
    for (i, name) in order_list.iter().enumerate() {
        pos_map.insert(name.clone(), i);
    }

    let max_pos = order_list.len();
    let mut node_pos: HashMap<NodeId, usize> = HashMap::new();

    // Compute positions bottom-up (postorder)
    if let Ok(ids) = tree.postorder(&root) {
        for &id in &ids {
            let mut pos = max_pos;
            if let Some(node) = tree.get_node(id) {
                // 1. Check self name
                if let Some(name) = &node.name {
                    if let Some(&p) = pos_map.get(name) {
                        pos = p;
                    }
                }

                // 2. Check children (if self didn't override, or maybe we want min of self+children?)
                // Usually for ladderize/sorting, if I am a leaf, I use my name.
                // If I am internal, I use my children's derived positions.
                // The newick.rs logic was: check descendants.
                // Here, since we are bottom-up, children are already processed.
                for &child in &node.children {
                    if let Some(&child_p) = node_pos.get(&child) {
                        if child_p < pos {
                            pos = child_p;
                        }
                    }
                }
            }
            node_pos.insert(id, pos);
        }
    }

    // Sort children
    // We can iterate all nodes or use levelorder.
    if let Ok(ids) = tree.levelorder(&root) {
        for id in ids {
            if let Some(node) = tree.get_node_mut(id) {
                if node.children.is_empty() {
                    continue;
                }
                node.children.sort_by(|a, b| {
                    let pos_a = node_pos.get(a).unwrap_or(&max_pos);
                    let pos_b = node_pos.get(b).unwrap_or(&max_pos);
                    if pos_a == pos_b {
                        a.cmp(b) // Stable tie-break by ID
                    } else {
                        pos_a.cmp(pos_b)
                    }
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::phylo::tree::Tree;

    #[test]
    fn test_sort_by_name() {
        let mut tree = Tree::new();
        let root = tree.add_node();
        tree.set_root(root);

        let c1 = tree.add_node();
        tree.get_node_mut(c1).unwrap().name = Some("C".to_string());
        let c2 = tree.add_node();
        tree.get_node_mut(c2).unwrap().name = Some("A".to_string());
        let c3 = tree.add_node();
        tree.get_node_mut(c3).unwrap().name = Some("B".to_string());

        tree.add_child(root, c1).unwrap();
        tree.add_child(root, c2).unwrap();
        tree.add_child(root, c3).unwrap();

        // Before sort: C, A, B
        sort_by_name(&mut tree, false);

        let children = &tree.get_node(root).unwrap().children;
        assert_eq!(children.len(), 3);
        assert_eq!(
            tree.get_node(children[0]).unwrap().name.as_deref(),
            Some("A")
        );
        assert_eq!(
            tree.get_node(children[1]).unwrap().name.as_deref(),
            Some("B")
        );
        assert_eq!(
            tree.get_node(children[2]).unwrap().name.as_deref(),
            Some("C")
        );
    }

    #[test]
    fn test_sort_by_descendants() {
        let mut tree = Tree::new();
        let root = tree.add_node();
        tree.set_root(root);

        // Child 1: Leaf (Size 1)
        let c1 = tree.add_node();

        // Child 2: Has 2 children (Size 3)
        let c2 = tree.add_node();
        let c2_1 = tree.add_node();
        let c2_2 = tree.add_node();
        tree.add_child(c2, c2_1).unwrap();
        tree.add_child(c2, c2_2).unwrap();

        tree.add_child(root, c1).unwrap();
        tree.add_child(root, c2).unwrap();

        // Sort ascending (smallest first) -> c1, c2
        sort_by_descendants(&mut tree, false);
        let children = &tree.get_node(root).unwrap().children;
        assert_eq!(children[0], c1);
        assert_eq!(children[1], c2);

        // Sort descending (largest first) -> c2, c1
        sort_by_descendants(&mut tree, true);
        let children = &tree.get_node(root).unwrap().children;
        assert_eq!(children[0], c2);
        assert_eq!(children[1], c1);
    }

    #[test]
    fn test_sort_by_list() {
        let mut tree = Tree::new();
        let root = tree.add_node();
        tree.set_root(root);

        // ( (A, B), C )
        // Let's create a structure:
        // root -> n1 (children A, B)
        // root -> C
        let n1 = tree.add_node();
        let c = tree.add_node();
        tree.get_node_mut(c).unwrap().name = Some("C".to_string());

        tree.add_child(root, n1).unwrap();
        tree.add_child(root, c).unwrap();

        let a = tree.add_node();
        tree.get_node_mut(a).unwrap().name = Some("A".to_string());
        let b = tree.add_node();
        tree.get_node_mut(b).unwrap().name = Some("B".to_string());

        tree.add_child(n1, a).unwrap();
        tree.add_child(n1, b).unwrap();

        // Current order of root children: n1, c
        // Current order of n1 children: a, b

        // Target list: ["C", "B", "A"]
        // Expected:
        // root children: C (pos 0), n1 (pos min(pos(B)=1, pos(A)=2) = 1) -> (C, n1)
        // n1 children: B (pos 1), A (pos 2) -> (B, A)

        let order = vec!["C".to_string(), "B".to_string(), "A".to_string()];
        sort_by_list(&mut tree, &order);

        let root_children = &tree.get_node(root).unwrap().children;
        assert_eq!(root_children[0], c);
        assert_eq!(root_children[1], n1);

        let n1_children = &tree.get_node(n1).unwrap().children;
        assert_eq!(n1_children[0], b);
        assert_eq!(n1_children[1], a);
    }

    #[test]
    fn test_sort_by_list_comprehensive() {
        // Case 1: Simple case with only leaf nodes
        let newick = "(A,B,C);";
        let mut tree = Tree::from_newick(newick).unwrap();
        sort_by_list(
            &mut tree,
            &["C".to_string(), "B".to_string(), "A".to_string()],
        );
        assert_eq!(tree.to_newick(), "(C,B,A);");

        // Case 2: Case with internal nodes
        let newick = "((A,B),(C,D));";
        let mut tree = Tree::from_newick(newick).unwrap();
        sort_by_list(
            &mut tree,
            &["C".to_string(), "B".to_string(), "A".to_string()],
        );
        assert_eq!(tree.to_newick(), "((C,D),(B,A));");

        // Case 3: Case with internal nodes and names
        let newick = "((A,B)X,(C,D)Y);";
        let mut tree = Tree::from_newick(newick).unwrap();
        sort_by_list(
            &mut tree,
            &["C".to_string(), "B".to_string(), "A".to_string()],
        );
        assert_eq!(tree.to_newick(), "((C,D)Y,(B,A)X);");

        // Case 4: Case with unlisted nodes
        let newick = "((A,B),(C,E));";
        let mut tree = Tree::from_newick(newick).unwrap();
        sort_by_list(&mut tree, &["C".to_string(), "B".to_string()]);
        assert_eq!(tree.to_newick(), "((C,E),(B,A));");
    }
}
