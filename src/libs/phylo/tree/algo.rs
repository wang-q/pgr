use super::Tree;
use crate::libs::phylo::node::NodeId;
use std::collections::{HashMap, HashSet};

/// Sort the children of each node by their name (label).
///
/// # Arguments
///
/// * `tree` - The tree to modify.
/// * `descending` - If true, sort in descending order (Z-A).
///
/// # Example
/// ```ignore
/// use pgr::libs::phylo::tree::Tree;
/// use pgr::libs::phylo::tree::algo;
///
/// let mut tree = Tree::from_newick("(B,A);").unwrap();
/// algo::sort_by_name(&mut tree, false);
/// assert_eq!(tree.to_newick(), "(A,B);");
/// ```
pub fn sort_by_name(tree: &mut Tree, descending: bool) {
    if tree.is_empty() {
        return;
    }

    let root = tree
        .get_root()
        .expect("internal: non-empty tree has a root");
    let ids = tree.postorder(&root);

    // Pre-collect names to avoid borrowing issues during sort
    let mut name_map: HashMap<NodeId, String> = HashMap::new();
    for &id in &ids {
        if let Some(node) = tree.get_node(id) {
            name_map.insert(id, node.name.clone().unwrap_or_default());
        }
    }

    for id in ids {
        let children = if let Some(node) = tree.get_node(id) {
            node.children.clone()
        } else {
            continue;
        };

        if children.is_empty() {
            continue;
        }

        let mut child_keys: HashMap<NodeId, String> = HashMap::new();
        for &child_id in &children {
            child_keys.insert(child_id, get_sort_key(&name_map, child_id));
        }

        if let Some(node) = tree.get_node_mut(id) {
            node.children.sort_by(|&a, &b| {
                let name_a = child_keys.get(&a).map(|s| s.as_str()).unwrap_or("");
                let name_b = child_keys.get(&b).map(|s| s.as_str()).unwrap_or("");

                if descending {
                    name_b.cmp(name_a)
                } else {
                    name_a.cmp(name_b)
                }
            });
        }

        if let Some(node) = tree.get_node(id) {
            if node.name.as_deref().unwrap_or("").is_empty() {
                if let Some(&first_child) = node.children.first() {
                    let child_key = get_sort_key(&name_map, first_child);
                    name_map.insert(id, child_key);
                }
            }
        }
    }
}

fn get_sort_key(name_map: &HashMap<NodeId, String>, id: NodeId) -> String {
    name_map.get(&id).cloned().unwrap_or_default()
}

/// Sort the children of each node by the number of descendants (also known as ladderize).
///
/// # Arguments
///
/// * `tree` - The tree to modify.
/// * `descending` - If true, nodes with more descendants come first.
///
/// # Example
/// ```ignore
/// use pgr::libs::phylo::tree::Tree;
/// use pgr::libs::phylo::tree::algo;
///
/// // ((A,B),C)
/// // (A,B) has 2 descendants (leaves), C has 1 descendant.
/// let mut tree = Tree::from_newick("((A,B),C);").unwrap();
///
/// // Ascending: C (1) < (A,B) (2)
/// algo::ladderize(&mut tree, false);
/// assert_eq!(tree.to_newick(), "(C,(A,B));");
/// ```
pub fn ladderize(tree: &mut Tree, descending: bool) {
    if tree.is_empty() {
        return;
    }

    let root = tree
        .get_root()
        .expect("internal: non-empty tree has a root");
    let ids = tree.levelorder(&root);

    let mut size_map: HashMap<NodeId, usize> = HashMap::new();

    let post_ids = tree.postorder(&root);
    for &id in &post_ids {
        let mut count = 0;
        if let Some(node) = tree.get_node(id) {
            if node.is_leaf() {
                count = 1;
            } else {
                count = 1;
                for child in &node.children {
                    count += size_map.get(child).unwrap_or(&0);
                }
            }
        }
        size_map.insert(id, count);
    }

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
/// ```ignore
/// use pgr::libs::phylo::tree::Tree;
/// use pgr::libs::phylo::tree::algo;
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

    let root = tree
        .get_root()
        .expect("internal: non-empty tree has a root");

    // Map name -> position
    let mut pos_map: HashMap<String, usize> = HashMap::new();
    for (i, name) in order_list.iter().enumerate() {
        pos_map.insert(name.clone(), i);
    }

    let max_pos = order_list.len();
    let mut node_pos: HashMap<NodeId, usize> = HashMap::new();

    let ids = tree.postorder(&root);
    for &id in &ids {
        let mut pos = max_pos;
        if let Some(node) = tree.get_node(id) {
            if let Some(name) = &node.name {
                if let Some(&p) = pos_map.get(name) {
                    pos = p;
                }
            }

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

    let ids = tree.levelorder(&root);
    for id in ids {
        if let Some(node) = tree.get_node_mut(id) {
            if node.children.is_empty() {
                continue;
            }
            node.children.sort_by(|a, b| {
                let pos_a = node_pos.get(a).unwrap_or(&max_pos);
                let pos_b = node_pos.get(b).unwrap_or(&max_pos);
                if pos_a == pos_b {
                    a.cmp(b)
                } else {
                    pos_a.cmp(pos_b)
                }
            });
        }
    }
}

/// Sort the children of each node by the number of descendants, alternating direction at each level.
///
/// This produces a "balanced" look for the tree.
/// Level 0 (Root children): Ascending (Light -> Heavy)
/// Level 1: Descending (Heavy -> Light)
/// ...
pub fn deladderize(tree: &mut Tree) {
    if tree.is_empty() {
        return;
    }

    let root = tree
        .get_root()
        .expect("internal: non-empty tree has a root");

    // 1. Calculate descendant counts (same as ladderize)
    let mut size_map: HashMap<NodeId, usize> = HashMap::new();
    let post_ids = tree.postorder(&root);
    for &id in &post_ids {
        let mut count = 0;
        if let Some(node) = tree.get_node(id) {
            if node.is_leaf() {
                count = 1;
            } else {
                count = 1;
                for child in &node.children {
                    count += size_map.get(child).unwrap_or(&0);
                }
            }
        }
        size_map.insert(id, count);
    }

    // 2. Traversal with state
    let mut queue = std::collections::VecDeque::new();
    queue.push_back((root, false)); // Start ascending

    while let Some((id, descending)) = queue.pop_front() {
        let children_ids = if let Some(node) = tree.get_node_mut(id) {
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
            node.children.clone()
        } else {
            continue;
        };

        for child in children_ids {
            queue.push_back((child, !descending));
        }
    }
}

/// Compute the set of node IDs to keep when inverting a prune around `targets`.
///
/// Returns the union of `targets`, all descendants of any target, and all
/// ancestors of any target. Nodes not in the returned set should be removed.
pub fn compute_keep_set<I>(tree: &Tree, targets: I) -> HashSet<NodeId>
where
    I: IntoIterator<Item = NodeId>,
{
    let mut keep = HashSet::new();
    let Some(root) = tree.get_root() else {
        return keep;
    };

    let target_set: HashSet<NodeId> = targets.into_iter().collect();
    let mut is_in_clade = HashSet::new();

    // Pass 1: Downward propagation (descendants of targets).
    // levelorder visits parents before children, so `is_in_clade` propagates.
    let all_nodes = tree.levelorder(&root);
    for &id in &all_nodes {
        let mut kept = target_set.contains(&id);
        if !kept {
            if let Some(node) = tree.get_node(id) {
                if let Some(parent) = node.parent {
                    if is_in_clade.contains(&parent) {
                        kept = true;
                    }
                }
            }
        }
        if kept {
            is_in_clade.insert(id);
            keep.insert(id);
        }
    }

    // Pass 2: Upward propagation (ancestors of kept nodes).
    // Iterate in reverse so children are visited before their parents.
    for &id in all_nodes.iter().rev() {
        if keep.contains(&id) {
            if let Some(node) = tree.get_node(id) {
                if let Some(parent) = node.parent {
                    keep.insert(parent);
                }
            }
        }
    }

    keep
}

/// Remove `to_remove` nodes from `tree`, then clean up internal nodes that
/// became leaves and collapse degree-2 nodes (single-child internals).
///
/// `to_remove` may include both leaves and internal nodes; removal is recursive
/// (subtrees are detached with their parent).
pub fn prune_nodes(tree: &mut Tree, to_remove: Vec<NodeId>) {
    // 1. Snapshot internal nodes before pruning, so we can detect those that
    //    become leaves after removal.
    let mut old_internals = Vec::new();
    if let Some(root) = tree.get_root() {
        let all_nodes = tree.levelorder(&root);
        for id in all_nodes {
            if let Some(node) = tree.get_node(id) {
                if !node.children.is_empty() {
                    old_internals.push(id);
                }
            }
        }
    }

    // 2. Remove the requested nodes.
    for id in to_remove {
        tree.remove_node(id, true);
    }

    // 3. Clean up internals that became leaves (reverse so deeper nodes first).
    for id in old_internals.into_iter().rev() {
        if let Some(node) = tree.get_node(id) {
            if node.children.is_empty() {
                tree.remove_node(id, true);
            }
        }
    }

    // 4. Collapse degree-2 nodes (post-order so children are visited first).
    if let Some(root) = tree.get_root() {
        let nodes = tree.postorder(&root);
        for id in nodes {
            if let Some(node) = tree.get_node(id) {
                if node.children.len() == 1 {
                    if tree.get_root() == Some(id) {
                        // Root with a single child: promote the child to root.
                        let child_id = node.children[0];
                        tree.set_root(child_id);
                        tree.remove_node(id, false);
                    } else {
                        let _ = tree.collapse_node(id);
                    }
                }
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
    fn test_ladderize() {
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
        ladderize(&mut tree, false);
        let children = &tree.get_node(root).unwrap().children;
        assert_eq!(children[0], c1);
        assert_eq!(children[1], c2);

        // Sort descending (largest first) -> c2, c1
        ladderize(&mut tree, true);
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

    #[test]
    fn test_deladderize() {
        let mut tree = Tree::new();
        let root = tree.add_node();
        tree.set_root(root);

        // Structure: ((A,B),(C,(D,E)),F)
        // Sizes:
        // A,B,C,D,E,F (leaves) = 1
        // (A,B) = 1 + 1 + 1 = 3
        // (D,E) = 1 + 1 + 1 = 3
        // (C,(D,E)) = 1 + 1 + 3 = 5
        // Root children:
        // 1. (A,B) - size 3
        // 2. (C,(D,E)) - size 5
        // 3. F - size 1

        let f = tree.add_node();
        tree.get_node_mut(f).unwrap().name = Some("F".to_string());

        let ab = tree.add_node();
        let a = tree.add_node();
        tree.get_node_mut(a).unwrap().name = Some("A".to_string());
        let b = tree.add_node();
        tree.get_node_mut(b).unwrap().name = Some("B".to_string());
        tree.add_child(ab, a).unwrap();
        tree.add_child(ab, b).unwrap();

        let cde = tree.add_node();
        let c = tree.add_node();
        tree.get_node_mut(c).unwrap().name = Some("C".to_string());
        let de = tree.add_node();
        let d = tree.add_node();
        tree.get_node_mut(d).unwrap().name = Some("D".to_string());
        let e = tree.add_node();
        tree.get_node_mut(e).unwrap().name = Some("E".to_string());
        tree.add_child(de, d).unwrap();
        tree.add_child(de, e).unwrap();

        tree.add_child(cde, c).unwrap(); // Add C first
        tree.add_child(cde, de).unwrap(); // Add (D,E) second

        tree.add_child(root, ab).unwrap();
        tree.add_child(root, cde).unwrap();
        tree.add_child(root, f).unwrap();

        // Run deladderize
        // Level 0 (Root children): Ascending -> F (1), (A,B) (3), (C,(D,E)) (5)
        // Level 1 (Children of Level 0 nodes): Descending
        // - F: no children
        // - (A,B): A(1), B(1) -> Equal, keep order (A,B)
        // - (C,(D,E)): C(1), (D,E)(3) -> Descending -> (D,E), C
        // Level 2 (Children of Level 1 nodes): Ascending
        // - (D,E): D(1), E(1) -> Equal, keep order (D,E)

        deladderize(&mut tree);

        let root_children = &tree.get_node(root).unwrap().children;
        assert_eq!(root_children.len(), 3);
        assert_eq!(root_children[0], f);
        assert_eq!(root_children[1], ab);
        assert_eq!(root_children[2], cde);

        let cde_children = &tree.get_node(cde).unwrap().children;
        assert_eq!(cde_children[0], de); // Larger one first (descending)
        assert_eq!(cde_children[1], c);
    }
}
