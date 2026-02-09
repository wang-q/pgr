use super::*;
// use crate::libs::phylo::tree::Tree;

#[test]
fn test_tree_traversals() {
    let mut tree = Tree::new();
    //    0
    //   / \
    //  1   2
    // / \   \
    //3   4   5
    let n0 = tree.add_node();
    let n1 = tree.add_node();
    let n2 = tree.add_node();
    let n3 = tree.add_node();
    let n4 = tree.add_node();
    let n5 = tree.add_node();

    tree.set_root(n0);
    tree.add_child(n0, n1).unwrap();
    tree.add_child(n0, n2).unwrap();
    tree.add_child(n1, n3).unwrap();
    tree.add_child(n1, n4).unwrap();
    tree.add_child(n2, n5).unwrap();

    // Preorder: 0, 1, 3, 4, 2, 5
    let pre = tree.preorder(&n0).unwrap();
    assert_eq!(pre, vec![n0, n1, n3, n4, n2, n5]);

    // Postorder: 3, 4, 1, 5, 2, 0
    let post = tree.postorder(&n0).unwrap();
    assert_eq!(post, vec![n3, n4, n1, n5, n2, n0]);

    // Levelorder: 0, 1, 2, 3, 4, 5
    let level = tree.levelorder(&n0).unwrap();
    assert_eq!(level, vec![n0, n1, n2, n3, n4, n5]);
}

#[test]
fn test_tree_basic_ops() {
    let mut tree = Tree::new();

    // Create nodes
    // 0(root) -> 1, 2
    // 1 -> 3
    let n0 = tree.add_node();
    let n1 = tree.add_node();
    let n2 = tree.add_node();
    let n3 = tree.add_node();

    tree.set_root(n0);

    assert_eq!(tree.add_child(n0, n1), Ok(()));
    assert_eq!(tree.add_child(n0, n2), Ok(()));
    assert_eq!(tree.add_child(n1, n3), Ok(()));

    assert_eq!(tree.len(), 4);

    // Check structure
    let root = tree.get_node(n0).unwrap();
    assert_eq!(root.children, vec![n1, n2]);

    let node1 = tree.get_node(n1).unwrap();
    assert_eq!(node1.parent, Some(n0));
    assert_eq!(node1.children, vec![n3]);
}

#[test]
fn test_tree_remove_and_compact() {
    let mut tree = Tree::new();
    // 0 -> 1 -> 2
    let n0 = tree.add_node();
    let n1 = tree.add_node();
    let n2 = tree.add_node();

    tree.add_child(n0, n1).unwrap();
    tree.add_child(n1, n2).unwrap();
    tree.set_root(n0);

    // Remove n1 (recursive=false), n2 becomes orphan
    tree.remove_node(n1, false);

    assert!(tree.get_node(n1).is_none()); // n1 is logically gone
    assert_eq!(tree.len(), 2); // 0 and 2 remain

    let node0 = tree.get_node(n0).unwrap();
    assert!(!node0.children.contains(&n1)); // 0 no longer points to 1

    let node2 = tree.get_node(n2).unwrap();
    assert_eq!(node2.parent, None); // 2 is orphaned

    // Compact
    // Before: [0:Valid, 1:Deleted, 2:Valid]
    // After:  [0':Old0, 1':Old2]
    tree.compact();

    assert_eq!(tree.len(), 2);
    // Old n0 should be at index 0
    let new_n0 = tree.get_node(0).unwrap();
    assert_eq!(new_n0.children.len(), 0);

    // Old n2 should be at index 1
    let new_n1 = tree.get_node(1).unwrap();
    assert_eq!(new_n1.parent, None);
}

#[test]
fn test_tree_paths_and_distances() {
    let mut tree = Tree::new();
    //    0
    //   / \
    //  1   2
    // / \
    //3   4
    let n0 = tree.add_node();
    let n1 = tree.add_node();
    let n2 = tree.add_node();
    let n3 = tree.add_node();
    let n4 = tree.add_node();

    tree.set_root(n0);
    tree.add_child(n0, n1).unwrap();
    tree.add_child(n0, n2).unwrap();
    tree.add_child(n1, n3).unwrap();
    tree.add_child(n1, n4).unwrap();

    // Set lengths
    tree.get_node_mut(n1).unwrap().length = Some(1.0);
    tree.get_node_mut(n2).unwrap().length = Some(2.0);
    tree.get_node_mut(n3).unwrap().length = Some(3.0);
    tree.get_node_mut(n4).unwrap().length = Some(4.0);

    // Paths
    assert_eq!(tree.get_path_from_root(&n3).unwrap(), vec![n0, n1, n3]);
    assert_eq!(tree.get_path_from_root(&n2).unwrap(), vec![n0, n2]);

    // LCA
    assert_eq!(tree.get_common_ancestor(&n3, &n4).unwrap(), n1);
    assert_eq!(tree.get_common_ancestor(&n3, &n2).unwrap(), n0);
    assert_eq!(tree.get_common_ancestor(&n1, &n3).unwrap(), n1);

    // Distance
    // n3 -> n4: n3(3.0)->n1 + n1->n4(4.0) = 7.0 (weighted). Steps: n3->n1->n4 = 2 edges.
    let (w, t) = tree.get_distance(&n3, &n4).unwrap();
    assert_eq!(w, 7.0);
    assert_eq!(t, 2);

    // n3 -> n2: n3(3.0)->n1(1.0)->n0 + n0->n2(2.0) = 6.0. Steps: n3->n1->n0->n2 = 3 edges.
    let (w, t) = tree.get_distance(&n3, &n2).unwrap();
    assert_eq!(w, 6.0);
    assert_eq!(t, 3);
}

#[test]
fn test_tree_queries() {
    let mut tree = Tree::new();
    //    0 (root, name="root")
    //   / \
    //  1   2 (leaf, name="leaf2")
    // /
    //3 (leaf, name="leaf3")
    let n0 = tree.add_node();
    tree.get_node_mut(n0).unwrap().name = Some("root".to_string());

    let n1 = tree.add_node();

    let n2 = tree.add_node();
    tree.get_node_mut(n2).unwrap().name = Some("leaf2".to_string());

    let n3 = tree.add_node();
    tree.get_node_mut(n3).unwrap().name = Some("leaf3".to_string());

    tree.set_root(n0);
    tree.add_child(n0, n1).unwrap();
    tree.add_child(n0, n2).unwrap();
    tree.add_child(n1, n3).unwrap();

    // Subtree
    // subtree(1) = [1, 3]
    let sub = tree.extract_subtree(&n1).unwrap();
    // IDs might change, but names should persist?
    // Subtree renumbers. Root of subtree is 0. Child is 1.
    // Names: None, leaf3.
    let sub_root = sub.get_root().unwrap();
    assert_eq!(sub.get_node(sub_root).unwrap().name, None); // Old n1
    let children = &sub.get_node(sub_root).unwrap().children;
    assert_eq!(children.len(), 1);
    assert_eq!(
        sub.get_node(children[0]).unwrap().name,
        Some("leaf3".to_string())
    );

    // Leaves
    // Leaves: 2, 3
    let leaves = tree.get_leaves();
    assert!(leaves.contains(&n2));
    assert!(leaves.contains(&n3));
    assert_eq!(leaves.len(), 2);

    // Find nodes
    let named_nodes = tree.find_nodes(|n| n.name.is_some());
    assert_eq!(named_nodes.len(), 3); // 0, 2, 3

    // Get by name
    assert_eq!(tree.get_node_by_name("root"), Some(n0));
    assert_eq!(tree.get_node_by_name("leaf2"), Some(n2));
    assert_eq!(tree.get_node_by_name("leaf3"), Some(n3));
    assert_eq!(tree.get_node_by_name("nonexistent"), None);
}

#[test]
fn test_tree_prune() {
    let mut tree = Tree::new();
    //    0
    //   / \
    //  1   2
    // /
    //3
    let n0 = tree.add_node();
    tree.get_node_mut(n0).unwrap().name = Some("root".to_string());

    let n1 = tree.add_node();
    tree.get_node_mut(n1).unwrap().name = Some("n1".to_string());

    let n2 = tree.add_node();
    tree.get_node_mut(n2).unwrap().name = Some("remove_me".to_string());

    let n3 = tree.add_node();
    tree.get_node_mut(n3).unwrap().name = Some("n3".to_string());

    tree.set_root(n0);
    tree.add_child(n0, n1).unwrap();
    tree.add_child(n0, n2).unwrap();
    tree.add_child(n1, n3).unwrap();

    assert_eq!(tree.len(), 4);

    // Prune node with name "remove_me"
    tree.prune_where(|n| n.name.as_deref() == Some("remove_me"));

    assert_eq!(tree.len(), 3);
    assert!(tree.get_node(n2).is_none());
    assert!(tree.get_node(n0).unwrap().children.contains(&n1));
    assert!(!tree.get_node(n0).unwrap().children.contains(&n2));

    // Prune n1, should also remove n3
    tree.prune_where(|n| n.id == n1);

    assert_eq!(tree.len(), 1); // Only root left
    assert!(tree.get_node(n1).is_none());
    assert!(tree.get_node(n3).is_none());
    assert!(tree.get_node(n0).unwrap().children.is_empty());
}

#[test]
fn test_get_partitions() {
    // Tree: ((A,B)C,D)E;
    let newick = "((A,B)C,D)E;";
    let tree = Tree::from_newick(newick).unwrap();

    let partitions = tree.get_partitions();

    // Helper to create string sets
    fn make_set(names: &[&str]) -> std::collections::BTreeSet<String> {
        names.iter().map(|s| s.to_string()).collect()
    }

    let mut expected = std::collections::HashSet::new();
    expected.insert(make_set(&["A"]));
    expected.insert(make_set(&["B"]));
    expected.insert(make_set(&["A", "B"]));
    expected.insert(make_set(&["D"]));
    expected.insert(make_set(&["A", "B", "D"]));

    assert_eq!(partitions.len(), 5);
    for p in &expected {
        assert!(partitions.contains(p), "Missing partition: {:?}", p);
    }
}

#[test]
fn test_is_binary() {
    let t1 = Tree::from_newick("((A,B),C);").unwrap();
    assert!(t1.is_binary());

    let t2 = Tree::from_newick("(A,B,C);").unwrap();
    assert!(!t2.is_binary());
}

#[test]
fn test_get_leaves() {
    let tree = Tree::from_newick("((A,B)C,D)E;").unwrap();
    let leaves = tree.get_leaf_names();
    let leaf_names: Vec<String> = leaves.into_iter().map(|n| n.unwrap()).collect();

    assert!(leaf_names.contains(&"A".to_string()));
    assert!(leaf_names.contains(&"B".to_string()));
    assert!(leaf_names.contains(&"D".to_string()));
    assert_eq!(leaf_names.len(), 3);
}

#[test]
fn test_diameter() {
    let newick = "((A:1,B:2):1,C:4);";
    let tree = Tree::from_newick(newick).unwrap();
    // Dist(A,B) = 3
    // Dist(A,C) = 1+1+4 = 6
    // Dist(B,C) = 2+1+4 = 7
    assert_eq!(tree.diameter().unwrap(), 7.0);
}

#[test]
fn test_robinson_foulds() {
    let t1 = Tree::from_newick("((A,B),C);").unwrap();
    let t2 = Tree::from_newick("((A,C),B);").unwrap();
    assert_eq!(t1.robinson_foulds(&t2).unwrap(), 2);

    let t3 = Tree::from_newick("((A,B),C);").unwrap();
    assert_eq!(t1.robinson_foulds(&t3).unwrap(), 0);
}

#[test]
fn test_deroot() {
    let mut tree = Tree::from_newick("((A:1,B:2)C:3,D:4)Root;").unwrap();
    tree.deroot().unwrap();

    let root = tree.get_root().unwrap();
    let children = &tree.get_node(root).unwrap().children;
    assert_eq!(children.len(), 3);

    // Check names of children
    let child_names: Vec<String> = children
        .iter()
        .map(|&id| tree.get_node(id).unwrap().name.clone().unwrap_or_default())
        .collect();

    assert!(child_names.contains(&"A".to_string()));
    assert!(child_names.contains(&"B".to_string()));
    assert!(child_names.contains(&"D".to_string()));
}

#[test]
fn test_reroot_support_values() {
    let mut tree = Tree::from_newick("(A,(B,C)Support)Root;").unwrap();
    let c_id = tree.get_node_by_name("C").unwrap();

    tree.reroot_at(c_id, true).unwrap();

    // C is root
    assert_eq!(tree.get_root(), Some(c_id));

    let root = tree.get_node(tree.get_root().unwrap()).unwrap();
    assert_eq!(root.name.as_deref(), Some("C"));

    let support_node_id = root.children[0]; // The old Support node
    let support_node = tree.get_node(support_node_id).unwrap();
    assert_eq!(support_node.name, None);

    let old_root_node_id = support_node
        .children
        .iter()
        .find(|&&id| {
            // Find the one that has A as child
            let n = tree.get_node(id).unwrap();
            n.children
                .iter()
                .any(|&child| tree.get_node(child).unwrap().name.as_deref() == Some("A"))
        })
        .unwrap();

    let old_root_node = tree.get_node(*old_root_node_id).unwrap();
    assert_eq!(old_root_node.name.as_deref(), Some("Support"));
}

#[test]
fn test_reroot_longest_branch() {
    // (A:1, B:2)Root;
    // Longest branch is B (len 2).
    // Reroot should pick B.
    let mut tree = Tree::from_newick("(A:1,B:2)Root;").unwrap();

    let target = tree.get_node_with_longest_edge().unwrap();
    let b_id = tree.get_node_by_name("B").unwrap();
    assert_eq!(target, b_id);

    tree.reroot_at(target, false).unwrap();
    assert_eq!(tree.get_root(), Some(b_id));
}
