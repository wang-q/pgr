use super::node::{Node, NodeId};
use super::writer;
use std::collections::{BTreeMap, HashMap};

#[derive(Debug, Default, Clone)]
pub struct Tree {
    /// Arena storage for all nodes
    nodes: Vec<Node>,

    /// Optional root ID (a tree might be empty or in construction)
    root: Option<NodeId>,
}

impl Tree {
    /// Create a new empty tree
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let tree = Tree::new();
    /// assert!(tree.is_empty());
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a new node to the tree. Returns the new node's ID.
    /// The node is initially detached (no parent).
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// let id = tree.add_node();
    /// assert_eq!(tree.len(), 1);
    /// ```
    pub fn add_node(&mut self) -> NodeId {
        let id = self.nodes.len();
        let node = Node::new(id);
        self.nodes.push(node);
        id
    }

    /// Get a reference to a node by ID. Returns None if ID is invalid or node is deleted.
    pub fn get_node(&self, id: NodeId) -> Option<&Node> {
        self.nodes.get(id).filter(|n| !n.deleted)
    }

    /// Get a mutable reference to a node by ID.
    pub fn get_node_mut(&mut self, id: NodeId) -> Option<&mut Node> {
        self.nodes.get_mut(id).filter(|n| !n.deleted)
    }

    /// Set a node as the root of the tree.
    pub fn set_root(&mut self, id: NodeId) {
        if self.get_node(id).is_some() {
            self.root = Some(id);
        }
    }

    /// Add a child to a parent node.
    /// Updates both parent's `children` list and child's `parent` field.
    ///
    /// # Errors
    /// Returns error if parent/child invalid, deleted, or cycle detected (basic check).
    pub fn add_child(&mut self, parent_id: NodeId, child_id: NodeId) -> Result<(), String> {
        // Validation
        if parent_id == child_id {
            return Err("Cannot add node as child of itself".to_string());
        }
        if self.get_node(parent_id).is_none() {
            return Err(format!("Parent node {} not found or deleted", parent_id));
        }
        if self.get_node(child_id).is_none() {
            return Err(format!("Child node {} not found or deleted", child_id));
        }

        // Check if child already has a parent
        let child_parent = self.nodes[child_id].parent;
        if let Some(old_parent) = child_parent {
            return Err(format!(
                "Node {} already has parent {}",
                child_id, old_parent
            ));
        }

        // Link
        self.nodes[child_id].parent = Some(parent_id);
        self.nodes[parent_id].children.push(child_id);

        Ok(())
    }

    /// Soft remove a node and its descendants (optional recursive).
    /// If recursive is false, children are orphaned (parent set to None).
    pub fn remove_node(&mut self, id: NodeId, recursive: bool) {
        if id >= self.nodes.len() || self.nodes[id].deleted {
            return;
        }

        // 1. Handle Parent Relation
        if let Some(parent_id) = self.nodes[id].parent {
            // Remove self from parent's children list
            // Note: This is O(N) for the children list, but usually small.
            if let Some(parent) = self.nodes.get_mut(parent_id) {
                parent.children.retain(|&child| child != id);
            }
        }

        // 2. Handle Children
        let children = self.nodes[id].children.clone();
        for child_id in children {
            if recursive {
                self.remove_node(child_id, true);
            } else {
                // Orphan the child
                if let Some(child) = self.nodes.get_mut(child_id) {
                    child.parent = None;
                }
            }
        }

        // 3. Mark as deleted
        if let Some(node) = self.nodes.get_mut(id) {
            node.deleted = true;
            node.children.clear();
            node.parent = None;
        }

        // 4. Update root if needed
        if self.root == Some(id) {
            self.root = None;
        }
    }

    /// Collapse a node, removing it and connecting its children to its parent.
    /// Edge lengths are summed (parent->node + node->child).
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// // 0 -> 1 -> 2
    /// let n0 = tree.add_node();
    /// let n1 = tree.add_node();
    /// let n2 = tree.add_node();
    /// tree.set_root(n0);
    ///
    /// tree.add_child(n0, n1).unwrap();
    /// tree.get_node_mut(n1).unwrap().length = Some(1.0);
    ///
    /// tree.add_child(n1, n2).unwrap();
    /// tree.get_node_mut(n2).unwrap().length = Some(2.0);
    ///
    /// // Collapse n1
    /// tree.collapse_node(n1).unwrap();
    ///
    /// // n0 -> n2 (len 3.0)
    /// let node2 = tree.get_node(n2).unwrap();
    /// assert_eq!(node2.parent, Some(n0));
    /// assert_eq!(node2.length, Some(3.0));
    /// assert!(tree.get_node(n0).unwrap().children.contains(&n2));
    /// assert!(!tree.get_node(n0).unwrap().children.contains(&n1));
    ///
    /// // From Newick
    /// let newick = "(A,(B)D);";
    /// let mut tree = Tree::from_newick(newick).unwrap();
    /// let id = *tree.get_name_id().get("D").unwrap();
    ///
    /// tree.collapse_node(id).unwrap();
    ///
    /// assert_eq!(tree.to_newick(), "(A,B);");
    ///
    /// // Case with edge lengths
    /// let newick = "(A:1,(B:1)D:2);";
    /// let mut tree = Tree::from_newick(newick).unwrap();
    /// let id = *tree.get_name_id().get("D").unwrap();
    ///
    /// tree.collapse_node(id).unwrap();
    ///
    /// assert_eq!(tree.to_newick(), "(A:1,B:3);");
    /// ```
    pub fn collapse_node(&mut self, id: NodeId) -> Result<(), String> {
        if self.get_node(id).is_none() {
            return Err(format!("Node {} not found", id));
        }
        if self.root == Some(id) {
            return Err("Cannot collapse root node".to_string());
        }

        // 1. Get info
        let (parent_id, parent_edge) = {
            let node = self.get_node(id).unwrap();
            // Safety: Checked root above, so parent must exist
            (node.parent.unwrap(), node.length)
        };

        let children_info: Vec<(NodeId, Option<f64>)> = {
            let node = self.get_node(id).unwrap();
            node.children
                .iter()
                .map(|&c| {
                    let child_node = self.nodes.get(c).unwrap();
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
            if let Some(child) = self.get_node_mut(child_id) {
                child.parent = Some(parent_id);
                child.length = new_edge;
            }
            new_children_ids.push(child_id);
        }

        // 3. Update parent's children list
        if let Some(parent) = self.get_node_mut(parent_id) {
            if let Some(pos) = parent.children.iter().position(|&x| x == id) {
                parent.children.splice(pos..pos + 1, new_children_ids);
            }
        }

        // 4. Mark deleted
        if let Some(node) = self.get_node_mut(id) {
            node.deleted = true;
            node.children.clear();
            node.parent = None;
        }

        Ok(())
    }

    /// Compact the tree by removing soft-deleted nodes and remapping IDs.
    /// This invalidates all existing NodeIds held outside!
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// let n0 = tree.add_node();
    /// let n1 = tree.add_node();
    /// tree.remove_node(n0, false);
    /// tree.compact();
    /// assert_eq!(tree.len(), 1);
    /// ```
    pub fn compact(&mut self) {
        let mut old_to_new = HashMap::new();
        let mut new_nodes = Vec::with_capacity(self.nodes.len());
        let mut new_idx = 0;

        // 1. Build mapping and new node list (without edges first)
        for old_node in &self.nodes {
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
        for (old_idx, node) in self.nodes.iter().enumerate() {
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
        if let Some(old_root) = self.root {
            self.root = old_to_new.get(&old_root).copied();
        }

        // 4. Swap
        self.nodes = new_nodes;
    }

    /// Perform a preorder traversal starting from a given node.
    /// Returns a vector of NodeIds in visitation order.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// let n0 = tree.add_node();
    /// let n1 = tree.add_node();
    /// let n2 = tree.add_node();
    /// tree.add_child(n0, n1);
    /// tree.add_child(n0, n2);
    /// let traversal = tree.preorder(&n0).unwrap();
    /// assert_eq!(traversal, vec![n0, n1, n2]);
    /// ```
    pub fn preorder(&self, start_node: &NodeId) -> Result<Vec<NodeId>, String> {
        if self.get_node(*start_node).is_none() {
            return Err(format!("Node {} not found", start_node));
        }

        let mut result = Vec::new();
        let mut stack = vec![*start_node];

        while let Some(curr) = stack.pop() {
            result.push(curr);
            // Push children in reverse order so they are popped in original order (left-to-right)
            if let Some(node) = self.get_node(curr) {
                for &child in node.children.iter().rev() {
                    stack.push(child);
                }
            }
        }
        Ok(result)
    }

    /// Get the path from the root to the specified node.
    /// Returns a vector of NodeIds starting from the root and ending at the target node.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// let n0 = tree.add_node();
    /// let n1 = tree.add_node();
    /// let n2 = tree.add_node();
    /// tree.add_child(n0, n1);
    /// tree.add_child(n1, n2);
    ///
    /// let path = tree.get_path_from_root(&n2).unwrap();
    /// assert_eq!(path, vec![n0, n1, n2]);
    ///
    /// // Error: Node not in tree (e.g., random large ID or deleted)
    /// assert!(tree.get_path_from_root(&9999).is_err());
    /// ```
    pub fn get_path_from_root(&self, target_node: &NodeId) -> Result<Vec<NodeId>, String> {
        if self.get_node(*target_node).is_none() {
            return Err(format!("Node {} not found", target_node));
        }

        let mut path = Vec::new();
        let mut curr = *target_node;

        // Traverse upwards
        loop {
            path.push(curr);
            if let Some(node) = self.get_node(curr) {
                if let Some(parent) = node.parent {
                    curr = parent;
                } else {
                    break;
                }
            } else {
                break;
            }
        }

        path.reverse();
        Ok(path)
    }

    /// Find the Lowest Common Ancestor (LCA) of two nodes.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// //    0
    /// //   / \
    /// //  1   2
    /// // / \
    /// //3   4
    /// let n0 = tree.add_node();
    /// let n1 = tree.add_node();
    /// let n2 = tree.add_node();
    /// let n3 = tree.add_node();
    /// let n4 = tree.add_node();
    /// tree.add_child(n0, n1);
    /// tree.add_child(n0, n2);
    /// tree.add_child(n1, n3);
    /// tree.add_child(n1, n4);
    ///
    /// assert_eq!(tree.get_common_ancestor(&n3, &n4).unwrap(), n1);
    /// assert_eq!(tree.get_common_ancestor(&n3, &n2).unwrap(), n0);
    ///
    /// // Error: Nodes in disjoint trees
    /// let mut tree2 = Tree::new();
    /// let m0 = tree2.add_node();
    /// // Since NodeIds are unique to the Tree instance's arena, we can't easily cross-reference
    /// // without merging. But if we simulate a disconnected graph within one Tree:
    /// let orphan = tree.add_node(); // Not connected to n0
    /// assert!(tree.get_common_ancestor(&n3, &orphan).is_err());
    /// ```
    pub fn get_common_ancestor(&self, a: &NodeId, b: &NodeId) -> Result<NodeId, String> {
        let path_a = self.get_path_from_root(a)?;
        let path_b = self.get_path_from_root(b)?;

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

    /// Calculate the distance between two nodes.
    /// Returns a tuple (weighted_distance, topological_distance).
    /// weighted_distance: Sum of edge lengths.
    /// topological_distance: Number of edges.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// let n0 = tree.add_node();
    /// let n1 = tree.add_node();
    /// let n2 = tree.add_node();
    /// tree.get_node_mut(n1).unwrap().length = Some(1.5);
    /// tree.get_node_mut(n2).unwrap().length = Some(2.5);
    /// tree.add_child(n0, n1);
    /// tree.add_child(n1, n2);
    ///
    /// let (w, t) = tree.get_distance(&n0, &n2).unwrap();
    /// assert_eq!(w, 4.0);
    /// assert_eq!(t, 2);
    ///
    /// // Error: Unreachable nodes
    /// let orphan = tree.add_node();
    /// assert!(tree.get_distance(&n0, &orphan).is_err());
    /// ```
    pub fn get_distance(&self, a: &NodeId, b: &NodeId) -> Result<(f64, usize), String> {
        let lca = self.get_common_ancestor(a, b)?;

        let dist_to_lca = |start: &NodeId, end: &NodeId| -> (f64, usize) {
            let mut weighted = 0.0;
            let mut topo = 0;
            let mut curr = *start;

            while curr != *end {
                if let Some(node) = self.get_node(curr) {
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

    /// Perform a postorder traversal starting from a given node.
    /// Returns a vector of NodeIds in visitation order (children before parent).
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// //    0
    /// //   / \
    /// //  1   2
    /// let n0 = tree.add_node();
    /// let n1 = tree.add_node();
    /// let n2 = tree.add_node();
    /// tree.add_child(n0, n1);
    /// tree.add_child(n0, n2);
    ///
    /// let traversal = tree.postorder(&n0).unwrap();
    /// assert_eq!(traversal, vec![n1, n2, n0]);
    /// ```
    pub fn postorder(&self, start_node: &NodeId) -> Result<Vec<NodeId>, String> {
        if self.get_node(*start_node).is_none() {
            return Err(format!("Node {} not found", start_node));
        }

        let mut result = Vec::new();
        let mut stack = vec![*start_node];

        // Using a second stack to reverse the order
        let mut output_stack = Vec::new();

        while let Some(curr) = stack.pop() {
            output_stack.push(curr);

            if let Some(node) = self.get_node(curr) {
                // Push children left-to-right, so they are popped right-to-left
                // But since we are building output_stack to be reversed later,
                // we want Right child to be pushed to output_stack BEFORE Left child?
                // Wait.
                // Preorder: Root, Left, Right. Stack: push Right, push Left. Pop Left, Pop Right.
                // Postorder: Left, Right, Root.
                // Reverse Postorder: Root, Right, Left.
                // So if we do Preorder but visit Right child before Left child:
                // Stack: push Left, push Right. Pop Right, Pop Left.
                // Sequence: Root, Right, Left.
                // Then reverse it -> Left, Right, Root.
                for &child in node.children.iter() {
                    stack.push(child);
                }
            }
        }

        while let Some(curr) = output_stack.pop() {
            result.push(curr);
        }

        Ok(result)
    }

    /// Perform a level-order (breadth-first) traversal starting from a given node.
    /// Returns a vector of NodeIds in visitation order.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// //    0
    /// //   / \
    /// //  1   2
    /// // /
    /// // 3
    /// let n0 = tree.add_node();
    /// let n1 = tree.add_node();
    /// let n2 = tree.add_node();
    /// let n3 = tree.add_node();
    /// tree.add_child(n0, n1);
    /// tree.add_child(n0, n2);
    /// tree.add_child(n1, n3);
    ///
    /// let traversal = tree.levelorder(&n0).unwrap();
    /// assert_eq!(traversal, vec![n0, n1, n2, n3]);
    /// ```
    pub fn levelorder(&self, start_node: &NodeId) -> Result<Vec<NodeId>, String> {
        if self.get_node(*start_node).is_none() {
            return Err(format!("Node {} not found", start_node));
        }

        let mut result = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(*start_node);

        while let Some(curr) = queue.pop_front() {
            result.push(curr);
            if let Some(node) = self.get_node(curr) {
                for &child in &node.children {
                    queue.push_back(child);
                }
            }
        }

        Ok(result)
    }

    /// Get the height of a node (maximum distance to its descendant leaves).
    /// Returns 0.0 if the node is a leaf.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// let n0 = tree.add_node();
    /// let n1 = tree.add_node();
    /// let n2 = tree.add_node();
    /// tree.add_child(n0, n1).unwrap();
    /// tree.add_child(n1, n2).unwrap();
    /// tree.get_node_mut(n1).unwrap().length = Some(1.0);
    /// tree.get_node_mut(n2).unwrap().length = Some(2.0);
    ///
    /// // n0 -> n1 (1.0) -> n2 (2.0)
    /// // Height of n0 = dist(n0, n2) = 3.0
    /// assert_eq!(tree.get_height(&n0).unwrap(), 3.0);
    /// assert_eq!(tree.get_height(&n2).unwrap(), 0.0);
    ///
    /// // Simple case
    /// let newick = "(A:1,(B:2)C:1);";
    /// let mut tree = Tree::from_newick(newick).unwrap();
    ///
    /// let id_c = tree.get_node_by_name("C").unwrap();
    /// assert_eq!(tree.get_height(&id_c).unwrap(), 2.0);
    ///
    /// let id_b = tree.get_node_by_name("B").unwrap();
    /// assert_eq!(tree.get_height(&id_b).unwrap(), 0.0);
    /// ```
    pub fn get_height(&self, id: &NodeId) -> Result<f64, String> {
        if self.get_node(*id).is_none() {
            return Err(format!("Node {} not found", id));
        }
        if self.get_node(*id).unwrap().children.is_empty() {
            return Ok(0.0);
        }
        let descendants = self.get_subtree(id)?;
        let mut max_dist = 0.0;
        for desc_id in descendants {
            if desc_id == *id {
                continue;
            }
            if let Some(node) = self.get_node(desc_id) {
                if node.children.is_empty() {
                    if let Ok((dist, _)) = self.get_distance(id, &desc_id) {
                        if dist > max_dist {
                            max_dist = dist;
                        }
                    }
                }
            }
        }
        Ok(max_dist)
    }

    // ###################
    // # TREE STATISTICS #
    // ###################

    /// Get names of all leaves.
    pub fn get_leaf_names(&self) -> Vec<Option<String>> {
        self.get_leaves()
            .iter()
            .map(|&id| self.get_node(id).unwrap().name.clone())
            .collect()
    }

    /// Check if the tree is binary.
    /// A binary tree is rooted and every internal node has at most 2 children.
    pub fn is_binary(&self) -> bool {
        for node in &self.nodes {
            if !node.deleted && node.children.len() > 2 {
                return false;
            }
        }
        true
    }

    /// Get the set of bipartitions defined by the tree.
    /// Each bipartition is represented by the set of leaf names in the clade defined by a node.
    ///
    /// # Errors
    /// Returns error if any leaf is unnamed.
    pub fn get_partitions(
        &self,
    ) -> Result<std::collections::HashSet<std::collections::BTreeSet<String>>, String> {
        let mut partitions = std::collections::HashSet::new();

        // 1. Check names
        for node in &self.nodes {
            if !node.deleted && node.children.is_empty() && node.name.is_none() {
                return Err("All leaves must be named to calculate partitions".to_string());
            }
        }

        if self.root.is_none() {
            return Ok(partitions);
        }

        // 2. Postorder traversal
        let traversal = self.postorder(&self.root.unwrap())?;
        let mut node_leaves: std::collections::HashMap<
            NodeId,
            std::collections::BTreeSet<String>,
        > = std::collections::HashMap::new();

        for &id in &traversal {
            let node = self.get_node(id).unwrap();
            let mut leaves = std::collections::BTreeSet::new();

            if node.children.is_empty() {
                // Leaf
                leaves.insert(node.name.clone().unwrap());
            } else {
                // Internal node: union of children's leaves
                for &child_id in &node.children {
                    if let Some(child_set) = node_leaves.remove(&child_id) {
                        leaves.extend(child_set);
                    }
                }
            }

            // Add to partitions (must clone here because we also insert into map for parent)
            partitions.insert(leaves.clone());

            // Store for parent
            node_leaves.insert(id, leaves);
        }

        Ok(partitions)
    }

    /// Calculate the diameter of the tree (maximum distance between any pair of leaves).
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let t = Tree::from_newick("((A:1,B:2):1,C:4);").unwrap();
    /// // Dist(A,B) = 3
    /// // Dist(A,C) = 1+1+4 = 6
    /// // Dist(B,C) = 2+1+4 = 7
    /// assert_eq!(t.diameter().unwrap(), 7.0);
    /// ```
    pub fn diameter(&self) -> Result<f64, String> {
        use itertools::Itertools;

        let leaves = self.get_leaves();
        if leaves.is_empty() {
            return Err("Tree is empty".to_string());
        }

        // If only 1 leaf, diameter is 0
        if leaves.len() == 1 {
            return Ok(0.0);
        }

        leaves
            .iter()
            .tuple_combinations()
            .map(|(a, b)| {
                self.get_distance(a, b).map(|(w, _)| w).unwrap_or(0.0)
            })
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .ok_or("Could not calculate diameter".to_string())
    }

    /// Calculate Robinson-Foulds distance to another tree.
    /// This is the symmetric difference of the sets of bipartitions (clades).
    /// RF = |A| + |B| - 2 * |A ∩ B|
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let t1 = Tree::from_newick("((A,B),C);").unwrap();
    /// let t2 = Tree::from_newick("((A,C),B);").unwrap();
    /// assert_eq!(t1.robinson_foulds(&t2).unwrap(), 2);
    /// ```
    pub fn robinson_foulds(&self, other: &Tree) -> Result<usize, String> {
        let p1 = self.get_partitions()?;
        let p2 = other.get_partitions()?;

        let get_leaves_from_partitions = |p: &std::collections::HashSet<std::collections::BTreeSet<String>>| -> std::collections::BTreeSet<String> {
            p.iter()
                .filter(|s| s.len() == 1)
                .flat_map(|s| s.iter())
                .cloned()
                .collect()
        };

        let l1 = get_leaves_from_partitions(&p1);
        let l2 = get_leaves_from_partitions(&p2);

        if l1 != l2 {
            return Err("Trees have different leaf sets".to_string());
        }

        let intersection = p1.intersection(&p2).count();
        let rf = p1.len() + p2.len() - 2 * intersection;
        Ok(rf)
    }

    /// Check if a set of nodes forms a monophyletic group (clade).
    ///
    /// A group is monophyletic if it includes a common ancestor and ALL of its descendants.
    /// In this implementation, we check if the set of leaf descendants of the
    /// lowest common ancestor (LCA) of the input nodes matches the input set exactly.
    /// Note: This assumes the input set consists of leaf nodes (tips).
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    ///
    /// //      Root
    /// //     /    \
    /// //    I1     C
    /// //   /  \
    /// //  A    B
    /// let newick = "((A,B)I1,C)Root;";
    /// let tree = Tree::from_newick(newick).unwrap();
    /// let map = tree.get_name_id();
    ///
    /// let id_a = *map.get("A").unwrap();
    /// let id_b = *map.get("B").unwrap();
    /// let id_c = *map.get("C").unwrap();
    ///
    /// let group = vec![id_a, id_b];
    /// assert!(tree.is_monophyletic(&group));
    ///
    /// let group_all = vec![id_a, id_b, id_c];
    /// // {A, B, C} is monophyletic (all descendants of Root)
    /// assert!(tree.is_monophyletic(&group_all));
    ///
    /// let group_ac = vec![id_a, id_c];
    /// // {A, C} is NOT monophyletic (missing B)
    /// assert!(!tree.is_monophyletic(&group_ac));
    /// ```
    pub fn is_monophyletic(&self, ids: &[NodeId]) -> bool {
        if ids.is_empty() {
            return false;
        }

        // Convert input slice to a BTreeSet for final comparison and duplicate handling
        let ids_set: std::collections::BTreeSet<NodeId> = ids.iter().cloned().collect();

        let mut nodes: Vec<NodeId> = ids.iter().cloned().collect();
        let mut sub_root = nodes.pop().unwrap(); // Safe due to is_empty check

        // 1. Find LCA of all nodes in the set
        for id in &nodes {
            match self.get_common_ancestor(&sub_root, id) {
                Ok(lca) => sub_root = lca,
                Err(_) => return false, // Disconnected or invalid
            }
        }

        // 2. Collect all leaf descendants of the LCA
        let mut descendants = std::collections::BTreeSet::new();
        // Use get_subtree (preorder) to find all descendants
        if let Ok(subtree_nodes) = self.get_subtree(&sub_root) {
            for id in subtree_nodes {
                if let Some(node) = self.get_node(id) {
                    if node.is_leaf() {
                        descendants.insert(id);
                    }
                }
            }
        } else {
            return false;
        }

        // 3. Compare sets
        // Note: The input `ids` should ideally contain only leaves for strict monophyly definition
        // in this context (based on original implementation logic).
        descendants.eq(&ids_set)
    }

    /// Get subtree nodes (including self) in preorder.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// let root = tree.add_node();
    /// assert_eq!(tree.get_subtree(&root).unwrap(), vec![root]);
    /// ```
    pub fn get_subtree(&self, root_id: &NodeId) -> Result<Vec<NodeId>, String> {
        self.preorder(root_id)
    }

    /// Get all leaf nodes in the tree.
    pub fn get_leaves(&self) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|n| !n.deleted && n.children.is_empty())
            .map(|n| n.id)
            .collect()
    }

    /// Find nodes that satisfy a predicate.
    pub fn find_nodes<F>(&self, predicate: F) -> Vec<NodeId>
    where
        F: Fn(&Node) -> bool,
    {
        self.nodes
            .iter()
            .filter(|n| !n.deleted && predicate(n))
            .map(|n| n.id)
            .collect()
    }

    /// Find the first node with the given name.
    pub fn get_node_by_name(&self, name: &str) -> Option<NodeId> {
        self.nodes
            .iter()
            .find(|n| !n.deleted && n.name.as_deref() == Some(name))
            .map(|n| n.id)
    }

    /// Get the root node ID
    pub fn get_root(&self) -> Option<NodeId> {
        self.root
    }

    /// Count the total number of descendant nodes (including children, grandchildren, etc.).
    /// Does not include self.
    pub fn count_descendants(&self, id: NodeId) -> usize {
        let mut count = 0;
        if let Some(node) = self.get_node(id) {
            for &child in &node.children {
                count += 1 + self.count_descendants(child);
            }
        }
        count
    }

    /// Deroot the tree by splicing out one of the root's children if the root is bifurcating.
    /// This effectively merges the two edges connected to the root into a single edge,
    /// removing the root node's structural role and making the tree multifurcating at the top level.
    /// The "heavier" child (with more descendants) is the one collapsed into the root.
    ///
    /// # Errors
    /// Returns error if the tree is empty or the root is not bifurcating (degree != 2).
    pub fn deroot(&mut self) -> Result<(), String> {
        let root = self.root.ok_or("Empty tree")?;
        let children = self.get_node(root).unwrap().children.clone();

        if children.len() != 2 {
            return Err("Root is not bifurcating (degree != 2)".to_string());
        }

        let c1 = children[0];
        let c2 = children[1];

        // Weight = 1 (self) + descendants
        let weight1 = 1 + self.count_descendants(c1);
        let weight2 = 1 + self.count_descendants(c2);

        // Collapse the heavier one. If equal, pick first (c1).
        let target = if weight1 >= weight2 { c1 } else { c2 };
        
        self.collapse_node(target)
    }

    /// Find the node with the longest parent edge.
    /// Used for "Longest Branch" rooting.
    pub fn get_node_with_longest_edge(&self) -> Option<NodeId> {
        self.nodes
            .iter()
            .filter(|n| !n.deleted && n.parent.is_some())
            .max_by(|a, b| {
                let len_a = a.length.unwrap_or(0.0);
                let len_b = b.length.unwrap_or(0.0);
                len_a.partial_cmp(&len_b).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|n| n.id)
    }

    /// Reroot the tree at the specified node.
    /// This reverses the direction of edges along the path from the old root to the new root.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// //    0 (root)
    /// //   / \
    /// //  1   2
    /// // /
    /// //3
    /// let n0 = tree.add_node();
    /// let n1 = tree.add_node();
    /// let n2 = tree.add_node();
    /// let n3 = tree.add_node();
    ///
    /// tree.set_root(n0);
    /// tree.add_child(n0, n1).unwrap();
    /// tree.add_child(n0, n2).unwrap();
    /// tree.add_child(n1, n3).unwrap();
    /// tree.get_node_mut(n1).unwrap().length = Some(10.0);
    /// tree.get_node_mut(n3).unwrap().length = Some(5.0);
    ///
    /// // Reroot at 3
    /// tree.reroot_at(n3, false).unwrap();
    ///
    /// // New structure:
    /// // 3 (root) -> 1 (len 5.0) -> 0 (len 10.0) -> 2 (len None/Default)
    /// assert_eq!(tree.get_root(), Some(n3));
    /// assert_eq!(tree.get_node(n3).unwrap().parent, None);
    /// assert_eq!(tree.get_node(n3).unwrap().length, None);
    ///
    /// let node1 = tree.get_node(n1).unwrap();
    /// assert_eq!(node1.parent, Some(n3));
    /// assert_eq!(node1.length, Some(5.0));
    /// assert!(node1.children.contains(&n0));
    ///
    /// let node0 = tree.get_node(n0).unwrap();
    /// assert_eq!(node0.parent, Some(n1));
    /// assert_eq!(node0.length, Some(10.0));
    /// assert!(node0.children.contains(&n2));
    ///
    /// // Error: Node not found
    /// assert!(tree.reroot_at(9999, false).is_err());
    /// ```
    pub fn reroot_at(&mut self, new_root_id: NodeId, process_support_values: bool) -> Result<(), String> {
        if self.get_node(new_root_id).is_none() {
            return Err(format!("Node {} not found", new_root_id));
        }

        let old_root_id = self.root.ok_or("Tree has no root")?;
        if old_root_id == new_root_id {
            return Ok(());
        }

        // 1. Get path from old root to new root
        let path = self.get_path_from_root(&new_root_id)?;

        // 1.5 Process Support Values (Labels)
        // Shift internal node labels along the path to align with edge reversals
        if process_support_values {
            let new_root_is_leaf = self.get_node(new_root_id).map(|n| n.children.is_empty()).unwrap_or(false);
            
            // Capture original names
            let names: Vec<Option<String>> = path.iter()
                .map(|&id| self.get_node(id).unwrap().name.clone())
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

                    if let Some(node) = self.get_node_mut(node_id) {
                        node.name = new_name;
                    }
                }
            }
        }

        // 2. Collect edge lengths along the path
        // path[i]'s length represents edge (path[i-1] -> path[i])
        let mut lengths = Vec::new();
        for &id in &path {
            lengths.push(self.get_node(id).unwrap().length);
        }

        // 3. Reverse edges
        for i in (1..path.len()).rev() {
            let child_id = path[i];
            let parent_id = path[i - 1];
            let length = lengths[i];

            // a. Remove child from parent's children
            if let Some(parent) = self.nodes.get_mut(parent_id) {
                parent.children.retain(|&x| x != child_id);
            }

            // b. Add parent to child's children
            if let Some(child) = self.nodes.get_mut(child_id) {
                child.children.push(parent_id);
            }

            // c. Update parent's parent pointer and length
            if let Some(parent) = self.nodes.get_mut(parent_id) {
                parent.parent = Some(child_id);
                parent.length = length;
            }
        }

        // 4. Finalize new root
        if let Some(new_root) = self.nodes.get_mut(new_root_id) {
            new_root.parent = None;
            new_root.length = None;
        }

        self.root = Some(new_root_id);

        Ok(())
    }

    /// Prune nodes that match a predicate.
    /// Warning: This removes the matching nodes AND their descendants.
    pub fn prune_where<F>(&mut self, predicate: F)
    where
        F: Fn(&Node) -> bool,
    {
        // We need to collect IDs first to avoid borrowing issues
        let to_remove: Vec<NodeId> = self
            .nodes
            .iter()
            .filter(|n| !n.deleted && predicate(n))
            .map(|n| n.id)
            .collect();

        for id in to_remove {
            self.remove_node(id, true);
        }
    }

    /// Get number of active nodes
    pub fn len(&self) -> usize {
        self.nodes.iter().filter(|n| !n.deleted).count()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get node names
    ///
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    ///
    /// let newick = "((A,B)D,C);";
    /// let tree = Tree::from_newick(newick).unwrap();
    /// assert_eq!(tree.get_names(), vec!["D".to_string(),"A".to_string(),"B".to_string(),"C".to_string(), ]);
    /// ```
    pub fn get_names(&self) -> Vec<String> {
        if let Some(root) = self.root {
            self.preorder(&root)
                .unwrap_or_default()
                .iter()
                .filter_map(|&id| self.get_node(id))
                .filter_map(|node| node.name.clone())
                .collect()
        } else {
            Vec::new()
        }
    }

    /// Get hash of name-id
    ///
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    ///
    /// let newick = "((A,B),C);";
    /// let tree = Tree::from_newick(newick).unwrap();
    /// let id_of = tree.get_name_id();
    /// assert_eq!(*id_of.get("A").unwrap(), 2usize);
    /// ```
    pub fn get_name_id(&self) -> BTreeMap<String, usize> {
        let mut id_of = BTreeMap::new();
        if let Some(root) = self.root {
            for id in self.preorder(&root).unwrap_or_default().iter() {
                if let Some(node) = self.get_node(*id) {
                    if let Some(name) = &node.name {
                        id_of.insert(name.clone(), *id);
                    }
                }
            }
        }
        id_of
    }

    /// Get all values for a specific property key across the tree.
    /// Returns a map of NodeId -> Property Value.
    ///
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    ///
    /// let newick = "((A[&&NHX:S=Human],B),C[&&NHX:S=Chimp]);";
    /// let tree = Tree::from_newick(newick).unwrap();
    /// let species_map = tree.get_property_values("S");
    ///
    /// let map = tree.get_name_id();
    /// let id_a = *map.get("A").unwrap();
    /// let id_c = *map.get("C").unwrap();
    ///
    /// assert_eq!(species_map.get(&id_a).unwrap(), "Human");
    /// assert_eq!(species_map.get(&id_c).unwrap(), "Chimp");
    /// assert!(!species_map.contains_key(map.get("B").unwrap()));
    /// ```
    pub fn get_property_values(&self, key: &str) -> BTreeMap<NodeId, String> {
        let mut values = BTreeMap::new();
        if let Some(root) = self.root {
            for id in self.preorder(&root).unwrap_or_default().iter() {
                if let Some(node) = self.get_node(*id) {
                    if let Some(val) = node.get_property(key) {
                        values.insert(*id, val.clone());
                    }
                }
            }
        }
        values
    }

    /// Serialize the tree to a Newick string (compact format).
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// let root = tree.add_node();
    /// tree.set_root(root);
    /// tree.get_node_mut(root).unwrap().set_name("A");
    /// assert_eq!(tree.to_newick(), "A;");
    /// ```
    pub fn to_newick(&self) -> String {
        writer::write_newick(self)
    }

    /// Serialize the tree to a Newick string with optional indentation.
    ///
    /// # Arguments
    /// * `indent` - The string to use for indentation (e.g., "  ", "\t").
    ///              If empty, output will be compact (no whitespace).
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    /// let mut tree = Tree::new();
    /// let root = tree.add_node();
    /// let child = tree.add_node();
    /// tree.set_root(root);
    /// tree.add_child(root, child);
    /// tree.get_node_mut(root).unwrap().set_name("Root");
    /// tree.get_node_mut(child).unwrap().set_name("Child");
    ///
    /// let expected = "(\n  Child\n)Root;";
    /// assert_eq!(tree.to_newick_with_format("  "), expected);
    /// ```
    pub fn to_newick_with_format(&self, indent: &str) -> String {
        writer::write_newick_with_format(self, indent)
    }

    /// Serialize the tree to Graphviz DOT format.
    pub fn to_dot(&self) -> String {
        writer::write_dot(self)
    }
}

impl Tree {
    // #####################
    // # TREE MANIPULATION #
    // #####################

    /// Insert a node in the middle of the desired node and its parent.
    /// Returns the new parent node ID.
    pub fn insert_parent(&mut self, id: NodeId) -> Result<NodeId, String> {
        let node = self.get_node(id).ok_or(format!("Node {} not found", id))?;
        let parent = node.parent.ok_or("Node has no parent")?;
        let length = node.length;
        let new_length = length.map(|l| l / 2.0);

        let new_node = self.add_node();

        // Link parent -> new_node
        self.add_child(parent, new_node)?;
        if let Some(n) = self.get_node_mut(new_node) {
            n.length = new_length;
        }

        // Unlink parent -> id
        if let Some(p_node) = self.get_node_mut(parent) {
            p_node.children.retain(|&c| c != id);
        }
        // Update id parent
        if let Some(node) = self.get_node_mut(id) {
            node.parent = None;
        }

        // Link new_node -> id
        self.add_child(new_node, id)?;
        if let Some(node) = self.get_node_mut(id) {
            node.length = new_length;
        }

        Ok(new_node)
    }

    /// Swap parent-child link of a node.
    /// This reverses the edge between the node and its parent.
    pub fn swap_parent(
        &mut self,
        id: NodeId,
        _prev_edge: Option<f64>,
    ) -> Result<Option<f64>, String> {
        let node = self.get_node(id).ok_or(format!("Node {} not found", id))?;
        let parent = node.parent.ok_or("Node has no parent")?;

        // Swap lengths
        let child_len = node.length;
        let parent_len = self.get_node(parent).ok_or("Parent not found")?.length;

        if let Some(p_node) = self.get_node_mut(parent) {
            p_node.length = child_len;
        }
        if let Some(node) = self.get_node_mut(id) {
            node.length = parent_len;
        }

        // Unlink parent -> id
        if let Some(p_node) = self.get_node_mut(parent) {
            p_node.children.retain(|&c| c != id);
        }
        // Unlink id -> parent
        if let Some(node) = self.get_node_mut(id) {
            node.parent = None;
        }

        // Link id -> parent (parent becomes child)
        // We must clear parent's parent pointer to satisfy add_child check (as parent is "moving down")
        if let Some(p_node) = self.get_node_mut(parent) {
            p_node.parent = None;
        }

        self.add_child(id, parent)?;

        Ok(None)
    }

    /// Insert a new parent node for a pair of nodes (LCA-based).
    /// Returns the new parent node ID.
    pub fn insert_parent_pair(&mut self, id1: NodeId, id2: NodeId) -> Result<NodeId, String> {
        let old = self.get_common_ancestor(&id1, &id2)?;

        // Get original edge lengths
        let edge1 = self.get_node(id1).ok_or("Node 1 not found")?.length;
        let edge2 = self.get_node(id2).ok_or("Node 2 not found")?.length;

        // New node with parent (old) has no edge length
        let new = self.add_node();
        self.add_child(old, new)?;

        // Move children to new node
        // 1. Unlink from their current parents
        let p1 = self.get_node(id1).and_then(|n| n.parent);
        if let Some(p) = p1 {
            if let Some(p_node) = self.get_node_mut(p) {
                p_node.children.retain(|&c| c != id1);
            }
        }
        
        let p2 = self.get_node(id2).and_then(|n| n.parent);
        if let Some(p) = p2 {
            if let Some(p_node) = self.get_node_mut(p) {
                p_node.children.retain(|&c| c != id2);
            }
        }

        if let Some(node) = self.get_node_mut(id1) {
            node.parent = None;
        }
        if let Some(node) = self.get_node_mut(id2) {
            node.parent = None;
        }

        // 2. Link to new
        self.add_child(new, id1)?;
        if let Some(node) = self.get_node_mut(id1) {
            node.length = edge1;
        }

        self.add_child(new, id2)?;
        if let Some(node) = self.get_node_mut(id2) {
            node.length = edge2;
        }

        Ok(new)
    }

    /// Remove nodes that have a parent and exactly one child (degree 2 nodes).
    /// This is often used after rerooting to clean up the tree.
    pub fn remove_degree_two_nodes(&mut self) {
        loop {
            // Find a node that is:
            // 1. Not deleted
            // 2. Has a parent (not root)
            // 3. Has exactly 1 child
            let to_remove = if let Some(_root) = self.get_root() {
                self.find_nodes(|n| {
                    n.parent.is_some() && // Not root
                    n.children.len() == 1 // Degree 2 (1 parent, 1 child)
                }).first().cloned()
            } else {
                None
            };
            
            if let Some(id) = to_remove {
                // Ignore result, just proceed
                let _ = self.collapse_node(id);
            } else {
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let sub = tree.get_subtree(&n1).unwrap();
        assert_eq!(sub, vec![n1, n3]);

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
        // Leaves: A, B, D
        // Partitions (clades):
        // - Leaf A: {A}
        // - Leaf B: {B}
        // - Node C: {A, B}
        // - Leaf D: {D}
        // - Root E: {A, B, D}

        let newick = "((A,B)C,D)E;";
        let tree = Tree::from_newick(newick).unwrap();

        let partitions = tree.get_partitions().unwrap();

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

        // Check strict equality of the set of sets
        // Since BTreeSet and HashSet implement Eq correctly based on content
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
        // (A:1,B:1)Root:1;
        // Children of Root: A (weight 1+0=1), B (weight 1+0=1)
        // Deroot should collapse A (since weight A >= weight B, and A comes first)
        let mut tree = Tree::from_newick("(A:1,B:1)Root:1;").unwrap();
        
        tree.deroot().unwrap();
        
        // After collapse A:
        // Root is gone/merged. Actually collapse_node removes the node and links children to parent.
        // But wait, collapse_node(target) removes 'target' and connects its children to 'target.parent'.
        // In deroot(), we collapse a child of the root.
        // If we collapse A (child of Root), A's children (none) become children of Root.
        // A is removed.
        // This effectively removes A? No.
        
        // Let's re-read deroot logic.
        // "The 'heavier' child (with more descendants) is the one collapsed into the root."
        // collapse_node(target): target is removed, target's children become children of target's parent (Root).
        // target's edge length is added to children's edge length.
        
        // Example: (A:1, B:1)Root
        // Weights: A=1, B=1. Target=A.
        // Collapse A.
        // A is a leaf. Children = [].
        // Root children: [B].
        // This seems wrong for "derooting" a bifurcating root to make it trifurcating?
        // Usually derooting means: (A,B,C); -> Root has A, B, C.
        // If we have ((A,B)C, D)Root; -> Remove C. Root has A, B, D.
        
        // If we have (A,B)Root; -> removing A makes Root have B. That's not derooting.
        // Derooting usually implies the root has 2 children, and one of them is an internal node.
        // We collapse that internal node to make the root have > 2 children.
        
        // Let's try: ((A,B)C,D)Root;
        // Root children: C, D.
        // C descendants: A, B. Weight(C) = 1 + 2 = 3.
        // D descendants: []. Weight(D) = 1.
        // Target = C.
        // Collapse C.
        // C is removed. C's children (A, B) become children of Root.
        // Root children: A, B, D.
        // Result: (A,B,D)Root;
        
        let mut tree = Tree::from_newick("((A:1,B:2)C:3,D:4)Root;").unwrap();
        tree.deroot().unwrap();
        
        let root = tree.get_root().unwrap();
        let children = &tree.get_node(root).unwrap().children;
        assert_eq!(children.len(), 3);
        
        // Check names of children
        let child_names: Vec<String> = children.iter()
            .map(|&id| tree.get_node(id).unwrap().name.clone().unwrap_or_default())
            .collect();
        
        // Order might depend on splice. splice replaces C with A,B.
        // Original: [C, D]. Replace C -> [A, B, D].
        assert!(child_names.contains(&"A".to_string()));
        assert!(child_names.contains(&"B".to_string()));
        assert!(child_names.contains(&"D".to_string()));
    }

    #[test]
    fn test_reroot_support_values() {
        // Tree: (A, (B, C)Support)Root;
        // Reroot at C.
        // Path: Root -> Support -> C.
        // New Root: C.
        // Old Root becomes child of Support?
        // Support label should move?
        
        // Let's look at `reroot_at` logic for support values.
        // It shifts labels along the path.
        // Path: [Root, Support, C]
        // i=0 (Root): new_name = names[1] (Support).
        // i=1 (Support): new_name = names[2] (C is leaf? yes). So None?
        // i=2 (C): is_leaf = true. No change.
        
        // Wait, "Only modify internal nodes (leaves keep Taxon names)".
        // C is new root. C was leaf.
        // new_root_is_leaf = true.
        
        // i=0 (Root): is_leaf = false.
        //   new_name = names[1] ("Support").
        //   Root name becomes "Support".
        
        // i=1 (Support): is_leaf = false?
        //   i < path.len()-1 (1 < 2).
        //   next_is_leaf = (2 == 2) && true = true.
        //   new_name = None.
        //   Support name becomes None.
        
        // i=2 (C): is_leaf = true. Skipped.
        
        // Result:
        // Root (now child of Support) -> name="Support"
        // Support (now child of C) -> name=None
        // C (root) -> name="C"
        
        let mut tree = Tree::from_newick("(A,(B,C)Support)Root;").unwrap();
        let c_id = tree.get_node_by_name("C").unwrap();
        
        tree.reroot_at(c_id, true).unwrap();
        
        // C is root
        assert_eq!(tree.get_root(), Some(c_id));
        
        // Old root should be named "Support"
        let old_root_id = tree.get_node_by_name("Support").unwrap();
        let _old_root = tree.get_node(old_root_id).unwrap();
        // Wait, get_node_by_name uses current names.
        // The node that WAS Root should now be named Support.
        // The node that WAS Support should now be named None (so not found by name "Support").
        
        // Let's find by ID if possible, but IDs are internal.
        // Let's check structure.
        // C -> SupportNode -> RootNode -> A
        //                  -> B
        
        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        assert_eq!(root.name.as_deref(), Some("C"));
        
        let support_node_id = root.children[0]; // The old Support node
        let support_node = tree.get_node(support_node_id).unwrap();
        assert_eq!(support_node.name, None);
        
        let old_root_node_id = support_node.children.iter().find(|&&id| {
             // Find the one that has A as child
             let n = tree.get_node(id).unwrap();
             n.children.iter().any(|&child| tree.get_node(child).unwrap().name.as_deref() == Some("A"))
        }).unwrap();
        
        let old_root_node = tree.get_node(*old_root_node_id).unwrap();
        assert_eq!(old_root_node.name.as_deref(), Some("Support"));
    }

    #[test]
    fn test_reroot_longest_branch() {
        // (A:1, B:2)Root;
        // Longest branch is B (len 2).
        // Reroot should pick B.
        let mut tree = Tree::from_newick("(A:1,B:2)Root;").unwrap();
        
        // We need to implement default reroot logic if we want to test it here,
        // but `reroot_at` takes an ID.
        // The logic for "default" is in CLI, but `get_node_with_longest_edge` is in Tree.
        
        let target = tree.get_node_with_longest_edge().unwrap();
        let b_id = tree.get_node_by_name("B").unwrap();
        assert_eq!(target, b_id);
        
        tree.reroot_at(target, false).unwrap();
        assert_eq!(tree.get_root(), Some(b_id));
    }
}
