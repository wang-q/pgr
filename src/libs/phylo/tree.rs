use super::node::{Node, NodeId};
use std::collections::HashMap;
use super::writer;

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
             return Err(format!("Node {} already has parent {}", child_id, old_parent));
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
            if node.deleted { continue; }
            
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

    /// Get all nodes in the subtree rooted at the specified node (inclusive).
    pub fn get_subtree(&self, root_id: &NodeId) -> Result<Vec<NodeId>, String> {
        self.preorder(root_id)
    }

    /// Get all leaf nodes in the tree.
    pub fn get_leaves(&self) -> Vec<NodeId> {
        self.nodes.iter()
            .filter(|n| !n.deleted && n.children.is_empty())
            .map(|n| n.id)
            .collect()
    }

    /// Find nodes that satisfy a predicate.
    pub fn find_nodes<F>(&self, predicate: F) -> Vec<NodeId> 
    where F: Fn(&Node) -> bool {
        self.nodes.iter()
            .filter(|n| !n.deleted && predicate(n))
            .map(|n| n.id)
            .collect()
    }

    /// Find the first node with the given name.
    pub fn get_node_by_name(&self, name: &str) -> Option<NodeId> {
        self.nodes.iter()
            .find(|n| !n.deleted && n.name.as_deref() == Some(name))
            .map(|n| n.id)
    }

    /// Get the root node ID
    pub fn get_root(&self) -> Option<NodeId> {
        self.root
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
    /// tree.reroot_at(n3).unwrap();
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
    /// assert!(tree.reroot_at(9999).is_err());
    /// ```
    pub fn reroot_at(&mut self, new_root_id: NodeId) -> Result<(), String> {
        if self.get_node(new_root_id).is_none() {
            return Err(format!("Node {} not found", new_root_id));
        }
        
        let old_root_id = self.root.ok_or("Tree has no root")?;
        if old_root_id == new_root_id {
            return Ok(());
        }

        // 1. Get path from old root to new root
        let path = self.get_path_from_root(&new_root_id)?;
        
        // 2. Collect edge lengths along the path
        // path[i]'s length represents edge (path[i-1] -> path[i])
        let mut lengths = Vec::new();
        for &id in &path {
            lengths.push(self.get_node(id).unwrap().length);
        }

        // 3. Reverse edges
        for i in (1..path.len()).rev() {
            let child_id = path[i];
            let parent_id = path[i-1];
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
    where F: Fn(&Node) -> bool {
        // We need to collect IDs first to avoid borrowing issues
        let to_remove: Vec<NodeId> = self.nodes.iter()
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

}
