use super::node::{Node, NodeId};
use std::collections::HashMap;

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
    
    /// Get the root node ID
    pub fn get_root(&self) -> Option<NodeId> {
        self.root
    }
    
    /// Get number of active nodes
    pub fn len(&self) -> usize {
        self.nodes.iter().filter(|n| !n.deleted).count()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
