pub mod io;
pub mod ops;
pub mod query;
pub mod stat;
pub mod support;
#[cfg(test)]
pub mod tests;
pub mod traversal;

use super::node::{Node, NodeId};
use std::collections::BTreeMap;

#[derive(Debug, Default, Clone)]
pub struct Tree {
    /// Arena storage for all nodes
    pub(super) nodes: Vec<Node>,

    /// Optional root ID (a tree might be empty or in construction)
    pub(super) root: Option<NodeId>,
}

impl Tree {
    /// Create a new empty tree
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a new node to the tree. Returns the new node's ID.
    pub fn add_node(&mut self) -> NodeId {
        let id = self.nodes.len();
        let node = Node::new(id);
        self.nodes.push(node);
        id
    }

    /// Get number of nodes
    pub fn len(&self) -> usize {
        self.nodes.iter().filter(|n| !n.deleted).count()
    }

    /// Check if tree is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get root ID
    pub fn get_root(&self) -> Option<NodeId> {
        self.root
    }

    /// Get a reference to a node by ID.
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

    // --- Delegation to ops ---

    pub fn add_child(&mut self, parent_id: NodeId, child_id: NodeId) -> Result<(), String> {
        ops::add_child(self, parent_id, child_id)
    }

    pub fn remove_node(&mut self, id: NodeId, recursive: bool) {
        ops::remove_node(self, id, recursive)
    }

    pub fn collapse_node(&mut self, id: NodeId) -> Result<(), String> {
        ops::collapse_node(self, id)
    }

    pub fn compact(&mut self) {
        ops::compact(self)
    }

    pub fn insert_parent(&mut self, child_id: NodeId) -> Result<NodeId, String> {
        ops::insert_parent(self, child_id)
    }

    pub fn swap_parent(
        &mut self,
        node_id: NodeId,
        prev_edge: Option<f64>,
    ) -> Result<Option<f64>, String> {
        ops::swap_parent(self, node_id, prev_edge)
    }

    pub fn insert_parent_pair(&mut self, id1: NodeId, id2: NodeId) -> Result<NodeId, String> {
        ops::insert_parent_pair(self, id1, id2)
    }

    pub fn remove_degree_two_nodes(&mut self) {
        ops::remove_degree_two_nodes(self)
    }

    pub fn deroot(&mut self) -> Result<(), String> {
        ops::deroot(self)
    }

    pub fn reroot_at(
        &mut self,
        new_root_id: NodeId,
        process_support_values: bool,
    ) -> Result<(), String> {
        ops::reroot_at(self, new_root_id, process_support_values)
    }

    pub fn prune_where<F>(&mut self, predicate: F)
    where
        F: Fn(&Node) -> bool + Copy,
    {
        ops::prune_where(self, predicate)
    }

    // --- Delegation to traversal ---

    pub fn preorder(&self, start_node: &NodeId) -> Result<Vec<NodeId>, String> {
        Ok(traversal::preorder(self, *start_node))
    }

    pub fn postorder(&self, start_node: &NodeId) -> Result<Vec<NodeId>, String> {
        Ok(traversal::postorder(self, *start_node))
    }

    pub fn levelorder(&self, start_node: &NodeId) -> Result<Vec<NodeId>, String> {
        Ok(traversal::levelorder(self, *start_node))
    }

    pub fn extract_subtree(&self, root_id: &NodeId) -> Result<Tree, String> {
        traversal::extract_subtree(self, *root_id)
    }

    pub fn get_subtree(&self, root_id: &NodeId) -> Result<Vec<NodeId>, String> {
        Ok(traversal::preorder(self, *root_id))
    }

    // --- Delegation to query ---

    pub fn get_path_from_root(&self, id: &NodeId) -> Result<Vec<NodeId>, String> {
        query::get_path_from_root(self, id)
    }

    pub fn get_common_ancestor(&self, a: &NodeId, b: &NodeId) -> Result<NodeId, String> {
        query::get_common_ancestor(self, a, b)
    }

    pub fn get_distance(&self, a: &NodeId, b: &NodeId) -> Result<(f64, usize), String> {
        query::get_distance(self, a, b)
    }

    pub fn find_nodes<F>(&self, predicate: F) -> Vec<NodeId>
    where
        F: Fn(&Node) -> bool,
    {
        query::find_nodes(self, predicate)
    }

    pub fn get_node_by_name(&self, name: &str) -> Option<NodeId> {
        query::get_node_by_name(self, name)
    }

    pub fn is_monophyletic(&self, nodes: &[NodeId]) -> bool {
        query::is_monophyletic(self, nodes)
    }

    pub fn get_height(&self, id: NodeId, weighted: bool) -> f64 {
        query::get_height(self, id, weighted)
    }

    // --- Delegation to stat ---

    pub fn get_leaves(&self) -> Vec<NodeId> {
        if let Some(root) = self.root {
            stat::get_leaves(self, root)
        } else {
            Vec::new()
        }
    }

    pub fn get_leaf_names(&self) -> Vec<Option<String>> {
        if let Some(root) = self.root {
            stat::get_leaf_names(self, root)
        } else {
            Vec::new()
        }
    }

    pub fn is_binary(&self) -> bool {
        stat::is_binary(self)
    }

    pub fn is_rooted(&self) -> bool {
        stat::is_rooted(self)
    }

    pub fn diameter(&self) -> Result<f64, String> {
        Ok(stat::diameter(self, true))
    }

    pub fn get_names(&self) -> Vec<String> {
        stat::get_names(self)
    }

    pub fn get_name_id(&self) -> BTreeMap<String, NodeId> {
        stat::get_name_id(self)
    }

    pub fn get_property_values(&self, key: &str) -> BTreeMap<NodeId, String> {
        stat::get_property_values(self, key)
    }

    pub fn get_node_with_longest_edge(&self) -> Option<NodeId> {
        stat::get_node_with_longest_edge(self)
    }

    // --- Delegation to io ---

    pub fn from_file(infile: &str) -> anyhow::Result<Vec<Tree>> {
        io::from_file(infile)
    }

    pub fn to_newick(&self) -> String {
        io::to_newick(self)
    }

    pub fn to_newick_subtree(&self, root: NodeId) -> String {
        io::to_newick_subtree(self, root, "")
    }

    pub fn to_newick_with_format(&self, indent: &str) -> String {
        io::to_newick_with_format(self, indent)
    }

    pub fn to_dot(&self) -> String {
        io::to_dot(self)
    }
}
