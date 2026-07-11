pub mod algo;
pub mod balance;
pub mod distance;
pub mod io;
pub mod ops;
pub mod query;
pub mod stat;
pub mod support;
#[cfg(test)]
pub mod tests;
pub mod traversal;

use super::node::{Node, NodeId};
use std::collections::{BTreeMap, BTreeSet};

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

    /// Add a child node under a parent.
    pub fn add_child(&mut self, parent_id: NodeId, child_id: NodeId) -> anyhow::Result<()> {
        ops::add_child(self, parent_id, child_id)
    }

    /// Remove a node, optionally recursively removing its subtree.
    pub fn remove_node(&mut self, id: NodeId, recursive: bool) {
        ops::remove_node(self, id, recursive)
    }

    /// Collapse an internal node by splicing its children to its parent.
    pub fn collapse_node(&mut self, id: NodeId) -> anyhow::Result<()> {
        ops::collapse_node(self, id)
    }

    /// Permanently remove soft-deleted nodes and reclaim arena slots.
    pub fn compact(&mut self) {
        ops::compact(self)
    }

    /// Insert a new parent above `child_id`, returning the new node's ID.
    pub fn insert_parent(&mut self, child_id: NodeId) -> anyhow::Result<NodeId> {
        ops::insert_parent(self, child_id)
    }

    /// Swap a node with its parent, returning the previous edge length.
    pub fn swap_parent(
        &mut self,
        node_id: NodeId,
        prev_edge: Option<f64>,
    ) -> anyhow::Result<Option<f64>> {
        ops::swap_parent(self, node_id, prev_edge)
    }

    /// Insert a common parent above two sibling nodes.
    pub fn insert_parent_pair(&mut self, id1: NodeId, id2: NodeId) -> anyhow::Result<NodeId> {
        ops::insert_parent_pair(self, id1, id2)
    }

    /// Remove all internal nodes with exactly one child.
    pub fn remove_degree_two_nodes(&mut self) {
        ops::remove_degree_two_nodes(self)
    }

    /// Deroot the tree by splicing out the root node.
    pub fn deroot(&mut self) -> anyhow::Result<()> {
        ops::deroot(self)
    }

    /// Reroot the tree at `new_root_id`, optionally shifting support values.
    pub fn reroot_at(
        &mut self,
        new_root_id: NodeId,
        process_support_values: bool,
    ) -> anyhow::Result<()> {
        ops::reroot_at(self, new_root_id, process_support_values)
    }

    /// Remove all nodes matching `predicate`.
    pub fn prune_where<F>(&mut self, predicate: F)
    where
        F: Fn(&Node) -> bool + Copy,
    {
        ops::prune_where(self, predicate)
    }

    /// Condense the subtree rooted at `sub_root_id` into a single named node.
    pub fn condense_subtree(
        &mut self,
        sub_root_id: NodeId,
        name: &str,
        member_count: usize,
    ) -> anyhow::Result<()> {
        ops::condense_subtree(self, sub_root_id, name, member_count)
    }

    // --- Delegation to traversal ---

    /// Pre-order traversal returning node IDs.
    pub fn preorder(&self, start_node: &NodeId) -> Vec<NodeId> {
        traversal::preorder(self, *start_node)
    }

    /// Post-order traversal returning node IDs.
    pub fn postorder(&self, start_node: &NodeId) -> Vec<NodeId> {
        traversal::postorder(self, *start_node)
    }

    /// Level-order (BFS) traversal returning node IDs.
    pub fn levelorder(&self, start_node: &NodeId) -> Vec<NodeId> {
        traversal::levelorder(self, *start_node)
    }

    /// Extract a subtree rooted at `root_id` as a new Tree.
    pub fn extract_subtree(&self, root_id: &NodeId) -> anyhow::Result<Tree> {
        traversal::extract_subtree(self, *root_id)
    }

    /// Return all node IDs in the subtree rooted at `root_id`.
    pub fn get_subtree(&self, root_id: &NodeId) -> Vec<NodeId> {
        traversal::preorder(self, *root_id)
    }

    // --- Delegation to query ---

    /// Return the path from root to `id` (inclusive).
    pub fn get_path_from_root(&self, id: &NodeId) -> anyhow::Result<Vec<NodeId>> {
        query::get_path_from_root(self, id)
    }

    /// Find the Lowest Common Ancestor (LCA) of two nodes.
    pub fn get_common_ancestor(&self, a: &NodeId, b: &NodeId) -> anyhow::Result<NodeId> {
        query::get_common_ancestor(self, a, b)
    }

    /// Find the Lowest Common Ancestor (LCA) of multiple nodes.
    pub fn get_lca(&self, nodes: &[NodeId]) -> anyhow::Result<NodeId> {
        query::get_lca(self, nodes)
    }

    /// Return (weighted distance, edge count) between two nodes.
    pub fn get_distance(&self, a: &NodeId, b: &NodeId) -> anyhow::Result<(f64, usize)> {
        query::get_distance(self, a, b)
    }

    /// Distance between two nodes. Uses branch lengths if non-zero, else edge count.
    pub fn node_distance(&self, a: &NodeId, b: &NodeId) -> anyhow::Result<f64> {
        query::node_distance(self, a, b)
    }

    /// Find all nodes matching a predicate.
    pub fn find_nodes<F>(&self, predicate: F) -> Vec<NodeId>
    where
        F: Fn(&Node) -> bool,
    {
        query::find_nodes(self, predicate)
    }

    /// Get node ID by name, returning the first match.
    pub fn get_node_by_name(&self, name: &str) -> Option<NodeId> {
        query::get_node_by_name(self, name)
    }

    /// Check if a set of nodes forms a monophyletic clade.
    pub fn is_monophyletic(&self, nodes: &[NodeId]) -> bool {
        query::is_monophyletic(self, nodes)
    }

    /// Collect IDs of all named leaves in the subtree rooted at `id`.
    pub fn get_named_leaves(&self, id: NodeId) -> BTreeSet<NodeId> {
        query::get_named_leaves(self, id)
    }

    /// Get height of a node (max distance to any leaf in its subtree).
    pub fn get_height(&self, id: NodeId, weighted: bool) -> f64 {
        query::get_height(self, id, weighted)
    }

    // --- Delegation to stat ---

    /// Return IDs of all leaves in the tree.
    pub fn get_leaves(&self) -> Vec<NodeId> {
        if let Some(root) = self.root {
            stat::get_leaves(self, root)
        } else {
            Vec::new()
        }
    }

    /// Return names of all leaves in the tree.
    pub fn get_leaf_names(&self) -> Vec<Option<String>> {
        if let Some(root) = self.root {
            stat::get_leaf_names(self, root)
        } else {
            Vec::new()
        }
    }

    /// Check if the tree is binary (all internal nodes have 2 children).
    pub fn is_binary(&self) -> bool {
        stat::is_binary(self)
    }

    /// Check if the tree is rooted (root has exactly 2 children).
    pub fn is_rooted(&self) -> bool {
        stat::is_rooted(self)
    }

    /// Calculate the tree diameter (longest path between any two nodes).
    pub fn diameter(&self) -> f64 {
        stat::diameter(self, true)
    }

    /// Return all node names in the tree.
    pub fn get_names(&self) -> Vec<String> {
        stat::get_names(self)
    }

    /// Return a map from node name to node ID.
    pub fn get_name_id(&self) -> BTreeMap<String, NodeId> {
        stat::get_name_id(self)
    }

    /// Return a map from node ID to property value for a given key.
    pub fn get_property_values(&self, key: &str) -> BTreeMap<NodeId, String> {
        stat::get_property_values(self, key)
    }

    /// Return the node ID with the longest edge length.
    pub fn get_node_with_longest_edge(&self) -> Option<NodeId> {
        stat::get_node_with_longest_edge(self)
    }

    // --- Delegation to io ---

    /// Read Newick trees from a file.
    pub fn from_file(infile: &str) -> anyhow::Result<Vec<Tree>> {
        io::from_file(infile)
    }

    /// Serialize the tree to a Newick string.
    pub fn to_newick(&self) -> String {
        io::to_newick(self)
    }

    /// Serialize a subtree to a Newick string.
    pub fn to_newick_subtree(&self, root: NodeId) -> String {
        io::to_newick_subtree(self, root, "")
    }

    /// Serialize the tree to a Newick string with indentation.
    pub fn to_newick_with_format(&self, indent: &str) -> String {
        io::to_newick_with_format(self, indent)
    }

    /// Serialize the tree to Graphviz DOT format.
    pub fn to_dot(&self) -> String {
        io::to_dot(self)
    }

    /// Serialize the tree to SVG format.
    pub fn to_svg(&self, height: f64, vskip: f64, width: f64) -> String {
        io::to_svg(self, height, vskip, width)
    }
}
