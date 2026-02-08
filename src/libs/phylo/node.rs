use std::collections::BTreeMap;

/// NodeId is an index into the Tree's node vector.
/// It is lightweight (Copy) and safe (no pointers).
pub type NodeId = usize;

#[derive(Debug, Clone)]
pub struct Node {
    /// Unique identifier for the node (index in the arena)
    pub id: NodeId,
    
    /// Parent node ID (None for root)
    pub parent: Option<NodeId>,
    
    /// List of child node IDs
    pub children: Vec<NodeId>,

    // --- Payload ---

    /// Node name/label (e.g., "human", "internal_1")
    pub name: Option<String>,
    
    /// Branch length to parent
    /// In rooted trees, edge length is an attribute of the child node.
    pub length: Option<f64>,
    
    /// Structured properties (e.g., NHX tags like [&&NHX:S=human])
    /// Using BTreeMap ensures deterministic output order.
    pub properties: Option<BTreeMap<String, String>>,
    
    /// Soft deletion flag.
    /// If true, this node is considered removed.
    /// Use Tree::compact() to permanently remove deleted nodes and reclaim memory.
    pub deleted: bool,
}

impl Node {
    /// Create a new empty node with a specific ID
    pub fn new(id: NodeId) -> Self {
        Self {
            id,
            parent: None,
            children: Vec::new(),
            name: None,
            length: None,
            properties: None,
            deleted: false,
        }
    }

    /// Set the name of the node
    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = Some(name.into());
    }

    /// Set the name of the node (builder pattern)
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Set the branch length
    pub fn with_length(mut self, length: f64) -> Self {
        self.length = Some(length);
        self
    }

    /// Add a property (key-value pair)
    pub fn add_property(&mut self, key: impl Into<String>, value: impl Into<String>) {
        if self.properties.is_none() {
            self.properties = Some(BTreeMap::new());
        }
        self.properties.as_mut().unwrap().insert(key.into(), value.into());
    }

    /// Check if the node is a leaf (no children)
    /// Note: This ignores soft-deleted children if checked externally,
    /// but strictly speaking a node is a leaf if `children` is empty.
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
}
