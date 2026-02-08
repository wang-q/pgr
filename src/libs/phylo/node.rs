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
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::node::Node;
    /// let node = Node::new(1);
    /// assert_eq!(node.id, 1);
    /// assert!(node.children.is_empty());
    /// assert!(node.name.is_none());
    /// ```
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
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::node::Node;
    /// let mut node = Node::new(1);
    /// node.set_name("Node1");
    /// assert_eq!(node.name, Some("Node1".to_string()));
    /// ```
    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = Some(name.into());
    }

    /// Set the name of the node (builder pattern)
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::node::Node;
    /// let node = Node::new(1).with_name("Node1");
    /// assert_eq!(node.name, Some("Node1".to_string()));
    /// ```
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Set the branch length
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::node::Node;
    /// let node = Node::new(1).with_length(0.5);
    /// assert_eq!(node.length, Some(0.5));
    /// ```
    pub fn with_length(mut self, length: f64) -> Self {
        self.length = Some(length);
        self
    }

    /// Add a property (key-value pair)
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::node::Node;
    /// let mut node = Node::new(1);
    /// node.add_property("color", "blue");
    /// assert_eq!(node.get_property("color"), Some(&"blue".to_string()));
    /// ```
    pub fn add_property(&mut self, key: impl Into<String>, value: impl Into<String>) {
        if self.properties.is_none() {
            self.properties = Some(BTreeMap::new());
        }
        self.properties.as_mut().unwrap().insert(key.into(), value.into());
    }

    /// Add properties from a string.
    /// Supports single "key=value" or multiple "key1=value1:key2=value2" (NHX style).
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::node::Node;
    /// let mut node = Node::new(0);
    /// node.add_property_from_str("color=red");
    /// assert_eq!(node.get_property("color"), Some(&"red".to_string()));
    ///
    /// // Multiple properties separated by colon
    /// node.add_property_from_str("S=Homo sapiens:T=9606");
    /// assert_eq!(node.get_property("S"), Some(&"Homo sapiens".to_string()));
    /// assert_eq!(node.get_property("T"), Some(&"9606".to_string()));
    /// ```
    pub fn add_property_from_str(&mut self, props_str: &str) {
        for part in props_str.split(':') {
            if let Some((key, value)) = part.split_once('=') {
                self.add_property(key, value);
            }
        }
    }

    /// Get the value of a property by key.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::node::Node;
    /// let mut node = Node::new(0);
    /// node.add_property("T", "9606");
    /// node.add_property("S", "Homo sapiens");
    /// 
    /// assert_eq!(node.get_property("T"), Some(&"9606".to_string()));
    /// assert_eq!(node.get_property("S"), Some(&"Homo sapiens".to_string()));
    /// assert_eq!(node.get_property("Missing"), None);
    /// ```
    pub fn get_property(&self, key: &str) -> Option<&String> {
        self.properties.as_ref().and_then(|p| p.get(key))
    }

    /// Check if the node is a leaf (no children)
    /// Note: This ignores soft-deleted children if checked externally,
    /// but strictly speaking a node is a leaf if `children` is empty.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::node::Node;
    /// let mut node = Node::new(1);
    /// assert!(node.is_leaf());
    /// 
    /// // Manually adding a child ID (simulating tree operation)
    /// node.children.push(2);
    /// assert!(!node.is_leaf());
    /// ```
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
}
