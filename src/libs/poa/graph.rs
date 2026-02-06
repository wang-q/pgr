use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::algo::toposort;
use super::align::Alignment;

#[derive(Clone, Debug, PartialEq)]
pub struct NodeData {
    pub base: u8,
    pub aligned_to: Option<NodeIndex>,
    pub weight: u32,
}

impl NodeData {
    pub fn new(base: u8) -> Self {
        Self { base, aligned_to: None, weight: 0 }
    }
}

pub type EdgeData = u32; // Weight

pub struct PoaGraph {
    pub graph: DiGraph<NodeData, EdgeData>,
}

impl PoaGraph {
    pub fn new() -> Self {
        Self {
            graph: DiGraph::new(),
        }
    }

    pub fn add_node(&mut self, base: u8) -> NodeIndex {
        self.graph.add_node(NodeData::new(base))
    }

    pub fn add_edge(&mut self, from: NodeIndex, to: NodeIndex, weight: u32) {
        if let Some(edge_idx) = self.graph.find_edge(from, to) {
            let w = self.graph.edge_weight_mut(edge_idx).unwrap();
            *w += weight;
        } else {
            self.graph.add_edge(from, to, weight);
        }
    }

    pub fn topological_sort(&self) -> Vec<NodeIndex> {
        toposort(&self.graph, None).unwrap_or_else(|_| panic!("Graph has cycles!"))
    }
    
    pub fn num_nodes(&self) -> usize {
        self.graph.node_count()
    }

    /// Adds an alignment to the graph, updating weights and adding new nodes/edges as needed.
    pub fn add_alignment(&mut self, alignment: &Alignment, sequence: &[u8]) {
        let mut prev_node: Option<NodeIndex> = None;
        
        for step in &alignment.path {
            match step {
                (Some(seq_idx), Some(graph_idx)) => {
                    let base = sequence[*seq_idx];
                    let node_base = self.graph[*graph_idx].base;
                    
                    let target_node = if base == node_base {
                        // Match: reuse node
                        self.graph[*graph_idx].weight += 1;
                        *graph_idx
                    } else {
                        // Mismatch: create new node, link to aligned graph node
                        let mut data = NodeData::new(base);
                        data.aligned_to = Some(*graph_idx);
                        data.weight = 1;
                        self.graph.add_node(data)
                    };
                    
                    if let Some(p) = prev_node {
                        self.add_edge(p, target_node, 1);
                    }
                    prev_node = Some(target_node);
                },
                (Some(seq_idx), None) => {
                    // Insertion
                    let base = sequence[*seq_idx];
                    let mut data = NodeData::new(base);
                    data.weight = 1;
                    let new_node = self.graph.add_node(data);
                    
                    if let Some(p) = prev_node {
                        self.add_edge(p, new_node, 1);
                    }
                    prev_node = Some(new_node);
                },
                (None, Some(_)) => {
                    // Deletion: skip
                    // prev_node remains same
                },
                (None, None) => {},
            }
        }
    }
}

impl Default for PoaGraph {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graph_creation() {
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        let n2 = graph.add_node(b'C');
        graph.add_edge(n1, n2, 1);
        
        assert_eq!(graph.graph.node_count(), 2);
        assert_eq!(graph.graph.edge_count(), 1);
    }

    #[test]
    fn test_add_alignment_linear() {
        // Test adding a linear sequence (all insertions relative to empty graph)
        let mut graph = PoaGraph::new();
        let sequence = b"ACGT";
        let alignment = Alignment {
            score: 0,
            path: vec![
                (Some(0), None),
                (Some(1), None),
                (Some(2), None),
                (Some(3), None),
            ],
        };
        
        graph.add_alignment(&alignment, sequence);
        
        assert_eq!(graph.num_nodes(), 4);
        assert_eq!(graph.graph.edge_count(), 3); // A->C, C->G, G->T
        
        let sorted = graph.topological_sort();
        assert_eq!(graph.graph[sorted[0]].base, b'A');
        assert_eq!(graph.graph[sorted[3]].base, b'T');
    }

    #[test]
    fn test_add_alignment_match() {
        // First seq: A C G
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        let n2 = graph.add_node(b'C');
        let n3 = graph.add_node(b'G');
        graph.add_edge(n1, n2, 1);
        graph.add_edge(n2, n3, 1);
        
        // Second seq: A C G (Identical)
        let sequence = b"ACG";
        let alignment = Alignment {
            score: 0,
            path: vec![
                (Some(0), Some(n1)),
                (Some(1), Some(n2)),
                (Some(2), Some(n3)),
            ],
        };
        
        graph.add_alignment(&alignment, sequence);
        
        assert_eq!(graph.num_nodes(), 3); // No new nodes
        // Weights should increase
        let e1 = graph.graph.find_edge(n1, n2).unwrap();
        assert_eq!(*graph.graph.edge_weight(e1).unwrap(), 2);
    }

    #[test]
    fn test_add_alignment_mismatch() {
        // First seq: A
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        
        // Second seq: C (Mismatch aligned to A)
        let sequence = b"C";
        let alignment = Alignment {
            score: 0,
            path: vec![
                (Some(0), Some(n1)),
            ],
        };
        
        graph.add_alignment(&alignment, sequence);
        
        assert_eq!(graph.num_nodes(), 2); // A and C
        let n2 = graph.graph.node_indices().find(|&n| n != n1).unwrap();
        assert_eq!(graph.graph[n2].base, b'C');
        assert_eq!(graph.graph[n2].aligned_to, Some(n1));
    }
}
