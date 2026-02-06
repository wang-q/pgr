use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::NodeIndexable;
use super::align::Alignment;

#[derive(Clone, Debug, PartialEq)]
pub struct NodeData {
    pub base: u8,
    pub aligned_nodes: Vec<NodeIndex>,
    pub weight: u32,
}

impl NodeData {
    pub fn new(base: u8) -> Self {
        Self { base, aligned_nodes: Vec::new(), weight: 0 }
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
        let mut sorted_nodes = Vec::with_capacity(self.graph.node_count());
        let mut marks = vec![0u8; self.graph.node_bound()]; // 0: Unvisited, 1: Visiting, 2: Visited
        let mut ignored = vec![false; self.graph.node_bound()];
        let mut stack = Vec::new();

        // Iterate over all nodes in index order (to match Spoa's iteration order 0..N)
        // Note: node_indices() in petgraph iterates indices.
        for node in self.graph.node_indices() {
            if marks[node.index()] != 0 {
                continue;
            }
            stack.push(node);

            while let Some(&curr) = stack.last() {
                let curr_idx = curr.index();
                
                // If already processed (visited), pop and continue
                if marks[curr_idx] == 2 {
                    stack.pop();
                    continue;
                }

                let mut is_valid = true;

                // Check incoming edges (dependencies)
                for neighbor in self.graph.neighbors_directed(curr, petgraph::Direction::Incoming) {
                    if marks[neighbor.index()] != 2 {
                        stack.push(neighbor);
                        is_valid = false;
                    }
                }

                // Check aligned nodes
                // If this node is not ignored, we must ensure its aligned clique is also ready
                if !ignored[curr_idx] {
                    let node_data = &self.graph[curr];
                    for &aligned in &node_data.aligned_nodes {
                        if marks[aligned.index()] != 2 {
                            stack.push(aligned);
                            ignored[aligned.index()] = true; // Mark aligned node as ignored (merged into clique leader)
                            is_valid = false;
                        }
                    }
                }

                if is_valid {
                    marks[curr_idx] = 2;
                    
                    if !ignored[curr_idx] {
                        sorted_nodes.push(curr);
                        // Add aligned nodes immediately to keep clique together
                        let node_data = &self.graph[curr];
                        for &aligned in &node_data.aligned_nodes {
                            // Only add if not already added? 
                            // Spoa logic: marks[aligned] becomes 2. 
                            // But we just checked marks[aligned] == 2 above?
                            // Wait, if marks[aligned] == 2, it means it was visited.
                            // But if it was ignored, it wasn't added to sorted_nodes yet?
                            // Yes, ignored nodes are NOT added in their own loop iteration.
                            // They are added here.
                            sorted_nodes.push(aligned);
                        }
                    }
                    stack.pop();
                } else {
                    // Dependencies not met. Mark as visiting.
                    if marks[curr_idx] == 1 {
                        // Cycle detected?
                        // If we are already visiting this node and we found an unsatisfied dependency,
                        // that dependency should have been pushed.
                        // If the dependency is THIS node (self-loop) or a loop back to THIS node,
                        // we would encounter a node with mark 1 in the dependency check.
                        // My dependency check loop doesn't check for mark 1 explicitly to panic, 
                        // but Spoa asserts `marks[curr] != 1`.
                    }
                    marks[curr_idx] = 1;
                }
            }
        }
        
        sorted_nodes
    }
    
    pub fn num_nodes(&self) -> usize {
        self.graph.node_count()
    }

    /// Adds an alignment to the graph, updating weights and adding new nodes/edges as needed.
    /// Returns the path of node indices corresponding to the sequence.
    pub fn add_alignment(&mut self, alignment: &Alignment, sequence: &[u8]) -> Vec<NodeIndex> {
        let mut prev_node: Option<NodeIndex> = None;
        let mut current_seq_idx = 0;
        let mut sequence_path = Vec::with_capacity(sequence.len());
        
        for step in &alignment.path {
            // Get sequence index if present
            let seq_idx = match step {
                (Some(idx), _) => *idx,
                _ => {
                    // Deletion (None, Some) or (None, None)
                    // Does not consume sequence, just skip
                    continue;
                }
            };

            // Fill unaligned sequence bases (prefix or gaps)
            while current_seq_idx < seq_idx {
                let base = sequence[current_seq_idx];
                let mut data = NodeData::new(base);
                data.weight = 1;
                let new_node = self.graph.add_node(data);
                sequence_path.push(new_node);
                
                if let Some(p) = prev_node {
                    self.add_edge(p, new_node, 1);
                }
                prev_node = Some(new_node);
                current_seq_idx += 1;
            }

            match step {
                (Some(idx), Some(graph_idx)) => {
                    let base = sequence[*idx];
                    let graph_node_base = self.graph[*graph_idx].base;
                    
                    let target_node_idx = if graph_node_base == base {
                        *graph_idx
                    } else {
                        // Check aligned nodes
                        let mut found = None;
                        for &aligned_idx in &self.graph[*graph_idx].aligned_nodes {
                            if self.graph[aligned_idx].base == base {
                                found = Some(aligned_idx);
                                break;
                            }
                        }
                        
                        if let Some(idx) = found {
                            idx
                        } else {
                            // Create new node
                            let mut data = NodeData::new(base);
                            data.weight = 0; // Will be incremented below
                            let new_node = self.graph.add_node(data);
                            
                            // Update cliques
                            let mut clique = self.graph[*graph_idx].aligned_nodes.clone();
                            clique.push(*graph_idx);
                            
                            for &peer in &clique {
                                self.graph[peer].aligned_nodes.push(new_node);
                            }
                            self.graph[new_node].aligned_nodes = clique;
                            
                            new_node
                        }
                    };
                    
                    sequence_path.push(target_node_idx);
                    
                    // Increment weight
                    self.graph[target_node_idx].weight += 1;
                    
                    if let Some(p) = prev_node {
                        self.add_edge(p, target_node_idx, 1);
                    }
                    prev_node = Some(target_node_idx);
                },
                (Some(idx), None) => {
                    // Insertion
                    let base = sequence[*idx];
                    let mut data = NodeData::new(base);
                    data.weight = 1;
                    let new_node = self.graph.add_node(data);
                    sequence_path.push(new_node);
                    
                    if let Some(p) = prev_node {
                        self.add_edge(p, new_node, 1);
                    }
                    prev_node = Some(new_node);
                },
                _ => {} // Handled above
            }
            
            current_seq_idx = seq_idx + 1;
        }

        // Fill remaining suffix
        while current_seq_idx < sequence.len() {
             let base = sequence[current_seq_idx];
             let mut data = NodeData::new(base);
             data.weight = 1;
             let new_node = self.graph.add_node(data);
             sequence_path.push(new_node);
             
             if let Some(p) = prev_node {
                 self.add_edge(p, new_node, 1);
             }
             prev_node = Some(new_node);
             current_seq_idx += 1;
        }

        sequence_path
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
        assert!(graph.graph[n2].aligned_nodes.contains(&n1));
        assert!(graph.graph[n1].aligned_nodes.contains(&n2));
    }
}
