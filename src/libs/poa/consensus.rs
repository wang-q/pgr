use super::graph::PoaGraph;
use petgraph::graph::NodeIndex;
use petgraph::Direction;
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

/// Generates a consensus sequence from the POA graph.
/// Uses a heaviest path algorithm (finding the path with maximum total weight).
/// Score[u] = NodeWeight[u] + max(Score[v] + EdgeWeight(v, u)) for all predecessors v.
pub fn generate_consensus(graph: &PoaGraph) -> Vec<u8> {
    let sorted_nodes = graph.topological_sort();
    
    if sorted_nodes.is_empty() {
        return Vec::new();
    }

    let mut scores: HashMap<NodeIndex, u32> = HashMap::new();
    let mut predecessors: HashMap<NodeIndex, NodeIndex> = HashMap::new();

    // Initialize scores
    // Since we process in topological order, predecessors are already processed.
    for &node_idx in &sorted_nodes {
        let node_weight = graph.graph[node_idx].weight;
        let mut max_prev_score = 0;
        let mut best_prev = None;

        // Collect and sort incoming edges to ensure deterministic behavior
        // Sort by source node index to consistently process S1 then S2
        let mut edges: Vec<_> = graph.graph.edges_directed(node_idx, Direction::Incoming).collect();
        edges.sort_by_key(|e| e.source().index());

        // Iterate over incoming edges
        for edge_ref in edges {
            let prev_node = edge_ref.source();
            let edge_weight = *edge_ref.weight();
            
            // If prev_node is not in scores, it might be unreachable or bug in sort?
            // Topo sort guarantees we visited prev_node unless graph has cycles (handled by panic in sort).
            if let Some(&prev_score) = scores.get(&prev_node) {
                let current_score = prev_score + edge_weight;
                // Use >= to favor later edges (usually S2/S3...) in ties
                // This matches spoa behavior where latest sequence wins in internal bubbles
                if current_score >= max_prev_score {
                    max_prev_score = current_score;
                    best_prev = Some(prev_node);
                }
            }
        }

        let total_score = node_weight + max_prev_score;
        scores.insert(node_idx, total_score);
        if let Some(prev) = best_prev {
            predecessors.insert(node_idx, prev);
        }
    }

    // Find the node with the highest score
    // Iterate sorted_nodes to ensure topological order
    let mut max_score = 0;
    let mut end_node = None;

    for &node in &sorted_nodes {
        if let Some(&score) = scores.get(&node) {
            if score > max_score {
                max_score = score;
                end_node = Some(node);
            } else if score == max_score {
                // Tie-breaker: Prefer lower node index (usually S1/Backbone)
                // This handles cases like Case 1 where S1 ends with GG and S2 ends with TT
                // and we want to preserve the S1 backbone ending.
                if let Some(curr) = end_node {
                    if node.index() < curr.index() {
                        end_node = Some(node);
                    }
                }
            }
        }
    }

    // Backtrack
    let mut consensus = Vec::new();
    if let Some(mut curr) = end_node {
        loop {
            consensus.push(graph.graph[curr].base);
            if let Some(&prev) = predecessors.get(&curr) {
                curr = prev;
            } else {
                break;
            }
        }
    }

    consensus.reverse();
    consensus
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::poa::align::{ScalarAlignmentEngine, AlignmentParams, AlignmentType, AlignmentEngine};

    #[test]
    fn test_consensus_simple() {
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        let n2 = graph.add_node(b'C');
        let n3 = graph.add_node(b'G');
        graph.add_edge(n1, n2, 10);
        graph.add_edge(n2, n3, 10);
        
        graph.graph[n1].weight = 10;
        graph.graph[n2].weight = 10;
        graph.graph[n3].weight = 10;

        let consensus = generate_consensus(&graph);
        assert_eq!(consensus, b"ACG");
    }

    #[test]
    fn test_consensus_branching() {
        // A -> C -> G (weight 10 path)
        // A -> T -> G (weight 5 path)
        let mut graph = PoaGraph::new();
        let n_a = graph.add_node(b'A');
        let n_c = graph.add_node(b'C');
        let n_t = graph.add_node(b'T');
        let n_g = graph.add_node(b'G');
        
        graph.graph[n_a].weight = 10;
        graph.graph[n_c].weight = 10;
        graph.graph[n_t].weight = 5;
        graph.graph[n_g].weight = 10;

        graph.add_edge(n_a, n_c, 10);
        graph.add_edge(n_c, n_g, 10);
        
        graph.add_edge(n_a, n_t, 5);
        graph.add_edge(n_t, n_g, 5);

        let consensus = generate_consensus(&graph);
        assert_eq!(consensus, b"ACG");
    }

    #[test]
    fn test_consensus_with_alignment() {
        let mut graph = PoaGraph::new();
        let engine = ScalarAlignmentEngine::new(AlignmentParams::default(), AlignmentType::Global);
        
        // Seq 1: ACGT
        let s1 = b"ACGT";
        let a1 = engine.align(s1, &graph);
        graph.add_alignment(&a1, s1);
        
        // Seq 2: ACGC
        let s2 = b"ACGC";
        let a2 = engine.align(s2, &graph);
        graph.add_alignment(&a2, s2);
        
        // Seq 3: ACGT
        let s3 = b"ACGT";
        let a3 = engine.align(s3, &graph);
        graph.add_alignment(&a3, s3);
        
        // Consensus should be ACGT (2 votes for T, 1 for C)
        let consensus = generate_consensus(&graph);
        assert_eq!(consensus, b"ACGT");
    }
}
