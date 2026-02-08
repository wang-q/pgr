use super::graph::PoaGraph;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use petgraph::Direction;
use std::collections::HashMap;

/// Generates a consensus sequence from the POA graph.
/// Uses a heaviest path algorithm (finding the path with maximum total weight).
/// Score[u] = NodeWeight[u] + max(Score[v] + EdgeWeight(v, u)) for all predecessors v.
pub fn generate_consensus(graph: &PoaGraph) -> Vec<u8> {
    let sorted_nodes = graph.topological_sort();

    if sorted_nodes.is_empty() {
        return Vec::new();
    }

    // Initialize scores
    // Spoa strategy: "Heaviest Bundle"
    // Score[u] = max_edge_weight(v->u) + Score[v]
    // Prioritizes heavy edges first, then predecessor score.
    // Score initialized to -1 (implicit in Spoa).
    // Here we use i64 to allow -1, though weights are u32.
    let mut scores: HashMap<NodeIndex, i64> = HashMap::new();
    let mut predecessors: HashMap<NodeIndex, NodeIndex> = HashMap::new();

    // Since we process in topological order, predecessors are already processed.
    for &node_idx in &sorted_nodes {
        let mut best_edge_weight = -1;
        let mut best_prev = None;

        // Collect and sort incoming edges to ensure deterministic behavior
        // Sort by source node index to consistently process S1 then S2
        let mut edges: Vec<_> = graph
            .graph
            .edges_directed(node_idx, Direction::Incoming)
            .collect();
        edges.sort_by_key(|e| e.source().index());

        // Iterate over incoming edges
        for edge_ref in edges {
            let prev_node = edge_ref.source();
            let edge_weight = *edge_ref.weight() as i64;

            // Get prev score, default to -1 (start node score)
            let prev_score = *scores.get(&prev_node).unwrap_or(&-1);

            // Spoa Logic:
            // if (scores[curr] < weight) || (scores[curr] == weight && scores[prev] <= scores[new_prev])
            // Here scores[curr] tracks the *best edge weight* seen so far for this node loop

            if best_edge_weight < edge_weight {
                best_edge_weight = edge_weight;
                best_prev = Some(prev_node);
            } else if best_edge_weight == edge_weight {
                // Tie-breaker: Check predecessor total scores
                // Spoa uses <= to swap. So it prefers the NEW predecessor if its score is >= OLD best.
                // We need to compare prev_score vs score_of_best_prev
                if let Some(curr_best_prev) = best_prev {
                    let best_prev_score = *scores.get(&curr_best_prev).unwrap_or(&-1);
                    if best_prev_score <= prev_score {
                        best_prev = Some(prev_node);
                    }
                } else {
                    best_prev = Some(prev_node);
                }
            }
        }

        // Calculate total score for this node
        let mut total_score = -1;
        if let Some(prev) = best_prev {
            predecessors.insert(node_idx, prev);
            let prev_score = *scores.get(&prev).unwrap_or(&-1);
            total_score = best_edge_weight + prev_score;
        }

        scores.insert(node_idx, total_score);
    }

    // Find the node with the highest score
    // Iterate sorted_nodes to ensure topological order
    let mut max_score = -1; // Spoa init max to nullptr/score -1
    let mut end_node = None;

    for &node in &sorted_nodes {
        if let Some(&score) = scores.get(&node) {
            if score > max_score {
                max_score = score;
                end_node = Some(node);
            } else if score == max_score {
                // Tie-breaker: Prefer lower node index (usually S1/Backbone)
                if let Some(curr) = end_node {
                    if node.index() < curr.index() {
                        end_node = Some(node);
                    }
                }
            }
        }
    }

    // If all scores are -1 (single node graph?), end_node might be None if we init max_score = -1.
    // But if graph has nodes, scores will be -1.
    // Spoa: if (!max || scores[max] < scores[it])
    // If scores[it] is -1.
    // If max is null, it sets max = it.
    // So it picks the first node.
    if end_node.is_none() && !sorted_nodes.is_empty() {
        end_node = Some(sorted_nodes[0]);
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
    use crate::libs::poa::align::{
        AlignmentEngine, AlignmentParams, AlignmentType, ScalarAlignmentEngine,
    };

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
