use super::graph::PoaGraph;
use petgraph::graph::NodeIndex;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlignmentType {
    Global,     // Needleman-Wunsch
    Local,      // Smith-Waterman
    SemiGlobal, // Overlap (Sequence fully aligned, Graph start/end free)
}

#[derive(Debug, Clone)]
pub struct AlignmentParams {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for AlignmentParams {
    fn default() -> Self {
        Self {
            match_score: 5,
            mismatch_score: -4,
            gap_open: -8,
            gap_extend: -6,
        }
    }
}

/// Represents an alignment between a sequence and a graph.
#[derive(Debug, Clone, Default)]
pub struct Alignment {
    pub score: i32,
    pub path: Vec<(Option<usize>, Option<NodeIndex>)>, 
}

pub trait AlignmentEngine {
    fn align(&self, sequence: &[u8], graph: &PoaGraph) -> Alignment;
}

pub struct ScalarAlignmentEngine {
    pub params: AlignmentParams,
    pub align_type: AlignmentType,
}

impl ScalarAlignmentEngine {
    pub fn new(params: AlignmentParams, align_type: AlignmentType) -> Self {
        Self { params, align_type }
    }
}

impl AlignmentEngine for ScalarAlignmentEngine {
    fn align(&self, sequence: &[u8], graph: &PoaGraph) -> Alignment {
        let sorted_nodes = graph.topological_sort();
        let n_nodes = sorted_nodes.len();
        let n_seq = sequence.len();
        
        if n_nodes == 0 {
             let mut path = Vec::new();
             for i in 0..n_seq {
                 path.push((Some(i), None));
             }
             return Alignment { score: 0, path };
        }

        let node_map: HashMap<NodeIndex, usize> = sorted_nodes
            .iter()
            .enumerate()
            .map(|(i, &n)| (n, i))
            .collect();

        let neg_inf = -1_000_000_000;
        
        // Matrices: [node_idx_linear][seq_idx]
        let mut m = vec![vec![neg_inf; n_seq + 1]; n_nodes];
        let mut e = vec![vec![neg_inf; n_seq + 1]; n_nodes];
        let mut f = vec![vec![neg_inf; n_seq + 1]; n_nodes];

        let is_local = self.align_type == AlignmentType::Local;
        let is_semi = self.align_type == AlignmentType::SemiGlobal;

        // Initialization
        for (i, &node_idx) in sorted_nodes.iter().enumerate() {
            let preds: Vec<NodeIndex> = graph.graph.neighbors_directed(node_idx, petgraph::Direction::Incoming).collect();
            let is_start_node = preds.is_empty();
            
            // 1. Initialize Column 0 (Sequence Empty)
            if is_local || is_semi {
                // Free start in graph: 0 cost to reach any node with empty sequence
                f[i][0] = 0;
                m[i][0] = neg_inf;
                e[i][0] = neg_inf;
            } else {
                if is_start_node {
                    f[i][0] = self.params.gap_open; 
                    m[i][0] = neg_inf;
                    e[i][0] = neg_inf;
                } else {
                    let mut max_prev = neg_inf;
                    for &pred in &preds {
                        let u = node_map[&pred];
                        // Spoa: penalty = max(penalty, F[pred])
                        // F[i] = penalty + e
                        if f[u][0] > neg_inf {
                             max_prev = max_prev.max(f[u][0]);
                        }
                    }
                    if max_prev > neg_inf {
                        f[i][0] = max_prev + self.params.gap_extend;
                    } else {
                        f[i][0] = neg_inf;
                    }
                    m[i][0] = neg_inf;
                    e[i][0] = neg_inf;
                }
            }
            
            // 2. Fill rest of columns
            let node_base = graph.graph[node_idx].base;
            
            for j in 1..=n_seq {
                 let seq_base = sequence[j-1];
                 let match_score = if seq_base == node_base { self.params.match_score } else { self.params.mismatch_score };
                 
                 // E[i][j]: Insertion
                 // Derived from M (gap open), E (gap extend), or F (gap open)
                 let from_m = if m[i][j-1] > neg_inf { m[i][j-1] + self.params.gap_open } else { neg_inf };
                 let from_e = if e[i][j-1] > neg_inf { e[i][j-1] + self.params.gap_extend } else { neg_inf };
                 let from_f = if f[i][j-1] > neg_inf { f[i][j-1] + self.params.gap_open } else { neg_inf };
                 
                 let mut max_e = from_m.max(from_e).max(from_f);
                 if is_local && max_e < 0 { max_e = neg_inf; } 
                 if max_e < neg_inf / 2 { max_e = neg_inf; }
                 e[i][j] = max_e;
                 
                 // M[i][j]: Match/Mismatch
                 let mut max_m = neg_inf;
                 
                 if is_start_node {
                     if j == 1 {
                         // Transition from virtual root
                         max_m = match_score;
                     } else {
                         // Gap from virtual root
                         // We have j-1 bases before current one.
                         // They are all insertions.
                         // Cost = Open + (cnt-1)*Extend.
                         // cnt = j-1.
                         // Cost = Open + (j-2)*Extend.
                         let ins_score = self.params.gap_open + (j as i32 - 2) * self.params.gap_extend;
                         max_m = ins_score + match_score;
                     }
                 } else {
                     for &pred in &preds {
                         let u = node_map[&pred];
                         let src = m[u][j-1].max(e[u][j-1]).max(f[u][j-1]);
                         if src > neg_inf {
                            max_m = max_m.max(src + match_score);
                         }
                     }
                     // For Semi/Local: Free start allows starting here from virtual root
                     if (is_local || is_semi) && j == 1 {
                         max_m = max_m.max(match_score);
                     }
                 }
                 
                 if is_local && max_m < 0 { max_m = 0; }
                 m[i][j] = max_m;
                 
                 // F[i][j]: Deletion
                 let mut max_f = neg_inf;
                 if !is_start_node {
                     for &pred in &preds {
                         let u = node_map[&pred];
                         let from_m = if m[u][j] > neg_inf { m[u][j] + self.params.gap_open } else { neg_inf };
                         let from_f = if f[u][j] > neg_inf { f[u][j] + self.params.gap_extend } else { neg_inf };
                         let from_e = if e[u][j] > neg_inf { e[u][j] + self.params.gap_open } else { neg_inf };
                         
                         max_f = max_f.max(from_m).max(from_f).max(from_e);
                     }
                 }
                 // For Semi/Local, F[i][j] coming from virtual root (skipping nodes) is handled by F[i][0] init and propagation?
                 // No, F[i][0]=0. F[i][j] means we have consumed seq[0..j].
                 // If we delete node i, we move from (pred, j) -> (i, j).
                 // So F[i][j] depends on predecessors at j.
                 
                 if is_local && max_f < 0 { max_f = neg_inf; }
                 if max_f < neg_inf / 2 { max_f = neg_inf; }
                 f[i][j] = max_f;
            }
        }
        
        // Find best end score
        let mut best_score = neg_inf;
        let mut best_node_idx = 0;
        let mut best_col = n_seq;
        let mut best_state = 0; 
        
        if is_local {
            // Check all cells
            best_score = 0; 
            for i in 0..n_nodes {
                for j in 1..=n_seq {
                    let score = m[i][j].max(e[i][j]).max(f[i][j]);
                    if score >= best_score { // Use >= to pick last occurrence? or >?
                        best_score = score;
                        best_node_idx = i;
                        best_col = j;
                        if score == m[i][j] { best_state = 0; }
                        else if score == e[i][j] { best_state = 1; }
                        else { best_state = 2; }
                    }
                }
            }
        } else if is_semi {
            // Check all nodes at last column
            for i in 0..n_nodes {
                let score = m[i][n_seq].max(e[i][n_seq]).max(f[i][n_seq]);
                if score > best_score {
                    best_score = score;
                    best_node_idx = i;
                    best_col = n_seq;
                    if score == m[i][n_seq] { best_state = 0; }
                    else if score == e[i][n_seq] { best_state = 1; }
                    else { best_state = 2; }
                }
            }
        } else {
            // Global: Check all nodes at last column (Free end in graph)
            // This allows the sequence to end before the graph ends without penalty (Semi-Global in Target)
            // which is consistent with Spoa behavior for consensus.
            for i in 0..n_nodes {
                let score = m[i][n_seq].max(e[i][n_seq]).max(f[i][n_seq]);
                if score > best_score {
                    best_score = score;
                    best_node_idx = i;
                    best_col = n_seq;
                    if score == m[i][n_seq] { best_state = 0; }
                    else if score == e[i][n_seq] { best_state = 1; }
                    else { best_state = 2; }
                }
            }
        }
        
        // Backtracking
        let mut path = Vec::new();
        let mut curr_i = best_node_idx;
        let mut curr_j = best_col;
        let mut curr_state = best_state;
        
        while curr_j > 0 || (curr_i > 0) { 
             let node_idx = sorted_nodes[curr_i];
             let preds: Vec<NodeIndex> = graph.graph.neighbors_directed(node_idx, petgraph::Direction::Incoming).collect();
             let is_start = preds.is_empty();
             
             // Stop conditions
             if is_local && best_score == 0 { break; } // If score 0, stop
             if is_local {
                 let s = match curr_state { 0 => m[curr_i][curr_j], 1 => e[curr_i][curr_j], _ => f[curr_i][curr_j] };
                 if s <= 0 { break; }
             }
             if is_semi && curr_j == 0 { break; } // Reached start of sequence
             
             if curr_j == 0 && is_start { break; }
             
             match curr_state {
                 0 => { // M
                     let match_score = if curr_j > 0 {
                         if sequence[curr_j-1] == graph.graph[node_idx].base { self.params.match_score } else { self.params.mismatch_score }
                     } else { 0 };

                     // If j=1 and Local/Semi, we could have started here
                     if (is_local || is_semi) && curr_j == 1 {
                         // Check if we started here (score == match_score)
                         if m[curr_i][curr_j] == match_score {
                             path.push((Some(curr_j-1), Some(node_idx)));
                             curr_j -= 1;
                             break;
                         }
                     }
                     
                     if is_start {
                         if curr_j > 0 {
                            path.push((Some(curr_j-1), Some(node_idx)));
                            curr_j -= 1;
                         }
                         break; 
                     } else {
                         let mut found = false;
                         for &pred in &preds {
                             let u = node_map[&pred];
                             let target = m[curr_i][curr_j] - match_score;
                             
                             if m[u][curr_j-1] == target {
                                 path.push((Some(curr_j-1), Some(node_idx)));
                                 curr_i = u; curr_j -= 1; curr_state = 0; found = true; break;
                             }
                             if e[u][curr_j-1] == target {
                                 path.push((Some(curr_j-1), Some(node_idx)));
                                 curr_i = u; curr_j -= 1; curr_state = 1; found = true; break;
                             }
                             if f[u][curr_j-1] == target {
                                 path.push((Some(curr_j-1), Some(node_idx)));
                                 curr_i = u; curr_j -= 1; curr_state = 2; found = true; break;
                             }
                         }
                         if !found { 
                             // Could be start of a branch from virtual source?
                             break;
                         }
                     }
                 },
                 1 => { // E
                     let target = e[curr_i][curr_j];
                     let score_e = e[curr_i][curr_j-1] + self.params.gap_extend;
                     
                     path.push((Some(curr_j-1), None));
                     
                     // Spoa checks H[i][j] == H[i][j-1] + g_ (gap open) for insertion start
                     // Spoa checks H[i][j] == E[i][j-1] + e_ (gap extend) for insertion extend
                     if target == score_e {
                         curr_j -= 1; curr_state = 1;
                     } else {
                            // Transition from M or F
                            let score_m = m[curr_i][curr_j-1] + self.params.gap_open;
                            if target == score_m {
                                curr_j -= 1; curr_state = 0;
                            } else {
                                // Must be F
                                // Verify for correctness/safety
                                // let score_f = f[curr_i][curr_j-1] + self.params.gap_open;
                                // if target == score_f { ... }
                                curr_j -= 1; curr_state = 2;
                            }
                        }
                 },
                 2 => { // F
                     let mut found = false;
                      for &pred in &preds {
                         let u = node_map[&pred];
                         let target = f[curr_i][curr_j];
                         if f[u][curr_j] + self.params.gap_extend == target {
                              path.push((None, Some(node_idx)));
                              curr_i = u; curr_state = 2; found = true; break;
                         }
                         if m[u][curr_j] + self.params.gap_open == target {
                              path.push((None, Some(node_idx)));
                              curr_i = u; curr_state = 0; found = true; break;
                         }
                         if e[u][curr_j] + self.params.gap_open == target {
                              path.push((None, Some(node_idx)));
                              curr_i = u; curr_state = 1; found = true; break;
                         }
                      }
                      if !found { 
                          if is_start {
                              path.push((None, Some(node_idx)));
                              break;
                          }
                          break;
                      }
                 }
                 _ => break,
             }
        }
        
        // Fill remaining sequence if any (for Global)
        // For Local/Semi, we stop.
        if !is_local && !is_semi {
            while curr_j > 0 {
                 path.push((Some(curr_j-1), None));
                 curr_j -= 1;
            }
        }

        path.reverse();
        
        Alignment {
            score: best_score,
            path,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::poa::graph::PoaGraph;

    #[test]
    fn test_align_exact_match() {
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        let n2 = graph.add_node(b'C');
        let n3 = graph.add_node(b'G');
        let n4 = graph.add_node(b'T');
        graph.add_edge(n1, n2, 1);
        graph.add_edge(n2, n3, 1);
        graph.add_edge(n3, n4, 1);

        let engine = ScalarAlignmentEngine::new(AlignmentParams::default(), AlignmentType::Global);
        let seq = b"ACGT";
        let alignment = engine.align(seq, &graph);

        assert_eq!(alignment.path.len(), 4);
        assert_eq!(alignment.path[0], (Some(0), Some(n1)));
        assert_eq!(alignment.path[1], (Some(1), Some(n2)));
        assert_eq!(alignment.path[2], (Some(2), Some(n3)));
        assert_eq!(alignment.path[3], (Some(3), Some(n4)));
        
        // Score: 4 * 5 = 20
        assert_eq!(alignment.score, 20);
    }

    #[test]
    fn test_align_mismatch() {
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        let n2 = graph.add_node(b'C'); // Mismatch here
        graph.add_edge(n1, n2, 1);

        let engine = ScalarAlignmentEngine::new(AlignmentParams::default(), AlignmentType::Global);
        let seq = b"AG"; // G vs C
        let alignment = engine.align(seq, &graph);

        assert_eq!(alignment.path.len(), 2);
        assert_eq!(alignment.path[0], (Some(0), Some(n1)));
        assert_eq!(alignment.path[1], (Some(1), Some(n2)));
        
        // Score: 5 (match) + -4 (mismatch) = 1
        assert_eq!(alignment.score, 1);
    }

    #[test]
    fn test_align_insertion_in_sequence() {
        // Graph: A -> C
        // Seq: A G C (G is inserted relative to graph)
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        let n2 = graph.add_node(b'C');
        graph.add_edge(n1, n2, 1);

        let engine = ScalarAlignmentEngine::new(AlignmentParams::default(), AlignmentType::Global);
        let seq = b"AGC";
        let alignment = engine.align(seq, &graph);

        // Expected: (A, A), (G, -), (C, C)
        assert_eq!(alignment.path.len(), 3);
        assert_eq!(alignment.path[0], (Some(0), Some(n1)));
        assert!(alignment.path[1].1.is_none()); // Insertion
        assert_eq!(alignment.path[2], (Some(2), Some(n2)));

        // Score: 5 (match) + (-8 gap open) + 5 (match) = 2
        // Note: Gap cost is gap_open + (k-1)*gap_extend. For k=1, cost is gap_open (-8).
        assert_eq!(alignment.score, 2);
    }

    #[test]
    fn test_align_deletion_in_sequence() {
        // Graph: A -> G -> C
        // Seq: A C (G is deleted relative to graph)
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        let n2 = graph.add_node(b'G');
        let n3 = graph.add_node(b'C');
        graph.add_edge(n1, n2, 1);
        graph.add_edge(n2, n3, 1);

        let engine = ScalarAlignmentEngine::new(AlignmentParams::default(), AlignmentType::Global);
        let seq = b"AC";
        let alignment = engine.align(seq, &graph);

        // Expected: (A, A), (-, G), (C, C)
        assert_eq!(alignment.path.len(), 3);
        assert_eq!(alignment.path[0], (Some(0), Some(n1)));
        assert!(alignment.path[1].0.is_none()); // Deletion
        assert_eq!(alignment.path[2], (Some(1), Some(n3)));
        
        // Score: 5 (match) - 8 (gap open) + 5 (match) = 2
        // Note: Gap cost is gap_open + (k-1)*gap_extend. For k=1, cost is gap_open (-8).
        assert_eq!(alignment.score, 2);
    }

    #[test]
    fn test_align_branching() {
        // Graph: A -> (G | T) -> C
        // Seq: A T C
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        let n2 = graph.add_node(b'G');
        let n3 = graph.add_node(b'T');
        let n4 = graph.add_node(b'C');
        graph.add_edge(n1, n2, 1);
        graph.add_edge(n1, n3, 1);
        graph.add_edge(n2, n4, 1);
        graph.add_edge(n3, n4, 1);

        let engine = ScalarAlignmentEngine::new(AlignmentParams::default(), AlignmentType::Global);
        
        // Align ATC -> should pick path A -> T -> C
        let seq = b"ATC";
        let alignment = engine.align(seq, &graph);

        assert_eq!(alignment.path.len(), 3);
        assert_eq!(alignment.path[0], (Some(0), Some(n1)));
        assert_eq!(alignment.path[1], (Some(1), Some(n3))); // Should pick T
        assert_eq!(alignment.path[2], (Some(2), Some(n4)));
        
        assert_eq!(alignment.score, 15);
    }

    #[test]
    fn test_align_semi_global() {
        // Graph: A -> C -> G -> T -> A -> C
        // Seq: G T (aligns to middle)
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        let n2 = graph.add_node(b'C');
        let n3 = graph.add_node(b'G');
        let n4 = graph.add_node(b'T');
        let n5 = graph.add_node(b'A');
        let n6 = graph.add_node(b'C');
        graph.add_edge(n1, n2, 1);
        graph.add_edge(n2, n3, 1);
        graph.add_edge(n3, n4, 1);
        graph.add_edge(n4, n5, 1);
        graph.add_edge(n5, n6, 1);

        let engine = ScalarAlignmentEngine::new(AlignmentParams::default(), AlignmentType::SemiGlobal);
        let seq = b"GT";
        let alignment = engine.align(seq, &graph);

        // Should align to G->T (n3 -> n4) with score 10
        // And path should be just (G, n3), (T, n4). 
        // Semi-Global implies no penalty for start/end of graph.
        
        assert_eq!(alignment.score, 10);
        assert_eq!(alignment.path.len(), 2);
        assert_eq!(alignment.path[0], (Some(0), Some(n3)));
        assert_eq!(alignment.path[1], (Some(1), Some(n4)));
    }

    #[test]
    fn test_align_local() {
        // Graph: A -> A -> A -> T -> T -> A -> A
        // Seq: C C T T G G
        // Should align T T to T T.
        let mut graph = PoaGraph::new();
        let n1 = graph.add_node(b'A');
        let n2 = graph.add_node(b'A');
        let n3 = graph.add_node(b'A');
        let n4 = graph.add_node(b'T');
        let n5 = graph.add_node(b'T');
        let n6 = graph.add_node(b'A');
        let n7 = graph.add_node(b'A');
        graph.add_edge(n1, n2, 1);
        graph.add_edge(n2, n3, 1);
        graph.add_edge(n3, n4, 1);
        graph.add_edge(n4, n5, 1);
        graph.add_edge(n5, n6, 1);
        graph.add_edge(n6, n7, 1);

        let engine = ScalarAlignmentEngine::new(AlignmentParams::default(), AlignmentType::Local);
        let seq = b"CCTTGG";
        let alignment = engine.align(seq, &graph);

        // Score: T-T (5), T-T (5) = 10.
        // C-A mismatch is -4. 
        // Local alignment should pick just TT.
        
        assert_eq!(alignment.score, 10);
        assert_eq!(alignment.path.len(), 2);
        assert_eq!(alignment.path[0], (Some(2), Some(n4))); // T
        assert_eq!(alignment.path[1], (Some(3), Some(n5))); // T
    }
}
