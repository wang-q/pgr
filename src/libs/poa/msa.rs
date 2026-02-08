use super::graph::PoaGraph;
use petgraph::graph::NodeIndex;
use std::collections::HashSet;

/// Generates a Multiple Sequence Alignment (MSA) from the POA graph.
/// Returns a vector of strings, where each string represents a sequence in the MSA.
/// Gaps are represented by '-'.
pub fn generate_msa(
    graph: &PoaGraph,
    sequences: &[Vec<u8>],
    paths: &[Vec<NodeIndex>],
) -> Vec<String> {
    if sequences.is_empty() {
        return Vec::new();
    }

    // 1. Topological sort to linearize the graph
    let sorted_nodes = graph.topological_sort();

    // 2. Identify columns (cliques)
    // We iterate through sorted_nodes. If a node is not visited, it starts a new column.
    // The column includes the node and all its aligned_nodes.
    let mut columns: Vec<HashSet<NodeIndex>> = Vec::new();
    let mut visited: HashSet<NodeIndex> = HashSet::new();

    for &node_idx in &sorted_nodes {
        if visited.contains(&node_idx) {
            continue;
        }

        let mut column = HashSet::new();
        column.insert(node_idx);
        visited.insert(node_idx);

        let node_data = &graph.graph[node_idx];
        for &aligned in &node_data.aligned_nodes {
            if visited.insert(aligned) {
                column.insert(aligned);
            }
        }

        columns.push(column);
    }

    // 3. Generate MSA
    let num_seqs = sequences.len();
    let mut msa = vec![String::new(); num_seqs];
    let mut current_indices = vec![0usize; num_seqs];

    for column in &columns {
        for i in 0..num_seqs {
            let current_idx = current_indices[i];
            let path = &paths[i];

            if current_idx < path.len() {
                let node_in_path = path[current_idx];
                if column.contains(&node_in_path) {
                    // Match: Sequence visits a node in this column
                    msa[i].push(sequences[i][current_idx] as char);
                    current_indices[i] += 1;
                } else {
                    // Gap: Sequence does not visit this column (path node is in a future column)
                    msa[i].push('-');
                }
            } else {
                // End of sequence, fill with gaps
                msa[i].push('-');
            }
        }
    }

    msa
}

#[cfg(test)]
mod tests {
    use crate::libs::poa::{AlignmentParams, AlignmentType, Poa};

    #[test]
    fn test_msa_empty() {
        let params = AlignmentParams::default();
        let poa = Poa::new(params, AlignmentType::Global);
        let msa = poa.msa();
        assert!(msa.is_empty());
    }

    #[test]
    fn test_msa_identical() {
        let params = AlignmentParams::default();
        let mut poa = Poa::new(params, AlignmentType::Global);
        poa.add_sequence(b"ACGT");
        poa.add_sequence(b"ACGT");

        let msa = poa.msa();
        assert_eq!(msa.len(), 2);
        assert_eq!(msa[0], "ACGT");
        assert_eq!(msa[1], "ACGT");
    }

    #[test]
    fn test_msa_mismatch() {
        let params = AlignmentParams::default();
        let mut poa = Poa::new(params, AlignmentType::Global);
        // A C G T
        // A T G T
        poa.add_sequence(b"ACGT");
        poa.add_sequence(b"ATGT");

        let msa = poa.msa();
        assert_eq!(msa.len(), 2);
        // Depending on alignment parameters, it might align C and T or gap them.
        // Default match=2, mismatch=-4, gap_open=-4, gap_extend=-2.
        // Mismatch (-4) is costly.
        // ACGT
        // ATGT
        // M M M M -> 2-4+2+2 = 2
        // A-CGT
        // AT-GT
        // M G G M M -> 2-4-4+2+2 = -2
        // So it should align C and T as mismatch.

        // Wait, POA topological sort might order them differently if they are distinct nodes.
        // If C and T are aligned in the graph (mismatch), they should be in the same column.
        // If they are not aligned (branched), they might be in different columns.
        // Let's check the output.
        // With standard scoring, C and T should be aligned as a mismatch node or separate nodes aligned to each other?
        // In SPOA/POA, if they align, they are edges.
        // In this implementation, `aligned_nodes` tracks vertical alignment (topological alignment).
        // If the alignment engine aligns C to T, `add_alignment` should mark them as aligned.

        // Let's verify what happens.
        assert_eq!(msa[0], "ACGT");
        assert_eq!(msa[1], "ATGT");
    }

    #[test]
    fn test_msa_gap() {
        let params = AlignmentParams::default();
        let mut poa = Poa::new(params, AlignmentType::Global);
        poa.add_sequence(b"ACGT");
        poa.add_sequence(b"ACT"); // Deletion of G

        let msa = poa.msa();
        assert_eq!(msa.len(), 2);
        assert_eq!(msa[0], "ACGT");
        assert_eq!(msa[1], "AC-T");
    }

    #[test]
    fn test_msa_insertion() {
        let params = AlignmentParams::default();
        let mut poa = Poa::new(params, AlignmentType::Global);
        poa.add_sequence(b"ACT");
        poa.add_sequence(b"ACGT"); // Insertion of G

        let msa = poa.msa();
        assert_eq!(msa.len(), 2);
        assert_eq!(msa[0], "AC-T");
        assert_eq!(msa[1], "ACGT");
    }

    #[test]
    fn test_msa_complex() {
        let params = AlignmentParams::default();
        let mut poa = Poa::new(params, AlignmentType::Global);

        // 1. ACGT
        poa.add_sequence(b"ACGT");
        // 2. AC-T (Deletion of G)
        poa.add_sequence(b"ACT");
        // 3. A-GT (Deletion of C)
        poa.add_sequence(b"AGT");

        let msa = poa.msa();
        assert_eq!(msa.len(), 3);

        // Expected:
        // ACGT
        // AC-T
        // A-GT
        // Columns: A, C, G, T
        // 1: A C G T
        // 2: A C - T
        // 3: A - G T

        assert_eq!(msa[0], "ACGT");
        assert_eq!(msa[1], "AC-T");
        assert_eq!(msa[2], "A-GT");
    }
}
