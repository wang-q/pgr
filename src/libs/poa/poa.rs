use super::align::{AlignmentEngine, AlignmentParams, AlignmentType, ScalarAlignmentEngine};
use super::consensus::generate_consensus;
use super::graph::PoaGraph;
use super::msa::generate_msa;
use petgraph::graph::NodeIndex;

pub struct Poa {
    graph: PoaGraph,
    engine: Box<dyn AlignmentEngine>,
    sequences: Vec<Vec<u8>>,
    paths: Vec<Vec<NodeIndex>>,
}

impl Poa {
    pub fn new(params: AlignmentParams, mode: AlignmentType) -> Self {
        Self {
            graph: PoaGraph::new(),
            engine: Box::new(ScalarAlignmentEngine::new(params, mode)),
            sequences: Vec::new(),
            paths: Vec::new(),
        }
    }

    pub fn add_sequence(&mut self, sequence: &[u8]) {
        let alignment = self.engine.align(sequence, &self.graph);
        let path = self.graph.add_alignment(&alignment, sequence);
        self.sequences.push(sequence.to_vec());
        self.paths.push(path);
    }

    pub fn consensus(&self) -> Vec<u8> {
        generate_consensus(&self.graph)
    }

    pub fn msa(&self) -> Vec<String> {
        generate_msa(&self.graph, &self.sequences, &self.paths)
    }

    pub fn num_nodes(&self) -> usize {
        self.graph.num_nodes()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poa_workflow() {
        let params = AlignmentParams::default();
        let mut poa = Poa::new(params, AlignmentType::Global);

        // 1. Add first sequence (creates linear graph)
        poa.add_sequence(b"ACGT");
        assert_eq!(poa.num_nodes(), 4);
        assert_eq!(poa.consensus(), b"ACGT");

        // 2. Add second sequence (identical)
        poa.add_sequence(b"ACGT");
        assert_eq!(poa.num_nodes(), 4);
        assert_eq!(poa.consensus(), b"ACGT");

        // 3. Add third sequence (mutation: A -> T at pos 0)
        // ACGT
        // ACGT
        // TCGT
        // Consensus should likely be ACGT (A:2, T:1)
        poa.add_sequence(b"TCGT");
        assert_eq!(poa.num_nodes(), 5); // New node for T
        assert_eq!(poa.consensus(), b"ACGT");
    }

    #[test]
    fn test_poa_insertion() {
        let params = AlignmentParams::default();
        let mut poa = Poa::new(params, AlignmentType::Global);

        poa.add_sequence(b"ACGT");
        // Add sequence with insertion: ACAGT
        poa.add_sequence(b"ACAGT");

        // Graph should have branched or added node
        // A -> C -> G -> T
        //      |
        //      A -> G
        // Ideally consensus should handle this.
        // A C A G T vs A C G T
        // If we just added 2 sequences, consensus might pick one path.
        // Let's add more weight to one.
        poa.add_sequence(b"ACAGT");

        assert_eq!(poa.consensus(), b"ACAGT");
    }

    #[test]
    fn test_poa_msa() {
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
        // MSA should be:
        // ACGT
        // AC-T
        // A-GT
        assert_eq!(msa[0], "ACGT");
        assert_eq!(msa[1], "AC-T");
        assert_eq!(msa[2], "A-GT");
    }
}
