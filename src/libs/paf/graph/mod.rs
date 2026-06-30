//! Coarse GFA graph induction from PAF (seqwish-style segment-level DSU).
//!
//! Splits each alignment at indels >= `min_var_len` into "match segments",
//! unions aligned segments via a disjoint-set union (transitive closure),
//! derives graph nodes (DSU classes) + edges (path adjacencies) + novel
//! segments (unaligned gaps), and emits GFA v1.0 (S/L/P).

mod builder;
mod dsu;
mod gfa;
mod report;
mod segment;
#[cfg(test)]
mod tests;

pub use report::GraphReport;

/// A graph edge (GFA L line).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Edge {
    pub from: u32,
    pub from_orient: char,
    pub to: u32,
    pub to_orient: char,
}

/// A path entry (one node visit in a GFA P line).
#[derive(Clone, Debug)]
pub struct PathStep {
    pub node: u32,
    pub orient: char,
}

/// Induced coarse graph, ready for GFA emission.
pub struct PafGraph {
    /// Node sequences (one per DSU class), indexed by node id.
    /// Empty when built without FASTA (topology-only mode).
    pub node_seqs: Vec<Vec<u8>>,
    /// Node lengths (bp). Populated even in topology-only mode from segment coords.
    pub node_lens: Vec<usize>,
    /// Per-node rGFA origin: (source sequence name, 0-based start offset).
    pub node_origins: Vec<(String, i32)>,
    /// Edges (deduplicated).
    pub edges: Vec<Edge>,
    /// Paths: (sequence name, steps).
    pub paths: Vec<(String, Vec<PathStep>)>,
}
