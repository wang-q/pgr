//! Chaining functionality for connecting alignment blocks.
//!
//! This module implements the chaining algorithm used to connect local alignments (blocks)
//! into larger chains, similar to the UCSC `axtChain` tool.
//!
//! # Core Components
//!
//! * [`kdtree`] - Data structures for efficient predecessor search (KD-tree).
//! * [`gap_calc`] - Gap cost calculation (linear and affine penalties).
//! * [`sub_matrix`] - DNA substitution matrices (e.g., HoxD55).
//! * [`connect`] - Core chaining logic (dynamic programming, overlap trimming).
//! * [`record`] - Data structures for reading/writing Chain format.
//!
//! Generic data structures such as [`BitMap`] live in [`crate::libs::ds`].
//!
//! # Algorithm Overview
//!
//! 1. **Input**: A set of alignment blocks (e.g., from PSL).
//! 2. **Indexing**: Blocks are indexed in a KD-tree by their start coordinates (query, target).
//! 3. **Dynamic Programming**:
//!    - For each block, find "predecessor" blocks that are strictly before it in both query and target.
//!    - Calculate score: `Score = BlockScore + Max(PredecessorScore - GapCost)`.
//!    - Gap costs are calculated based on the distance between blocks.
//! 4. **Overlap Handling**: Overlaps between adjacent blocks are trimmed to maximize the total score.
//! 5. **Output**: Chains passing a minimum score threshold are output.

pub mod anti_repeat;
pub mod connect;
pub mod gap_calc;
pub mod kdtree;
pub mod net;
pub mod pre_net;
pub mod psl_chain;
pub mod record;
pub mod sort;
pub mod stitch;
pub mod sub_matrix;

pub use connect::{calc_block_score, chain_blocks, ChainableBlock, ScoreContext};
pub use gap_calc::GapCalc;
pub use kdtree::{ChainItem, KdTree};
pub use pre_net::{is_haplotype, pre_net, PreNetOptions};
pub use psl_chain::{chain_psl, group_psl_blocks, GroupData, GroupKey};
pub use record::{read_chains, Block, Chain, ChainData, ChainHeader, ChainReader};
pub use sort::sort_chains;
pub use stitch::stitch_chains;
pub use sub_matrix::SubMatrix;

/// Derive a 3-digit lump bucket name from a sequence name.
///
/// Scans `name` for the first ASCII-digit run and returns `val % lump` formatted
/// as 3 digits. If no digits are found, falls back to a stable hash of the name.
pub fn lump_name(name: &str, lump: usize) -> String {
    // Look for integer part of name
    let mut s = name;
    while let Some(idx) = s.find(|c: char| c.is_ascii_digit()) {
        s = &s[idx..];
        let end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
        let digits = &s[..end];
        if let Ok(val) = digits.parse::<usize>() {
            return format!("{:03}", val % lump);
        }
        s = &s[end..];
    }

    // If no digits found, hash it with fxhash for stability across Rust versions.
    let hash = fxhash::hash64(name.as_bytes());
    format!("{:03}", (hash as usize) % lump)
}
