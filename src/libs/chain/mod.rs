//! Chaining functionality for connecting alignment blocks.
//!
//! This module implements the chaining algorithm used to connect local alignments (blocks)
//! into larger chains, similar to the UCSC `axtChain` tool.
//!
//! # Core Components
//!
//! * [`algo`] - Data structures for efficient predecessor search (KD-tree).
//! * [`gap_calc`] - Gap cost calculation (linear and affine penalties).
//! * [`sub_matrix`] - DNA substitution matrices (e.g., HoxD55).
//! * [`connect`] - Core chaining logic (dynamic programming, overlap trimming).
//! * [`record`] - Data structures for reading/writing Chain format.
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

pub mod algo;
pub mod connect;
pub mod gap_calc;
pub mod record;
pub mod sub_matrix;

pub use algo::{ChainItem, KdTree};
pub use connect::{calc_block_score, chain_blocks, ChainableBlock, ScoreContext};
pub use gap_calc::GapCalc;
pub use record::{read_chains, Block, Chain, ChainData, ChainHeader, ChainReader};
pub use sub_matrix::SubMatrix;
