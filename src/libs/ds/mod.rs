//! General-purpose data structures used across the library.

pub mod bitmap;
pub mod crossover;
pub mod dupe_tree;
pub mod gap_calc;
pub mod interval;
pub mod kdtree;
pub mod top_k_purity;

pub use bitmap::BitMap;
pub use crossover::best_crossover;
pub use dupe_tree::{DupeTree, Segment};
pub use gap_calc::GapCalc;
pub use interval::merge_intervals;
pub use kdtree::{KdTree, KdTreeItem};
pub use top_k_purity::TopKPurity;
