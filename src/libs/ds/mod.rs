//! General-purpose data structures used across the library.

pub mod bitmap;
pub mod dupe_tree;
pub mod gap_calc;
pub mod kdtree;

pub use bitmap::BitMap;
pub use dupe_tree::{DupeTree, Segment};
pub use gap_calc::GapCalc;
pub use kdtree::{KdTree, KdTreeItem};
