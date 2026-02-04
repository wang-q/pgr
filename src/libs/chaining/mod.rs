pub mod gap_calc;
pub mod sub_matrix;
pub mod seq;
pub mod record;
pub mod algo;

pub use gap_calc::GapCalc;
pub use sub_matrix::SubMatrix;
pub use record::{Chain, ChainHeader, ChainData, Block, ChainReader, read_chains};
pub use algo::{KdTree, ChainItem};
pub use seq::{chain_blocks, ChainableBlock, ScoreContext, calc_block_score};
