//! Chain sorting helpers.

use super::record::Chain;

/// Sort chains in place by score descending. If `renumber`, reassign ids
/// starting from 1.
pub fn sort_chains(chains: &mut [Chain], renumber: bool) {
    chains.sort_by(|a, b| b.header.score.total_cmp(&a.header.score));

    if renumber {
        for (i, chain) in chains.iter_mut().enumerate() {
            chain.header.id = (i + 1) as u64;
        }
    }
}
