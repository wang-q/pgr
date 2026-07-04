//! Chain sorting helpers.

use std::cmp::Ordering;

use super::record::Chain;

/// Sort chains in place by score descending. If `renumber`, reassign ids
/// starting from 1.
pub fn sort_chains(chains: &mut [Chain], renumber: bool) {
    chains.sort_by(|a, b| {
        b.header
            .score
            .partial_cmp(&a.header.score)
            .unwrap_or(Ordering::Equal)
    });

    if renumber {
        for (i, chain) in chains.iter_mut().enumerate() {
            chain.header.id = (i + 1) as u64;
        }
    }
}
