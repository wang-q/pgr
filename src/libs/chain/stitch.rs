//! Stitch chain fragments sharing the same chain ID into a single chain.

use super::record::{Chain, ChainReader};
use anyhow::Result;
use std::collections::HashMap;
use std::io::{BufRead, Write};

/// Read chains from `reader`, merge fragments with the same ID, and write stitched chains to `writer`.
///
/// Fragments are merged by converting to blocks, sorting by t_start, and rebuilding.
/// Output is sorted by score descending.
pub fn stitch_chains<R: BufRead, W: Write>(reader: R, mut writer: W) -> Result<()> {
    let chain_reader = ChainReader::new(reader);
    let mut chains: HashMap<u64, Chain> = HashMap::new();

    for res in chain_reader {
        let chain = res?;
        chains
            .entry(chain.header.id)
            .and_modify(|existing| {
                // Merge logic: tackOnFrag
                // Check consistency
                if existing.header.t_name != chain.header.t_name
                    || existing.header.q_name != chain.header.q_name
                    || existing.header.q_strand != chain.header.q_strand
                {
                    log::warn!(
                        "Inconsistent chain info for ID {}: skipping fragment",
                        chain.header.id
                    );
                    return;
                }

                // Convert both to blocks
                let mut blocks = existing.to_blocks();
                let new_blocks = chain.to_blocks();

                // Append new blocks
                blocks.extend(new_blocks);

                // Sort blocks by t_start, then q_start for deterministic ordering
                blocks.sort_by_key(|a| (a.t_start, a.q_start));

                // Reconstruct data from blocks
                // This updates header ranges automatically
                existing.data = Chain::from_blocks(&mut existing.header, &blocks);

                // Sum score
                existing.header.score += chain.header.score;
            })
            .or_insert(chain);
    }

    // Collect and sort by score (descending)
    let mut chain_list: Vec<Chain> = chains.into_values().collect();
    chain_list.sort_by(|a, b| b.header.score.total_cmp(&a.header.score));

    for chain in chain_list {
        chain.write(&mut writer)?;
    }

    Ok(())
}
