//! chainPreNet: filter chains that cannot be netted, marking used target/query ranges.

use super::record::ChainReader;
use crate::libs::ds::BitMap;
use anyhow::{bail, Result};
use std::collections::HashMap;
use std::io::{BufRead, Write};

/// Return true if the sequence name looks like a haplotype/alt contig.
///
/// Matches the UCSC `haplotype()` heuristic: contains "_hap" or "_alt".
pub fn is_haplotype(name: &str) -> bool {
    name.contains("_hap") || name.contains("_alt")
}

/// Options for `pre_net`.
pub struct PreNetOptions {
    pub pad: u64,
    pub incl_hap: bool,
    pub dots: Option<usize>,
}

/// Run chainPreNet: filter chains, mark used ranges in target/query bitmaps.
///
/// Reads chains from `reader`, writes passing chains to `writer`. `t_hash` and
/// `q_hash` are mutated in place to track used ranges.
pub fn pre_net<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
    t_hash: &mut HashMap<String, BitMap>,
    q_hash: &mut HashMap<String, BitMap>,
    opts: &PreNetOptions,
) -> Result<()> {
    let chain_reader = ChainReader::new(reader);
    let mut last_score = f64::MAX;

    for (count, res) in chain_reader.enumerate() {
        let chain = res?;

        // Check sort order
        if chain.header.score > last_score {
            bail!(
                "Input not sorted by score: {} > {}",
                chain.header.score,
                last_score
            );
        }
        last_score = chain.header.score;

        if let Some(d) = opts.dots {
            if count > 0 && count % d == 0 {
                eprint!(".");
            }
        }

        if !opts.incl_hap && is_haplotype(&chain.header.q_name) {
            continue;
        }

        let t_chrom = t_hash.get_mut(&chain.header.t_name).ok_or_else(|| {
            anyhow::anyhow!("Target sequence {} not found in sizes", chain.header.t_name)
        })?;
        let q_chrom = q_hash.get_mut(&chain.header.q_name).ok_or_else(|| {
            anyhow::anyhow!("Query sequence {} not found in sizes", chain.header.q_name)
        })?;

        let blocks = chain.to_blocks();
        let mut any_open = false;
        for b in &blocks {
            if !q_chrom.is_fully_set(b.q_start, b.q_end - b.q_start) {
                any_open = true;
                break;
            }
            if !t_chrom.is_fully_set(b.t_start, b.t_end - b.t_start) {
                any_open = true;
                break;
            }
        }

        if any_open {
            chain.write(&mut writer)?;
            for b in &blocks {
                let q_s = b.q_start.saturating_sub(opts.pad);
                let q_len = (b.q_end + opts.pad).min(q_chrom.size) - q_s;
                q_chrom.set_range(q_s, q_len);

                let t_s = b.t_start.saturating_sub(opts.pad);
                let t_len = (b.t_end + opts.pad).min(t_chrom.size) - t_s;
                t_chrom.set_range(t_s, t_len);
            }
        }
    }

    if opts.dots.is_some() {
        eprintln!();
    }
    Ok(())
}
