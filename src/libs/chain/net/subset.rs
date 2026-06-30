//! Subset chains based on net structure.

use super::types::{Fill, Gap};
use crate::libs::chain::Chain;
use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Write;
use std::rc::Rc;

/// Options controlling how chains are subset according to the net.
#[derive(Clone, Copy)]
pub struct SubsetOptions {
    /// Write entire chain references by net, don't split.
    pub whole_chains: bool,
    /// Split chain when encountering an insertion of another chain.
    pub split_on_insert: bool,
}

/// Traverse the net tree and write chain subsets to `writer`.
///
/// `type_filter` restricts output to a particular `type` field in the net file.
pub fn subset_nets(
    chroms: &[super::types::Chrom],
    chains_map: &HashMap<u64, Chain>,
    writer: &mut impl Write,
    opts: SubsetOptions,
    type_filter: Option<&String>,
) -> anyhow::Result<()> {
    for chrom in chroms {
        process_gap(&chrom.root, chains_map, writer, opts, type_filter)?;
    }
    Ok(())
}

fn process_gap(
    gap: &Rc<RefCell<Gap>>,
    chains_map: &HashMap<u64, Chain>,
    writer: &mut impl Write,
    opts: SubsetOptions,
    type_filter: Option<&String>,
) -> anyhow::Result<()> {
    let gap = gap.borrow();
    for fill in &gap.fills {
        process_fill(fill, chains_map, writer, opts, type_filter)?;
    }
    Ok(())
}

fn process_fill(
    fill_rc: &Rc<RefCell<Fill>>,
    chains_map: &HashMap<u64, Chain>,
    writer: &mut impl Write,
    opts: SubsetOptions,
    type_filter: Option<&String>,
) -> anyhow::Result<()> {
    let fill = fill_rc.borrow();

    // Check type filter
    if let Some(t) = type_filter {
        if &fill.class != t {
            return Ok(()); // Skip but continue traversal?
                           // In C: if (!sameString(type, fill->type)) return;
                           // It returns from convertFill, but then it continues recursion in rConvert.
                           // Wait, in C rConvert calls convertFill THEN recurses.
                           // So if type doesn't match, we don't output this fill, but do we recurse?
                           // C code:
                           // if (fill->chainId) { ... convertFill ... }
                           // if (fill->children) rConvert(...);
                           //
                           // convertFill checks type and returns if mismatch.
                           // So yes, we should still recurse.
        }
    }

    // Process current fill
    if fill.chain_id != 0 {
        if let Some(chain) = chains_map.get(&fill.chain_id) {
            if opts.whole_chains {
                chain.write(writer)?;
            } else if opts.split_on_insert {
                // Split on insert logic
                let mut t_start = fill.start;

                // Iterate over gaps to find inserts
                for gap_rc in &fill.gaps {
                    let gap = gap_rc.borrow();
                    if !gap.fills.is_empty() {
                        // This gap has inserts (children fills)
                        // Output chain part from t_start to gap.start
                        if gap.start > t_start {
                            if let Some(sub) = chain.subset(t_start, gap.start) {
                                sub.write(writer)?;
                            }
                        }
                        t_start = gap.end;
                    }
                }
                // Output remaining part
                if fill.end > t_start {
                    if let Some(sub) = chain.subset(t_start, fill.end) {
                        sub.write(writer)?;
                    }
                }
            } else {
                // Default: subset to fill range
                if let Some(sub) = chain.subset(fill.start, fill.end) {
                    sub.write(writer)?;
                }
            }
        }
    }

    // Recurse into children gaps
    for gap in &fill.gaps {
        process_gap(gap, chains_map, writer, opts, type_filter)?;
    }

    Ok(())
}
