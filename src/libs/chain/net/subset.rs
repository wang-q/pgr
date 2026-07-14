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

    // Process current fill only if type matches (when a filter is set).
    if type_filter.is_none_or(|t| &fill.class == t) && fill.chain_id != 0 {
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

    // Recurse into children gaps regardless of type filter.
    for gap in &fill.gaps {
        process_gap(gap, chains_map, writer, opts, type_filter)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::chain::net::types::{Chrom, Fill, Gap};
    use crate::libs::chain::record::{Chain, ChainData, ChainHeader};
    use std::cell::RefCell;
    use std::rc::Rc;

    fn make_chain(id: u64, t_start: u64, t_end: u64, q_start: u64, q_end: u64) -> Chain {
        Chain {
            header: ChainHeader {
                score: 100.0,
                t_name: "chr1".to_string(),
                t_size: 1000,
                t_strand: '+',
                t_start,
                t_end,
                q_name: "chr2".to_string(),
                q_size: 1000,
                q_strand: '+',
                q_start,
                q_end,
                id,
            },
            data: vec![ChainData {
                size: t_end - t_start,
                dt: 0,
                dq: 0,
            }],
        }
    }

    #[test]
    fn test_subset_type_filter_recurses_into_children() {
        // Top fill (type "top") contains a gap with a nested fill (type "syn").
        let chrom = Chrom::new("chr1", 1000);
        let nested_fill = Rc::new(RefCell::new(Fill {
            start: 10,
            end: 20,
            o_start: 10,
            o_end: 20,
            o_chrom: "chr2".to_string(),
            o_strand: '+',
            chain_id: 2,
            score: 50.0,
            ali: 10,
            class: "syn".to_string(),
            q_dup: None,
            q_over: None,
            q_far: None,
            chain: None,
            gaps: Vec::new(),
            t_n: None,
            q_n: None,
            t_r: None,
            q_r: None,
            t_trf: None,
            q_trf: None,
        }));
        let gap = Rc::new(RefCell::new(Gap {
            start: 10,
            end: 20,
            o_start: 10,
            o_end: 20,
            fills: vec![nested_fill],
            t_n: None,
            q_n: None,
            t_r: None,
            q_r: None,
            t_trf: None,
            q_trf: None,
        }));
        let top_fill = Rc::new(RefCell::new(Fill {
            start: 0,
            end: 100,
            o_start: 0,
            o_end: 100,
            o_chrom: "chr2".to_string(),
            o_strand: '+',
            chain_id: 1,
            score: 100.0,
            ali: 100,
            class: "top".to_string(),
            q_dup: None,
            q_over: None,
            q_far: None,
            chain: None,
            gaps: vec![gap],
            t_n: None,
            q_n: None,
            t_r: None,
            q_r: None,
            t_trf: None,
            q_trf: None,
        }));
        chrom.root.borrow_mut().fills.push(top_fill);

        let mut chains_map = HashMap::new();
        chains_map.insert(1, make_chain(1, 0, 100, 0, 100));
        chains_map.insert(2, make_chain(2, 10, 20, 10, 20));

        let mut buf = Vec::new();
        let opts = SubsetOptions {
            whole_chains: false,
            split_on_insert: false,
        };
        let type_filter = Some("syn".to_string());
        subset_nets(&[chrom], &chains_map, &mut buf, opts, type_filter.as_ref()).unwrap();
        let output = String::from_utf8(buf).unwrap();

        // The top fill should be skipped, but the nested "syn" fill must still be emitted.
        assert!(!output.contains("chain 100 chr1 1000 + 0 100 chr2 1000 + 0 100 1"));
        assert!(output.contains("chain 100 chr1 1000 + 10 20 chr2 1000 + 10 20 2"));
    }
}
