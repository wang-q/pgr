//! ChainNet builder: inserts chains into chromosome gap trees.
//!
//! Each chain's alignment blocks are placed into the target (or query)
//! chromosome's gap tree, producing a hierarchical Net structure.

use super::types::{Chrom, Fill, Gap, Space};
use crate::libs::alignment::coords::reverse_range;
use crate::libs::chain::record::{Block, Chain};
use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap};
use std::rc::Rc;

/// A collection of per-chromosome net trees built from chains.
pub struct ChainNet {
    pub chroms: HashMap<String, RefCell<Chrom>>,
    pub chains: Vec<Rc<Chain>>,
}

impl ChainNet {
    pub fn new(target_sizes: &BTreeMap<String, u64>) -> Self {
        let mut chroms = HashMap::new();
        for (name, size) in target_sizes {
            chroms.insert(name.clone(), RefCell::new(Chrom::new(name, *size)));
        }
        Self {
            chroms,
            chains: Vec::new(),
        }
    }

    pub fn add_chain(&mut self, chain: Chain, min_space: u64, min_fill: u64, min_score: f64) {
        if chain.header.score < min_score {
            return;
        }
        let chain_rc = Rc::new(chain);
        self.chains.push(chain_rc.clone());

        // Add to target net
        if let Some(chrom) = self.chroms.get(&chain_rc.header.t_name) {
            let mut chrom = chrom.borrow_mut();
            let blocks = chain_rc.to_blocks();
            add_chain_core(
                &mut chrom,
                chain_rc.clone(),
                blocks,
                false,
                min_space,
                min_fill,
            );
        }
    }

    pub fn add_chain_as_q(&mut self, chain: Chain, min_space: u64, min_fill: u64, min_score: f64) {
        if chain.header.score < min_score {
            return;
        }
        let chain_rc = Rc::new(chain);
        self.chains.push(chain_rc.clone());

        if let Some(chrom) = self.chroms.get(&chain_rc.header.q_name) {
            let mut chrom = chrom.borrow_mut();
            let mut blocks = chain_rc.to_blocks();

            if chain_rc.header.q_strand == '-' {
                reverse_blocks_q(&mut blocks, chain_rc.header.q_size);
            }

            add_chain_core(
                &mut chrom,
                chain_rc.clone(),
                blocks,
                true,
                min_space,
                min_fill,
            );
        }
    }
}

fn reverse_blocks_q(blocks: &mut [Block], size: u64) {
    blocks.reverse();
    for b in blocks {
        reverse_range(&mut b.q_start, &mut b.q_end, size);
    }
}

fn add_chain_core(
    chrom: &mut Chrom,
    chain: Rc<Chain>,
    blocks: Vec<Block>,
    is_q: bool,
    min_space: u64,
    min_fill: u64,
) {
    let (start, end) = if is_q {
        let mut s = chain.header.q_start;
        let mut e = chain.header.q_end;
        if chain.header.q_strand == '-' {
            reverse_range(&mut s, &mut e, chain.header.q_size);
        }
        (s, e)
    } else {
        (chain.header.t_start, chain.header.t_end)
    };

    let spaces = chrom.find_spaces(start, end);
    let mut start_block_idx = 0;

    for space in spaces {
        let mut first_idx = None;
        let mut last_idx = None;
        let mut s = u64::MAX;
        let mut e = 0;

        for (i, b) in blocks.iter().enumerate().skip(start_block_idx) {
            let (b_start, b_end) = if is_q {
                (b.q_start, b.q_end)
            } else {
                (b.t_start, b.t_end)
            };

            if b_end <= space.start {
                continue;
            }
            if b_start >= space.end {
                break;
            }

            if first_idx.is_none() {
                first_idx = Some(i);
            }
            last_idx = Some(i);

            let curr_s: u64 = b_start.max(space.start);
            let curr_e: u64 = b_end.min(space.end);

            if curr_s < s {
                s = curr_s;
            }
            if curr_e > e {
                e = curr_e;
            }
        }

        if let Some(idx) = first_idx {
            start_block_idx = idx;
        } else {
            continue;
        }

        if s >= e || (e - s) < min_fill {
            continue;
        }

        fill_space(
            chrom,
            space,
            chain.clone(),
            &blocks,
            first_idx.unwrap(),
            last_idx.unwrap(),
            s,
            e,
            min_space,
            is_q,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_space(
    chrom: &mut Chrom,
    space: Space,
    chain: Rc<Chain>,
    blocks: &[Block],
    first_idx: usize,
    last_idx: usize,
    fill_start: u64,
    fill_end: u64,
    min_space: u64,
    is_q: bool,
) {
    // Remove old space
    chrom.spaces.remove(&space.start);

    // Calculate other side coords for the fill
    let (o_start, o_end) = if !is_q {
        let b1 = &blocks[first_idx];
        let offset1 = fill_start - b1.t_start;
        let mut q1 = b1.q_start + offset1;

        let b2 = &blocks[last_idx];
        let offset2 = fill_end - b2.t_start;
        let mut q2 = b2.q_start + offset2;

        if chain.header.q_strand == '-' {
            reverse_range(&mut q1, &mut q2, chain.header.q_size);
        }

        (q1, q2)
    } else {
        let b1 = &blocks[first_idx];
        let offset1 = fill_start - b1.q_start;

        let t1 = if chain.header.q_strand == '-' {
            b1.t_end - offset1
        } else {
            b1.t_start + offset1
        };

        let b2 = &blocks[last_idx];
        let offset2 = fill_end - b2.q_start;
        let t2 = if chain.header.q_strand == '-' {
            b2.t_end - offset2
        } else {
            b2.t_start + offset2
        };

        if t1 > t2 {
            (t2, t1)
        } else {
            (t1, t2)
        }
    };

    let o_chrom = if is_q {
        &chain.header.t_name
    } else {
        &chain.header.q_name
    };
    let o_strand = chain.header.q_strand;

    // Create Fill
    let fill = Rc::new(RefCell::new(Fill {
        start: fill_start,
        end: fill_end,
        o_start,
        o_end,
        o_chrom: o_chrom.clone(),
        o_strand,
        chain_id: chain.header.id,
        score: 0.0,
        ali: 0,
        class: String::new(),
        q_dup: None,
        q_over: None,
        q_far: None,
        chain: Some(chain.clone()),
        gaps: Vec::new(),
        t_n: None,
        q_n: None,
        t_r: None,
        q_r: None,
        t_trf: None,
        q_trf: None,
    }));

    // Add Left Space
    if fill_start > space.start && (fill_start - space.start) >= min_space {
        chrom.spaces.insert(
            space.start,
            Space {
                start: space.start,
                end: fill_start,
                gap: space.gap.clone(),
            },
        );
    }

    // Add Right Space
    if fill_end < space.end && (space.end - fill_end) >= min_space {
        chrom.spaces.insert(
            fill_end,
            Space {
                start: fill_end,
                end: space.end,
                gap: space.gap.clone(),
            },
        );
    }

    // Internal gaps
    for i in first_idx..last_idx {
        let b1 = &blocks[i];
        let b2 = &blocks[i + 1];

        let (gap_start, gap_end) = if is_q {
            (b1.q_end, b2.q_start)
        } else {
            (b1.t_end, b2.t_start)
        };

        if gap_start > fill.borrow().start
            && gap_end < fill.borrow().end
            && (gap_end - gap_start) >= min_space
        {
            let (mut os, mut oe) = if !is_q {
                (b1.q_end, b2.q_start)
            } else if chain.header.q_strand == '-' {
                (b2.t_start, b1.t_end)
            } else {
                (b1.t_end, b2.t_start)
            };

            if !is_q && chain.header.q_strand == '-' {
                reverse_range(&mut os, &mut oe, chain.header.q_size);
            }

            let new_gap = Rc::new(RefCell::new(Gap {
                start: gap_start,
                end: gap_end,
                o_start: os,
                o_end: oe,
                fills: Vec::new(),
                t_n: None,
                q_n: None,
                t_r: None,
                q_r: None,
                t_trf: None,
                q_trf: None,
            }));

            chrom.spaces.insert(
                gap_start,
                Space {
                    start: gap_start,
                    end: gap_end,
                    gap: new_gap.clone(),
                },
            );

            fill.borrow_mut().gaps.push(new_gap);
        }
    }

    // Add fill to parent gap
    space.gap.borrow_mut().fills.push(fill);
}
