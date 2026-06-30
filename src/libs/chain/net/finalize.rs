//! Net finalization: sort fills/gaps by start and recompute o_start/o_end.
//!
//! After all chains have been inserted, `finalize_net` sorts each level of the
//! tree by `start` and recomputes the other-side coordinates (`o_start`/
//! `o_end`) for each fill by walking its attached chain's alignment data.

use super::types::Gap;
use crate::libs::alignment::coords::reverse_range;
use std::cell::RefCell;
use std::rc::Rc;

/// Sort the net tree and recompute fill o_start/o_end from chain data.
pub fn finalize_net(chrom: &mut super::types::Chrom, is_q: bool) {
    sort_net(&chrom.root);
    calc_other_fill(&chrom.root, is_q);
}

fn sort_net(gap: &Rc<RefCell<Gap>>) {
    let mut gap_borrow = gap.borrow_mut();
    gap_borrow.fills.sort_by_key(|f| f.borrow().start);

    for fill in &gap_borrow.fills {
        let mut fill_borrow = fill.borrow_mut();
        fill_borrow.gaps.sort_by_key(|g| g.borrow().start);
        for g in &fill_borrow.gaps {
            sort_net(g);
        }
    }
}

fn calc_other_fill(gap: &Rc<RefCell<Gap>>, is_q: bool) {
    let gap_borrow = gap.borrow();
    for fill in &gap_borrow.fills {
        let mut fill_borrow = fill.borrow_mut();

        if let Some(chain) = fill_borrow.chain.clone() {
            let clip_start = fill_borrow.start;
            let clip_end = fill_borrow.end;

            if !is_q {
                let mut q_min = u64::MAX;
                let mut q_max = 0;

                let mut t_curr = chain.header.t_start;
                let mut q_curr = chain.header.q_start;

                for d in &chain.data {
                    let t_s = t_curr;
                    let t_e = t_curr + d.size;
                    let q_s = q_curr;

                    let start = t_s.max(clip_start);
                    let end = t_e.min(clip_end);

                    if start < end {
                        let offset = start - t_s;
                        let len = end - start;
                        let qs = q_s + offset;
                        let qe = qs + len;

                        if qs < q_min {
                            q_min = qs;
                        }
                        if qe > q_max {
                            q_max = qe;
                        }
                    }

                    t_curr += d.size + d.dt;
                    q_curr += d.size + d.dq;
                }

                if q_min < q_max {
                    if chain.header.q_strand == '-' {
                        reverse_range(&mut q_min, &mut q_max, chain.header.q_size);
                    }
                    fill_borrow.o_start = q_min;
                    fill_borrow.o_end = q_max;
                }
            } else {
                let mut t_min = u64::MAX;
                let mut t_max = 0;

                let mut t_curr = chain.header.t_start;
                let mut q_curr = chain.header.q_start;

                for d in &chain.data {
                    let t_s = t_curr;
                    let q_s = q_curr;
                    let q_e = q_curr + d.size;

                    let (c_start, c_end) = (clip_start, clip_end);

                    let (mut fq_s, mut fq_e) = (q_s, q_e);
                    if chain.header.q_strand == '-' {
                        reverse_range(&mut fq_s, &mut fq_e, chain.header.q_size);
                    }

                    let start = fq_s.max(c_start);
                    let end = fq_e.min(c_end);

                    if start < end {
                        let len = end - start;
                        let (ts, te) = if chain.header.q_strand == '-' {
                            let rq_s = chain.header.q_size - end;
                            let offset = rq_s - q_s;
                            let ts = t_s + offset;
                            (ts, ts + len)
                        } else {
                            let offset = start - q_s;
                            let ts = t_s + offset;
                            (ts, ts + len)
                        };

                        if ts < t_min {
                            t_min = ts;
                        }
                        if te > t_max {
                            t_max = te;
                        }
                    }

                    t_curr += d.size + d.dt;
                    q_curr += d.size + d.dq;
                }

                if t_min < t_max {
                    fill_borrow.o_start = t_min;
                    fill_borrow.o_end = t_max;
                }
            }
        }

        drop(fill_borrow);
        for g in &fill.borrow().gaps {
            calc_other_fill(g, is_q);
        }
    }
}
