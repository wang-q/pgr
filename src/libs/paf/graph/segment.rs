//! Segment-level representation and alignment-splitting helpers for
//! [`super::PafGraph`].
//!
//! [`Segment`] is a forward-orientation region on a sequence. Alignments are
//! split at indels >= `min_var_len` into `Segment` pairs (one on each side of
//! the alignment), linked bidirectionally via [`AlignmentLink`].

use crate::libs::paf::cigar::CigarOp;
use std::collections::HashMap;

/// A forward-orientation region on a sequence.
#[derive(Clone, Debug)]
pub(super) struct Segment {
    pub seq_id: u32,
    pub start: i32, // 0-based, inclusive
    pub end: i32,   // exclusive
}

/// A bidirectional alignment link between two segments.
pub(super) struct AlignmentLink {
    pub a: usize, // index into `segments`
    pub b: usize, // index into `segments`
    // True if the alignment between a and b is reverse-strand.
    // Currently unused — coarse GFA always emits '+' orientation; reserved
    // for future rGFA support to tag segment orientation per link.
    #[allow(dead_code)]
    pub reverse: bool,
}

/// Walk a CIGAR and split at indels >= `min_var_len`, emitting segment pairs.
#[allow(clippy::too_many_arguments)]
pub(super) fn split_alignment(
    tid: u32,
    t_start: i32,
    t_end: i32,
    qid: u32,
    q_start: i32,
    q_end: i32,
    q_size: i32,
    reverse: bool,
    cigar: &[CigarOp],
    min_var_len: i32,
    segments: &mut Vec<Segment>,
    links: &mut Vec<AlignmentLink>,
) {
    if cigar.is_empty() {
        // No CIGAR: treat whole record as one segment pair.
        let (qs_fwd, qe_fwd) = fwd_query_coords(q_start, q_end, q_size, reverse);
        let a = segments.len();
        segments.push(Segment {
            seq_id: tid,
            start: t_start,
            end: t_end,
        });
        let b = segments.len();
        segments.push(Segment {
            seq_id: qid,
            start: qs_fwd,
            end: qe_fwd,
        });
        links.push(AlignmentLink { a, b, reverse });
        return;
    }

    let mut ct = t_start;
    let mut cq = q_start; // in alignment-orientation coords (forward if '+', reverse if '-')
    let mut seg_start_t = ct;
    let mut seg_start_q = cq;
    let mut have_seg = false;

    let flush = |seg_start_t: i32,
                 ct: i32,
                 seg_start_q: i32,
                 cq: i32,
                 segments: &mut Vec<Segment>,
                 links: &mut Vec<AlignmentLink>,
                 have_seg: &mut bool| {
        if !*have_seg || ct <= seg_start_t {
            *have_seg = false;
            return;
        }
        let (qs_fwd, qe_fwd) = fwd_query_coords(seg_start_q, cq, q_size, reverse);
        if qs_fwd >= qe_fwd {
            *have_seg = false;
            return;
        }
        let a = segments.len();
        segments.push(Segment {
            seq_id: tid,
            start: seg_start_t,
            end: ct,
        });
        let b = segments.len();
        segments.push(Segment {
            seq_id: qid,
            start: qs_fwd,
            end: qe_fwd,
        });
        links.push(AlignmentLink { a, b, reverse });
        *have_seg = false;
    };

    for op in cigar {
        let td = op.target_delta() as i32;
        let qd = op.query_delta() as i32;
        let is_indel = (op.op() == 'I' || op.op() == 'D') && op.len() as i32 >= min_var_len;

        if is_indel {
            // Close current segment before the large indel.
            flush(
                seg_start_t,
                ct,
                seg_start_q,
                cq,
                segments,
                links,
                &mut have_seg,
            );
            ct += td;
            cq += qd;
            seg_start_t = ct;
            seg_start_q = cq;
        } else {
            // Extend current segment (small indels stay within segment as variations).
            if !have_seg {
                seg_start_t = ct;
                seg_start_q = cq;
                have_seg = true;
            }
            ct += td;
            cq += qd;
        }
    }
    // Flush trailing segment.
    flush(
        seg_start_t,
        ct,
        seg_start_q,
        cq,
        segments,
        links,
        &mut have_seg,
    );
}

/// Convert alignment-orientation query coords to forward-strand coords.
/// For '+' strand, coords are already forward. For '-' strand, flip.
fn fwd_query_coords(qs: i32, qe: i32, q_size: i32, reverse: bool) -> (i32, i32) {
    if reverse {
        crate::libs::alignment::coords::reverse_range_pair(qs, qe, q_size)
    } else {
        (qs, qe)
    }
}

/// Create a novel (unaligned) node for a gap region, return its node id.
pub(super) fn novel_node_for(
    node_seqs: &mut Vec<Vec<u8>>,
    node_origins: &mut Vec<(String, i32)>,
    sid: u32,
    start: i32,
    end: i32,
    seqs: Option<&HashMap<String, Vec<u8>>>,
    id_to_name: &[String],
) -> u32 {
    let name = id_to_name
        .get(sid as usize)
        .map(|s| s.as_str())
        .unwrap_or("?");
    let bytes = if let Some(seqs_map) = seqs {
        let seq_bytes = seqs_map.get(name).map(|v| v.as_slice()).unwrap_or(&[]);
        let s = start.max(0) as usize;
        let e = (end as usize).min(seq_bytes.len());
        if s < e {
            seq_bytes[s..e].to_vec()
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };
    let node_id = node_seqs.len() as u32;
    node_seqs.push(bytes);
    node_origins.push((name.to_string(), start));
    node_id
}
