//! Net-to-AXT conversion: emit AXT alignments from a net tree + chains + 2bit.
//!
//! Walks each fill in the net tree, looks up its source chain, and emits AXT
//! records by slicing target/query sequences from 2bit files and interleaving
//! aligned blocks with gap-derived dashes.

use super::types::{Chrom, Fill, Gap};
use crate::libs::chain::record::Chain;
use crate::libs::chain::sub_matrix::SubMatrix;
use crate::libs::io::SequenceReader;
use crate::libs::nt;
use anyhow::anyhow;
use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Write;
use std::rc::Rc;

/// Convert all nets to AXT records written to `writer`. Returns the number of
/// AXT records emitted.
pub fn net_to_axt<S: SequenceReader, W: Write>(
    nets: &[Chrom],
    chains: &HashMap<u64, Chain>,
    t_2bit: &mut S,
    q_2bit: &mut S,
    matrix: &SubMatrix,
    writer: &mut W,
) -> anyhow::Result<usize> {
    // Write header comments from the first net (if any).
    if let Some(first_net) = nets.first() {
        for comment in &first_net.comments {
            writeln!(writer, "{}", comment)?;
        }
    }

    let mut counter = 0;
    for net in nets {
        r_convert(
            &net.root,
            chains,
            t_2bit,
            q_2bit,
            matrix,
            writer,
            &mut counter,
        )?;
    }
    Ok(counter)
}

fn r_convert<S: SequenceReader, W: Write>(
    gap: &Rc<RefCell<Gap>>,
    chains: &HashMap<u64, Chain>,
    t_2bit: &mut S,
    q_2bit: &mut S,
    matrix: &SubMatrix,
    writer: &mut W,
    counter: &mut usize,
) -> anyhow::Result<()> {
    let g = gap.borrow();
    for fill in &g.fills {
        let f = fill.borrow();
        if f.chain_id != 0 {
            if let Some(chain) = chains.get(&f.chain_id) {
                convert_fill(&f, chain, chains, t_2bit, q_2bit, matrix, writer, counter)?;
            }
        } else {
            // If no chain, just recurse into gaps
            for gap_rc in &f.gaps {
                r_convert(gap_rc, chains, t_2bit, q_2bit, matrix, writer, counter)?;
            }
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn convert_fill<S: SequenceReader, W: Write>(
    fill: &Fill,
    chain: &Chain,
    chains: &HashMap<u64, Chain>,
    t_2bit: &mut S,
    q_2bit: &mut S,
    matrix: &SubMatrix,
    writer: &mut W,
    counter: &mut usize,
) -> anyhow::Result<()> {
    // UCSC netToAxt splitWrite + rConvert: first emit ALL segments of this
    // fill (splitting at gaps with inserts), THEN recurse into children.
    // This produces parent AXT records before child AXT records (pre-order).
    let mut cur = fill.start;
    let mut child_gaps: Vec<Rc<RefCell<Gap>>> = Vec::new();

    // First pass: emit all segments (UCSC splitWrite behavior).
    for gap_rc in &fill.gaps {
        let (g_start, g_end, has_children) = {
            let g = gap_rc.borrow();
            (g.start, g.end, !g.fills.is_empty())
        };

        // Mirror UCSC splitWrite: split only at gaps that carry an inserted
        // child fill; plain gaps are covered by the surrounding segment and
        // rendered as dashes (chainToAxt then splits single-sided gaps that
        // exceed maxGap).
        let should_split = has_children;

        if should_split {
            if g_start > cur {
                convert_segment(cur, g_start, chain, t_2bit, q_2bit, matrix, writer, counter)?;
            }
            cur = cur.max(g_end);
            if has_children {
                child_gaps.push(gap_rc.clone());
            }
        }
    }

    // Tail segment.
    if cur < fill.end {
        convert_segment(
            cur, fill.end, chain, t_2bit, q_2bit, matrix, writer, counter,
        )?;
    }

    // Second pass: recurse into children (UCSC rConvert behavior).
    for gap_rc in child_gaps {
        r_convert(&gap_rc, chains, t_2bit, q_2bit, matrix, writer, counter)?;
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn convert_segment<S: SequenceReader, W: Write>(
    t_start: u64,
    t_end: u64,
    chain: &Chain,
    t_2bit: &mut S,
    q_2bit: &mut S,
    matrix: &SubMatrix,
    writer: &mut W,
    counter: &mut usize,
) -> anyhow::Result<()> {
    let blocks = chain.to_blocks();

    // Mirror UCSC chainToAxt: break the segment wherever a block gap is
    // double-sided (dq > 0 && dt > 0) or a single-sided gap exceeds maxGap.
    // Splits happen at block boundaries; the first split point may be before
    // t_start when the segment begins mid-gap, so clamp to t_start.
    const MAX_GAP: i64 = 100;
    let mut segments: Vec<(u64, u64)> = Vec::new();
    let mut seg_start = t_start;

    for pair in blocks.windows(2) {
        let a = &pair[0];
        let b = &pair[1];
        if b.t_start <= t_start {
            continue;
        }
        if a.t_end >= t_end {
            break;
        }
        let dq = b.q_start as i64 - a.q_end as i64;
        let dt = b.t_start as i64 - a.t_end as i64;
        if (dq > 0 && dt > 0) || dt > MAX_GAP || dq > MAX_GAP {
            // Mirror UCSC chainToAxt/axtFromBlocks: the segment ends at the
            // last block before the split gap (a.t_end) and the next segment
            // starts at the following block (b.t_start); the gap itself is
            // skipped.
            if a.t_end > seg_start && a.t_end < t_end {
                segments.push((seg_start, a.t_end));
                seg_start = b.t_start;
            }
        }
    }
    segments.push((seg_start, t_end));

    for (s, e) in segments {
        emit_axt_segment(s, e, chain, t_2bit, q_2bit, matrix, writer, counter)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn emit_axt_segment<S: SequenceReader, W: Write>(
    t_start: u64,
    t_end: u64,
    chain: &Chain,
    t_2bit: &mut S,
    q_2bit: &mut S,
    matrix: &SubMatrix,
    writer: &mut W,
    counter: &mut usize,
) -> anyhow::Result<()> {
    // Get subset of chain
    let blocks = chain.to_blocks();

    // Find blocks overlapping [t_start, t_end)
    let mut idx_start = None;
    let mut idx_end = None;

    for (i, block) in blocks.iter().enumerate() {
        if block.t_end > t_start && block.t_start < t_end {
            if idx_start.is_none() {
                idx_start = Some(i);
            }
            idx_end = Some(i);
        }
    }

    if idx_start.is_none() {
        return Ok(());
    }

    if idx_start.is_none() {
        return Ok(());
    }

    let idx_start = idx_start.unwrap();
    let idx_end = idx_end.unwrap();

    let mut t_seq_all = String::new();
    let mut q_seq_all = String::new();

    // Helper to read Q sequence considering strand
    let read_q = |start: u64, end: u64, q_2bit: &mut S| -> anyhow::Result<String> {
        let (r_start, r_end) = if chain.header.q_strand == '-' {
            crate::libs::alignment::coords::reverse_range_pair(start, end, chain.header.q_size)
        } else {
            (start, end)
        };
        let mut seq = q_2bit.read_sequence(
            &chain.header.q_name,
            Some(r_start as usize),
            Some(r_end as usize),
        )?;
        if chain.header.q_strand == '-' {
            let rev = nt::rev_comp(seq.as_bytes()).collect();
            seq = String::from_utf8(rev)
                .map_err(|e| anyhow!("reverse-complemented query is not valid UTF-8: {}", e))?;
        }
        Ok(seq)
    };

    // Calculate initial q_start for the AXT record
    let q_start_out_base = if idx_start > 0 && t_start < blocks[idx_start].t_start {
        // We start in the gap before block[idx_start]
        let prev = &blocks[idx_start - 1];
        let gap_start_t = prev.t_end;
        if t_start <= gap_start_t {
            prev.q_end // We include dq
        } else {
            blocks[idx_start].q_start // prev.q_end + dq
        }
    } else {
        // We start in block[idx_start]
        let b = &blocks[idx_start];
        let offset = t_start.saturating_sub(b.t_start);
        b.q_start + offset
    };

    for i in idx_start..=idx_end {
        let block = &blocks[i];

        // 1. Handle gap BEFORE this block, but only when the previous block is
        // inside this segment (i > idx_start).  UCSC axtFromBlocks starts at
        // the first block of the segment, so gaps before idx_start are not
        // part of the AXT.
        if i > idx_start {
            let prev = &blocks[i - 1];
            // Gap range on T: [prev.t_end, block.t_start)
            let gap_start_t = prev.t_end;
            let gap_end_t = block.t_start;

            let overlap_start = gap_start_t.max(t_start);
            let overlap_end = gap_end_t.min(t_end);

            if overlap_start < overlap_end {
                // There is overlap with dt (T gap)
                let t_chunk = t_2bit.read_sequence(
                    &chain.header.t_name,
                    Some(overlap_start as usize),
                    Some(overlap_end as usize),
                )?;
                t_seq_all.push_str(&t_chunk);

                // Append Q dashes
                for _ in 0..(overlap_end - overlap_start) {
                    q_seq_all.push('-');
                }
            }

            // Handle dq (Q gap)
            if t_start <= gap_start_t && gap_start_t < t_end {
                let dq_len = block.q_start.saturating_sub(prev.q_end);
                if dq_len > 0 {
                    let q_chunk = read_q(prev.q_end, block.q_start, q_2bit)?;
                    q_seq_all.push_str(&q_chunk);

                    for _ in 0..dq_len {
                        t_seq_all.push('-');
                    }
                }
            }
        }

        // 2. Handle Block
        let start = block.t_start.max(t_start);
        let end = block.t_end.min(t_end);

        if start < end {
            let t_offset = start - block.t_start;
            let len = end - start;

            let t_chunk = t_2bit.read_sequence(
                &chain.header.t_name,
                Some(start as usize),
                Some(end as usize),
            )?;
            t_seq_all.push_str(&t_chunk);

            let q_start_seg = block.q_start + t_offset;
            let q_end_seg = q_start_seg + len;
            let q_chunk = read_q(q_start_seg, q_end_seg, q_2bit)?;
            q_seq_all.push_str(&q_chunk);
        }
    }

    // Calculate final q_end based on q_seq content (bases only)
    let q_bases_count = q_seq_all.chars().filter(|c| *c != '-').count() as u64;
    let q_end_out = q_start_out_base + q_bases_count;

    // Calculate score
    let score = calculate_score(&t_seq_all, &q_seq_all, matrix);

    writeln!(
        writer,
        "{} {} {} {} {} {} {} {} {}",
        *counter,
        chain.header.t_name,
        t_start + 1, // AXT 1-based
        t_end,
        chain.header.q_name,
        q_start_out_base + 1, // AXT 1-based
        q_end_out,
        chain.header.q_strand,
        score
    )?;
    *counter += 1;

    writeln!(writer, "{}", t_seq_all)?;
    writeln!(writer, "{}", q_seq_all)?;
    writeln!(writer)?;

    Ok(())
}

fn calculate_score(t_seq: &str, q_seq: &str, matrix: &SubMatrix) -> i32 {
    let mut score = 0;
    let t_chars: Vec<char> = t_seq.chars().collect();
    let q_chars: Vec<char> = q_seq.chars().collect();
    let len = t_chars.len();

    let mut in_gap_t = false;
    let mut in_gap_q = false;

    for i in 0..len {
        let t = t_chars[i];
        let q = q_chars[i];

        if t == '-' {
            // Gap in T (insertion in Q)
            if !in_gap_t {
                score -= matrix.gap_open;
                in_gap_t = true;
            }
            score -= matrix.gap_extend;
            in_gap_q = false;
        } else if q == '-' {
            // Gap in Q (deletion in Q)
            if !in_gap_q {
                score -= matrix.gap_open;
                in_gap_q = true;
            }
            score -= matrix.gap_extend;
            in_gap_t = false;
        } else {
            // Match/Mismatch
            score += matrix.get_score(t, q);
            in_gap_t = false;
            in_gap_q = false;
        }
    }
    score
}
