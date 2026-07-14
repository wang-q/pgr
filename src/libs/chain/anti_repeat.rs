//! Chain repeat and degeneracy filtering logic (ported from UCSC chainAntiRepeat).
//!
//! Note: `nt_val` here intentionally uses a T=0,C=1,A=2,G=3 mapping (with -1 for
//! invalid) rather than `crate::libs::nt::NT_VAL` (A=0,C=1,G=2,T=3, 255 for
//! invalid), because the reverse-strand complement logic `(v + 2) % 4` relies on
//! this specific ordering (T<->A, C<->G). Switching to `NT_VAL` would break the
//! complement calculation.

use crate::libs::alignment::coords::reverse_range_pair;
use crate::libs::chain::{Block, Chain};
use crate::libs::ds::TopKPurity;
use crate::libs::fmt::twobit::TwoBitFile;
use crate::libs::nt::is_lower;

/// Check a chain against both the degeneracy and repeat filters.
pub fn check_chain<R: std::io::Read + std::io::Seek>(
    chain: &Chain,
    t_2bit: &mut TwoBitFile<R>,
    q_2bit: &mut TwoBitFile<R>,
    min_score: f64,
) -> bool {
    // Check if sequences exist
    if !t_2bit.sequence_offsets.contains_key(&chain.header.t_name)
        || !q_2bit.sequence_offsets.contains_key(&chain.header.q_name)
    {
        return false;
    }

    let blocks = chain.to_blocks();

    // 1. Degeneracy Filter (Low complexity check)
    if !check_degeneracy(chain, &blocks, t_2bit, q_2bit, min_score) {
        return false;
    }

    // 2. Repeat Filter (Lowercase check)
    check_repeat(chain, &blocks, t_2bit, q_2bit, min_score)
}

/// Read target and query slices for a block, accounting for query strand.
pub fn get_slices<R: std::io::Read + std::io::Seek>(
    block: &Block,
    t_2bit: &mut TwoBitFile<R>,
    q_2bit: &mut TwoBitFile<R>,
    t_name: &str,
    q_name: &str,
    q_strand: char,
    q_size: u64,
) -> Option<(Vec<u8>, Vec<u8>)> {
    // Read Target Slice
    let t_seq = match t_2bit.read_sequence(
        t_name,
        Some(block.t_start as usize),
        Some(block.t_end as usize),
        false, // include soft masks
    ) {
        Ok(s) => s.into_bytes(),
        Err(_) => return None,
    };

    // Calculate Query Range
    let (q_start, q_end) = if q_strand == '+' {
        (block.q_start as usize, block.q_end as usize)
    } else {
        let (s, e) = reverse_range_pair(block.q_start as i32, block.q_end as i32, q_size as i32);
        (s as usize, e as usize)
    };

    // Read Query Slice
    let q_seq = match q_2bit.read_sequence(
        q_name,
        Some(q_start),
        Some(q_end),
        false, // include soft masks
    ) {
        Ok(s) => s.into_bytes(),
        Err(_) => return None,
    };

    Some((t_seq, q_seq))
}

/// Low-complexity (degeneracy) filter: penalize chains dominated by 1-2 bases.
pub fn check_degeneracy<R: std::io::Read + std::io::Seek>(
    chain: &Chain,
    blocks: &[Block],
    t_2bit: &mut TwoBitFile<R>,
    q_2bit: &mut TwoBitFile<R>,
    min_score: f64,
) -> bool {
    let mut detector = TopKPurity::new(4, 2, 0.80);

    for block in blocks {
        if let Some((t_slice, q_slice)) = get_slices(
            block,
            t_2bit,
            q_2bit,
            &chain.header.t_name,
            &chain.header.q_name,
            chain.header.q_strand,
            chain.header.q_size,
        ) {
            for i in 0..t_slice.len() {
                let t_base = t_slice[i];
                let q_base_raw = if chain.header.q_strand == '+' {
                    q_slice[i]
                } else {
                    q_slice[q_slice.len() - 1 - i]
                };

                let t_val = nt_val(t_base);
                let mut q_val = nt_val(q_base_raw);

                if chain.header.q_strand == '-' && q_val >= 0 {
                    q_val = (q_val + 2) % 4;
                }

                if t_val >= 0 && t_val == q_val {
                    detector.increment(t_val as usize);
                }
            }
        }
    }

    if detector.total() == 0 {
        return false;
    }

    match detector.penalty_factor() {
        None => true,
        Some(factor) => {
            let adjusted_score = chain.header.score * factor;
            if adjusted_score < min_score {
                log::info!(
                    "Chain {} filtered by degeneracy: score {} -> {}",
                    chain.header.id,
                    chain.header.score,
                    adjusted_score
                );
                false
            } else {
                true
            }
        }
    }
}

/// Repeat filter: penalize chains with many soft-masked (lowercase) bases.
pub fn check_repeat<R: std::io::Read + std::io::Seek>(
    chain: &Chain,
    blocks: &[Block],
    t_2bit: &mut TwoBitFile<R>,
    q_2bit: &mut TwoBitFile<R>,
    min_score: f64,
) -> bool {
    let mut rep_count = 0;
    let mut total = 0;

    for block in blocks {
        if let Some((t_slice, q_slice)) = get_slices(
            block,
            t_2bit,
            q_2bit,
            &chain.header.t_name,
            &chain.header.q_name,
            chain.header.q_strand,
            chain.header.q_size,
        ) {
            for i in 0..t_slice.len() {
                let t_base = t_slice[i];
                let q_base = if chain.header.q_strand == '+' {
                    q_slice[i]
                } else {
                    q_slice[q_slice.len() - 1 - i]
                };

                if is_lower(t_base) || is_lower(q_base) {
                    rep_count += 1;
                }
            }
            total += t_slice.len();
        }
    }

    if total == 0 {
        return false;
    }

    let adjusted_score = chain.header.score * 2.0 * ((total - rep_count) as f64) / (total as f64);
    if adjusted_score < min_score {
        log::info!(
            "Chain {} filtered by repeat: score {} -> {} (rep {}/{})",
            chain.header.id,
            chain.header.score,
            adjusted_score,
            rep_count,
            total
        );
        false
    } else {
        true
    }
}

/// Map a nucleotide byte to T=0, C=1, A=2, G=3, or -1 for non-ACGT.
///
/// This mapping is intentionally distinct from `crate::libs::nt::NT_VAL` because
/// the reverse-strand complement logic in `check_degeneracy` relies on the
/// `(v + 2) % 4` identity (T<->A, C<->G) that this ordering provides.
pub fn nt_val(base: u8) -> i8 {
    match base {
        b'T' | b't' => 0,
        b'C' | b'c' => 1,
        b'A' | b'a' => 2,
        b'G' | b'g' => 3,
        _ => -1,
    }
}
