use std::cmp::Ordering;
use std::io::{Read, Seek};
use crate::libs::chaining::record::{Chain, ChainData, ChainHeader};
use crate::libs::chaining::algo::{ChainItem, KdTree};
use crate::libs::chaining::gap_calc::GapCalc;
use crate::libs::chaining::sub_matrix::SubMatrix;
use crate::libs::twobit::TwoBitFile;
use crate::libs::nt;

/// Represents a single alignment block that can be chained.
///
/// Contains coordinates in both target and query sequences, as well as a score.
#[derive(Clone, Debug)]
pub struct ChainableBlock {
    pub t_start: u64,
    pub t_end: u64,
    pub q_start: u64,
    pub q_end: u64,
    pub score: f64,
}

impl ChainItem for ChainableBlock {
    fn q_start(&self) -> u64 { self.q_start }
    fn q_end(&self) -> u64 { self.q_end }
    fn t_start(&self) -> u64 { self.t_start }
    fn t_end(&self) -> u64 { self.t_end }
    fn score(&self) -> f64 { self.score }
}

/// Context required for scoring chains based on actual sequence data.
///
/// Holds references to 2bit files for target and query, and a substitution matrix.
pub struct ScoreContext<'a, R> {
    pub t_2bit: &'a mut TwoBitFile<R>,
    pub q_2bit: &'a mut TwoBitFile<R>,
    pub matrix: &'a SubMatrix,
}

struct DpEntry {
    best_pred: Option<usize>,
    total_score: f64,
    hit: bool,
}

/// Chains a set of alignment blocks into optimal chains using a KD-tree and dynamic programming.
///
/// # Arguments
///
/// * `blocks` - A slice of `ChainableBlock`s to chain.
/// * `gap_calc` - Calculator for gap costs.
/// * `score_ctx` - Optional context for recalculating scores and trimming overlaps using sequence data.
/// * `q_name` - Query sequence name.
/// * `q_size` - Query sequence size.
/// * `q_strand` - Query sequence strand.
/// * `t_name` - Target sequence name.
/// * `t_size` - Target sequence size.
/// * `chain_id_counter` - Counter for generating unique chain IDs.
///
/// # Returns
///
/// A vector of `Chain`s sorted by score.
pub fn chain_blocks<R: Read + Seek>(
    blocks: &[ChainableBlock],
    gap_calc: &GapCalc,
    score_ctx: &mut Option<ScoreContext<R>>,
    q_name: &str,
    q_size: u64,
    q_strand: char,
    t_name: &str,
    t_size: u64,
    chain_id_counter: &mut usize,
) -> Vec<Chain> {
    if blocks.is_empty() {
        return Vec::new();
    }

    // 1. Create DP entries
    let mut dp_entries: Vec<DpEntry> = blocks
        .iter()
        .map(|b| DpEntry {
            best_pred: None,
            total_score: b.score,
            hit: false,
        })
        .collect();

    // 2. Build KD-tree
    let mut leaf_indices: Vec<usize> = (0..dp_entries.len()).collect();
    let mut tree = KdTree::build(&mut leaf_indices, blocks);

    // 3. Find best predecessors
    for i in 0..dp_entries.len() {
        let current_score = dp_entries[i].total_score;
        let cost_func = |cand_idx: usize, target_idx: usize| -> Option<f64> {
            let cand = &blocks[cand_idx];
            let target = &blocks[target_idx];
            
            // Ensure monotonic order for chaining
            if cand.t_start > target.t_start {
                return None;
            }
            // For q_start, strictly increasing is required for a linear chain
            if cand.q_start > target.q_start {
                return None;
            }

            let dt = target.t_start as i64 - cand.t_end as i64;
            let dq = target.q_start as i64 - cand.q_end as i64;
            
            let mut overlap_penalty = 0.0;

            // Handle overlaps (negative distance)
            if dt < 0 || dq < 0 {
                let ov_t = if dt < 0 { -dt } else { 0 };
                let ov_q = if dq < 0 { -dq } else { 0 };
                let overlap_len = std::cmp::max(ov_t, ov_q) as f64;
                
                // Estimate overlap penalty using score density
                // We use the maximum density of the two blocks as a conservative estimate
                let cand_len = (cand.t_end - cand.t_start) as f64;
                let target_len = (target.t_end - target.t_start) as f64;
                
                let cand_density = if cand_len > 0.0 { cand.score / cand_len } else { 0.0 };
                let target_density = if target_len > 0.0 { target.score / target_len } else { 0.0 };
                
                let density = cand_density.max(target_density);
                overlap_penalty = overlap_len * density;
            }

            let cost = gap_calc.calc(dq as i32, dt as i32) as f64;
            Some(dp_entries[cand_idx].total_score + target.score - cost - overlap_penalty)
        };
        let lower_bound_func = |dq: u64, dt: u64| -> f64 {
            gap_calc.calc(dq as i32, dt as i32) as f64
        };

        let (best_score, best_pred) = tree.best_predecessor(
            i,
            current_score,
            blocks,
            &cost_func,
            &lower_bound_func,
        );

        if best_score > dp_entries[i].total_score {
            dp_entries[i].total_score = best_score;
            dp_entries[i].best_pred = best_pred;
        }
        tree.update_scores(i, dp_entries[i].total_score, blocks);
    }

    // 4. Peel chains
    let mut sorted_indices: Vec<usize> = (0..dp_entries.len()).collect();
    sorted_indices.sort_by(|&a, &b| {
        dp_entries[b].total_score.partial_cmp(&dp_entries[a].total_score).unwrap_or(Ordering::Equal)
    });

    let mut chains = Vec::new();

    for &leaf_idx in &sorted_indices {
        if dp_entries[leaf_idx].hit {
            continue;
        }

        let mut chain_blocks_rev = Vec::new();
        let mut curr_idx = leaf_idx;
        
        loop {
            dp_entries[curr_idx].hit = true;
            chain_blocks_rev.push(blocks[curr_idx].clone());
            
            if let Some(pred_idx) = dp_entries[curr_idx].best_pred {
                curr_idx = pred_idx;
                if dp_entries[curr_idx].hit {
                    break; 
                }
            } else {
                break;
            }
        }

        chain_blocks_rev.reverse();
        
        remove_exact_overlaps(&mut chain_blocks_rev);
        merge_abutting_blocks(&mut chain_blocks_rev);

        // Trim overlaps if we have score context
        if let Some(ctx) = score_ctx {
             trim_overlaps(&mut chain_blocks_rev, ctx, q_name, t_name, q_size, q_strand);
        }

        let first = &chain_blocks_rev[0];
        let last = &chain_blocks_rev[chain_blocks_rev.len() - 1];

        // Recalculate score exactly
        let score = score_chain(&chain_blocks_rev, gap_calc, score_ctx, q_name, t_name, q_size, q_strand);

        if score <= 0.0 {
            continue; 
        }

        let mut chain_data = Vec::new();
        for i in 0..chain_blocks_rev.len() {
            let b = &chain_blocks_rev[i];
            let size = b.t_end - b.t_start;
            let (dt, dq) = if i < chain_blocks_rev.len() - 1 {
                let next = &chain_blocks_rev[i + 1];
                (next.t_start - b.t_end, next.q_start - b.q_end)
            } else {
                (0, 0)
            };
            chain_data.push(ChainData { size, dt, dq });
        }

        chains.push(Chain {
            header: ChainHeader {
                score,
                t_name: t_name.to_string(),
                t_size,
                t_strand: '+',
                t_start: first.t_start,
                t_end: last.t_end,
                q_name: q_name.to_string(),
                q_size,
                q_strand,
                q_start: first.q_start,
                q_end: last.q_end,
                id: *chain_id_counter as u64,
            },
            data: chain_data,
        });
        *chain_id_counter += 1;
    }

    chains
}

/// Trims overlaps between adjacent blocks in a chain using sequence data.
///
/// Adjusts the boundaries of overlapping blocks to maximize the score.
fn trim_overlaps<R: Read + Seek>(
    blocks: &mut Vec<ChainableBlock>,
    ctx: &mut ScoreContext<R>,
    q_name: &str,
    t_name: &str,
    q_size: u64,
    q_strand: char,
) {
    if blocks.len() < 2 {
        return;
    }
    
    let mut i = 0;
    while i < blocks.len() - 1 {
        let curr = &blocks[i];
        let next = &blocks[i+1];
        
        let overlap = if curr.t_end > next.t_start {
            (curr.t_end - next.t_start) as i64
        } else {
            0
        };
        
        if overlap > 0 {
            let overlap = overlap as usize;
            let (cut_pos, _) = find_crossover(
                &blocks[i], &blocks[i+1], overlap, ctx, q_name, t_name, q_size, q_strand
            );
            
            let trim_left = overlap as i64 - cut_pos as i64;
            let trim_right = cut_pos as i64;
            
            blocks[i].t_end -= trim_left as u64;
            blocks[i].q_end -= trim_left as u64;
            
            blocks[i+1].t_start += trim_right as u64;
            blocks[i+1].q_start += trim_right as u64;
        }
        i += 1;
    }
}

/// Finds the optimal crossover point for two overlapping blocks.
///
/// Returns the best cut position within the overlap and the score adjustment.
fn find_crossover<R: Read + Seek>(
    left: &ChainableBlock,
    right: &ChainableBlock,
    overlap: usize,
    ctx: &mut ScoreContext<R>,
    q_name: &str,
    t_name: &str,
    q_size: u64,
    q_strand: char,
) -> (usize, f64) {
    let l_t_seq = ctx.t_2bit.read_sequence(t_name, Some((left.t_end - overlap as u64) as usize), Some(left.t_end as usize), false).unwrap();
    let r_t_seq = ctx.t_2bit.read_sequence(t_name, Some(right.t_start as usize), Some((right.t_start + overlap as u64) as usize), false).unwrap();
    
    let l_q_seq = if q_strand == '+' {
        ctx.q_2bit.read_sequence(q_name, Some((left.q_end - overlap as u64) as usize), Some(left.q_end as usize), false).unwrap()
    } else {
        let start = (q_size - left.q_end) as usize;
        let end = (q_size - (left.q_end - overlap as u64)) as usize;
        let s = ctx.q_2bit.read_sequence(q_name, Some(start), Some(end), false).unwrap();
        let rc: Vec<u8> = nt::rev_comp(s.as_bytes()).collect();
        String::from_utf8(rc).unwrap()
    };
    
    let r_q_seq = if q_strand == '+' {
        ctx.q_2bit.read_sequence(q_name, Some(right.q_start as usize), Some((right.q_start + overlap as u64) as usize), false).unwrap()
    } else {
        let start = (q_size - (right.q_start + overlap as u64)) as usize;
        let end = (q_size - right.q_start) as usize;
        let s = ctx.q_2bit.read_sequence(q_name, Some(start), Some(end), false).unwrap();
        let rc: Vec<u8> = nt::rev_comp(s.as_bytes()).collect();
        String::from_utf8(rc).unwrap()
    };
    
    let mut best_pos = 0;
    let mut best_score = -1e9;
    
    let mut r_score = 0.0;
    let mut l_score = 0.0;
    
    let l_t_chars: Vec<char> = l_t_seq.chars().collect();
    let l_q_chars: Vec<char> = l_q_seq.chars().collect();
    let r_t_chars: Vec<char> = r_t_seq.chars().collect();
    let r_q_chars: Vec<char> = r_q_seq.chars().collect();
    
    for i in 0..overlap {
        l_score += ctx.matrix.get_score(l_t_chars[i], l_q_chars[i]) as f64;
        r_score += ctx.matrix.get_score(r_t_chars[i], r_q_chars[i]) as f64;
    }
    
    let mut current_l = 0.0;
    let mut current_r = r_score;
    
    for i in 0..=overlap {
        let score = current_l + current_r;
        if score > best_score {
            best_score = score;
            best_pos = i;
        }
        
        if i < overlap {
            current_l += ctx.matrix.get_score(l_t_chars[i], l_q_chars[i]) as f64;
            current_r -= ctx.matrix.get_score(r_t_chars[i], r_q_chars[i]) as f64;
        }
    }
    
    let adjustment = (r_score + l_score) - best_score;
    (best_pos, adjustment)
}

/// Removes duplicate blocks that have exact same coordinates.
fn remove_exact_overlaps(blocks: &mut Vec<ChainableBlock>) {
    if blocks.is_empty() {
        return;
    }
    
    let mut write_idx = 0;
    for read_idx in 1..blocks.len() {
        let is_duplicate = {
            let prev = &blocks[write_idx];
            let curr = &blocks[read_idx];
            curr.t_start == prev.t_start && curr.q_start == prev.q_start &&
            curr.t_end == prev.t_end && curr.q_end == prev.q_end
        };
        
        if is_duplicate {
            continue;
        }
        
        write_idx += 1;
        if write_idx != read_idx {
            blocks[write_idx] = blocks[read_idx].clone();
        }
    }
    blocks.truncate(write_idx + 1);
}

/// Merges adjacent blocks that abut perfectly.
fn merge_abutting_blocks(blocks: &mut Vec<ChainableBlock>) {
    if blocks.len() < 2 {
        return;
    }
    
    let mut write_idx = 0;
    for read_idx in 1..blocks.len() {
        let should_merge = {
            let prev = &blocks[write_idx];
            let curr = &blocks[read_idx];
            curr.t_start == prev.t_end && curr.q_start == prev.q_end
        };
        
        if should_merge {
             let curr_t_end = blocks[read_idx].t_end;
             let curr_q_end = blocks[read_idx].q_end;
             let curr_score = blocks[read_idx].score;
             
             blocks[write_idx].t_end = curr_t_end;
             blocks[write_idx].q_end = curr_q_end;
             blocks[write_idx].score += curr_score;
        } else {
            write_idx += 1;
            if write_idx != read_idx {
                blocks[write_idx] = blocks[read_idx].clone();
            }
        }
    }
    blocks.truncate(write_idx + 1);
}

/// Calculates the total score of a chain.
///
/// If `score_ctx` is provided, it recalculates block scores and gap costs using sequence data.
/// Otherwise, it uses the pre-calculated block scores and standard gap costs.
fn score_chain<R: Read + Seek>(
    blocks: &[ChainableBlock], 
    gap_calc: &GapCalc, 
    score_ctx: &mut Option<ScoreContext<R>>,
    q_name: &str,
    t_name: &str,
    q_size: u64,
    q_strand: char,
) -> f64 {
    let mut score = 0.0;
    for i in 0..blocks.len() {
        let block_score = if let Some(ctx) = score_ctx {
            calc_block_score(&blocks[i], ctx, q_name, t_name, q_size, q_strand).unwrap_or(0.0)
        } else {
            blocks[i].score
        };
        score += block_score;

        if i > 0 {
            let prev = &blocks[i - 1];
            let curr = &blocks[i];
            
            let dt = if curr.t_start >= prev.t_end {
                curr.t_start - prev.t_end
            } else {
                0 
            };
            let dq = if curr.q_start >= prev.q_end {
                curr.q_start - prev.q_end
            } else {
                0
            };
            
            if let Some(_ctx) = score_ctx {
                // If trimmed, dt >= 0.
                score -= gap_calc.calc(dq as i32, dt as i32) as f64;
            } else {
                score -= gap_calc.calc(dq as i32, dt as i32) as f64;
            }
        }
    }
    
    if let Some(ctx) = score_ctx {
         let mut exact_score = 0.0;
         for b in blocks {
             exact_score += calc_block_score(b, ctx, q_name, t_name, q_size, q_strand).unwrap_or(0.0);
         }
         for i in 1..blocks.len() {
             let prev = &blocks[i - 1];
             let curr = &blocks[i];
             let dt = (curr.t_start - prev.t_end) as i32;
             let dq = (curr.q_start - prev.q_end) as i32;
             exact_score -= gap_calc.calc(dq, dt) as f64;
         }
         return exact_score;
    }

    score
}

/// Calculates the score of a single block using sequence data and the substitution matrix.
pub fn calc_block_score<R: Read + Seek>(
    b: &ChainableBlock,
    ctx: &mut ScoreContext<R>,
    q_name: &str,
    t_name: &str,
    q_size: u64,
    q_strand: char,
) -> Option<f64> {
    let t_seq_res = ctx.t_2bit.read_sequence(
        t_name, 
        Some(b.t_start as usize), 
        Some(b.t_end as usize), 
        false
    );
    
    let q_seq_res = if q_strand == '+' {
        ctx.q_2bit.read_sequence(
            q_name, 
            Some(b.q_start as usize), 
            Some(b.q_end as usize), 
            false
        )
    } else {
        let start_pos = (q_size - b.q_end) as usize;
        let end_pos = (q_size - b.q_start) as usize;
        
        ctx.q_2bit.read_sequence(q_name, Some(start_pos), Some(end_pos), false)
            .map(|s| {
                let rc_bytes: Vec<u8> = nt::rev_comp(s.as_bytes()).collect();
                String::from_utf8(rc_bytes).unwrap()
            })
    };

    if let (Ok(t_seq), Ok(q_seq)) = (t_seq_res, q_seq_res) {
        let mut exact_score = 0.0;
        for (t, q) in t_seq.chars().zip(q_seq.chars()) {
            let val = ctx.matrix.get_score(t, q);
            exact_score += val as f64;
        }
        Some(exact_score)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_chain_blocks_basic() {
        // Gap costs in GapCalc::medium() are quite high (e.g., ~750 for length 10).
        // We need high block scores to justify chaining.
        let blocks = vec![
            ChainableBlock { t_start: 0, t_end: 10, q_start: 0, q_end: 10, score: 1000.0 },
            ChainableBlock { t_start: 20, t_end: 30, q_start: 20, q_end: 30, score: 1000.0 },
            ChainableBlock { t_start: 5, t_end: 15, q_start: 5, q_end: 15, score: 50.0 }, // Overlapping/conflicting
        ];

        let gap_calc = GapCalc::medium();
        // Use Cursor<Vec<u8>> as a dummy reader for ScoreContext since we pass None
        let mut score_ctx: Option<ScoreContext<Cursor<Vec<u8>>>> = None;
        let mut chain_id = 0;

        let chains = chain_blocks(
            &blocks,
            &gap_calc,
            &mut score_ctx,
            "chr1",
            100,
            '+',
            "chr1",
            100,
            &mut chain_id,
        );

        assert!(!chains.is_empty());
        let best_chain = &chains[0];
        
        // Should pick block 1 and block 2 (indices 0 and 1)
        // Score = 1000 + 1000 - gap_cost (~750) ≈ 1250
        assert!(best_chain.header.score > 1000.0);
        
        // Verify structure
        // ChainData stores blocks. 
        // Logic:
        // Block 1: size=10, dt=10, dq=10 (gap to next)
        // Block 2: size=10, dt=0, dq=0 (last block)
        assert_eq!(best_chain.data.len(), 2);
        assert_eq!(best_chain.data[0].size, 10);
        assert_eq!(best_chain.data[0].dt, 10);
        assert_eq!(best_chain.data[0].dq, 10);
        assert_eq!(best_chain.data[1].size, 10);
    }

    #[test]
    fn test_remove_exact_overlaps() {
        let mut blocks = vec![
            ChainableBlock { t_start: 0, t_end: 10, q_start: 0, q_end: 10, score: 100.0 },
            ChainableBlock { t_start: 0, t_end: 10, q_start: 0, q_end: 10, score: 100.0 },
            ChainableBlock { t_start: 20, t_end: 30, q_start: 20, q_end: 30, score: 100.0 },
        ];
        
        remove_exact_overlaps(&mut blocks);
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].t_start, 0);
        assert_eq!(blocks[1].t_start, 20);
    }

    #[test]
    fn test_merge_abutting_blocks() {
        let mut blocks = vec![
            ChainableBlock { t_start: 0, t_end: 10, q_start: 0, q_end: 10, score: 100.0 },
            ChainableBlock { t_start: 10, t_end: 20, q_start: 10, q_end: 20, score: 100.0 },
            ChainableBlock { t_start: 25, t_end: 30, q_start: 25, q_end: 30, score: 100.0 },
        ];
        
        merge_abutting_blocks(&mut blocks);
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].t_start, 0);
        assert_eq!(blocks[0].t_end, 20);
        assert_eq!(blocks[0].score, 200.0);
        assert_eq!(blocks[1].t_start, 25);
    }
}
