//! Optimal cut-point search for overlapping alignment blocks.

/// Find the best cut position inside an overlap between two aligned sequence
/// pairs, mirroring UCSC `cBlockFindCrossover`.
///
/// For each cut `i` in `0..=overlap`, the score is the left prefix of the
/// first pair plus the right suffix of the second pair; the cut maximizing
/// that score is returned together with the overlap adjustment
/// (`r_score + l_score - best_score`), which is what UCSC subtracts from the
/// chained score for the overlapping region.
///
/// All four slices must have the same length (the overlap size).
pub fn best_crossover(
    l_t: &[u8],
    l_q: &[u8],
    r_t: &[u8],
    r_q: &[u8],
    score: impl Fn(u8, u8) -> f64,
) -> (usize, f64) {
    let overlap = l_t.len();
    debug_assert_eq!(l_q.len(), overlap);
    debug_assert_eq!(r_t.len(), overlap);
    debug_assert_eq!(r_q.len(), overlap);

    let mut best_pos = 0;
    let mut best_score = f64::NEG_INFINITY;

    let mut r_score = 0.0;
    for i in 0..overlap {
        r_score += score(r_t[i], r_q[i]);
    }

    let mut current_l = 0.0;
    let mut current_r = r_score;

    for i in 0..=overlap {
        let total = current_l + current_r;
        if total > best_score {
            best_score = total;
            best_pos = i;
        }

        if i < overlap {
            current_l += score(l_t[i], l_q[i]);
            current_r -= score(r_t[i], r_q[i]);
        }
    }

    let overlap_adjustment = r_score + current_l - best_score;
    (best_pos, overlap_adjustment)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn score(a: u8, b: u8) -> f64 {
        if a == b {
            1.0
        } else {
            -1.0
        }
    }

    #[test]
    fn test_cut_prefers_left_prefix_and_right_suffix() {
        // Left pair matches only in the first half, right pair only in the
        // second half; the best cut should fall in the middle.
        let l_t = b"AAAAACCCCC";
        let l_q = b"AAAAAGGGGG";
        let r_t = b"CCCCCGGGGG";
        let r_q = b"TTTTTGGGGG";
        let (pos, _) = best_crossover(l_t, l_q, r_t, r_q, score);
        assert_eq!(pos, 5);
    }

    #[test]
    fn test_cut_at_edges() {
        // Left pair always better: best cut at the end (all overlap kept by
        // the left pair).
        let l_t = b"AAAAAAAAAA";
        let l_q = b"AAAAAAAAAA";
        let r_t = b"CCCCCCCCCC";
        let r_q = b"GGGGGGGGGG";
        let (pos, _) = best_crossover(l_t, l_q, r_t, r_q, score);
        assert_eq!(pos, 10);
    }

    #[test]
    fn test_adjustment_is_r_plus_l_minus_best() {
        let l_t = b"ACGT";
        let l_q = b"ACGT";
        let r_t = b"TGCG";
        let r_q = b"TGCG";
        let (pos, adj) = best_crossover(l_t, l_q, r_t, r_q, score);
        let _ = pos;
        let l_score = (0..4).map(|i| score(l_t[i], l_q[i])).sum::<f64>();
        let r_score = (0..4).map(|i| score(r_t[i], r_q[i])).sum::<f64>();
        let best = (0..=4)
            .map(|i| {
                (0..i).map(|k| score(l_t[k], l_q[k])).sum::<f64>()
                    + (i..4).map(|k| score(r_t[k], r_q[k])).sum::<f64>()
            })
            .fold(f64::NEG_INFINITY, f64::max);
        assert_eq!(adj, r_score + l_score - best);
    }
}
