//! Myers wavefront extension, ported from FastGA `align.c` `forward_wave`.
//!
//! A local alignment is extended from an anchor point in both directions
//! using the unit-cost edit-distance wavefront (V[k] = furthest reaching
//! anti-diagonal on diagonal k), with the three-branch update and match-snake
//! from the original algorithm. The wave expands one diagonal per edit from
//! the anchor (WFA-style), and the exact path is reconstructed from the
//! per-wave predecessor trace.

/// Result of a bidirectional wave extension from an anchor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WaveAlign {
    /// Number of matched bases.
    pub matches: usize,
    /// Query bases with `-` for gaps.
    pub q_aln: Vec<u8>,
    /// Target bases with `-` for gaps (same length as `q_aln`).
    pub t_aln: Vec<u8>,
    /// Offset of the first aligned query base.
    pub q_start: usize,
    /// Offset of the first aligned target base.
    pub t_start: usize,
}

/// One wavefront trace entry: predecessor diagonal and the query position
/// right after the edit (before the match snake).
#[derive(Clone, Copy, Default)]
struct Trace {
    pred: i64,
    pre_x: i64,
}

const TRIM_MLAG: i64 = 250; // FastGA: stop when the best lags this many edits

/// Forward wavefront from an anchor point toward increasing coordinates.
///
/// Returns the wave history (V per wave/diagonal, traces, best point) or
/// `None` when nothing extends.
fn forward_wave(q: &[u8], t: &[u8], anchor: (i64, i64), band: usize) -> Option<WaveHistory> {
    let n = q.len() as i64;
    let m = t.len() as i64;
    let (aq, at) = anchor;
    let k0 = aq - at;
    let k_lo = (k0 - band as i64).max(-m);
    let k_hi = (k0 + band as i64).min(n);
    if k_lo > k_hi {
        return None;
    }
    let width = (k_hi - k_lo + 1) as usize;
    let d_cap = (n + m) as usize + 8;

    let mut history = WaveHistory {
        v: vec![-1i64; d_cap * width],
        trace: vec![Trace::default(); d_cap * width],
        width,
        k_lo,
        best_d: 0,
        best_k: k0,
        best_c: -1,
    };
    let off = |k: i64| (k - k_lo) as usize;

    // 0-wave: snake from the anchor point on its diagonal.
    let mut x = aq;
    while x < n && x - k0 < m && q[x as usize] == t[(x - k0) as usize] {
        x += 1;
    }
    let c = (x << 1) - k0;
    history.v[off(k0)] = c;
    history.best_c = c;
    if c <= aq + at {
        return None; // no match even at the anchor itself
    }

    // Successive waves: expand the diagonal range by one per edit.
    let mut besta = c;
    let mut lasta = c;
    let mut d = 1usize;
    while lasta >= besta - TRIM_MLAG {
        let lo = (k0 - d as i64).max(k_lo);
        let hi = (k0 + d as i64).min(k_hi);
        if lo > hi {
            break;
        }
        let prev = (d - 1) * width;
        let mut new_besta = -1i64;
        for k in lo..=hi {
            let am = if k > k_lo {
                history.v[prev + off(k - 1)]
            } else {
                -1
            };
            let ac = history.v[prev + off(k)];
            let ap = if k < k_hi {
                history.v[prev + off(k + 1)]
            } else {
                -1
            };
            // Three-branch update; ties prefer ap > am > ac (FastGA order).
            let (cand, pred) = if ac < am {
                if am < ap {
                    (ap + 1, k + 1)
                } else {
                    (am + 1, k - 1)
                }
            } else if ac < ap {
                (ap + 1, k + 1)
            } else {
                (ac + 2, k)
            };
            if cand < 0 {
                continue;
            }
            let mut x = (cand + k) >> 1;
            if x > n {
                x = n;
            }
            if x - k > m {
                x = m + k;
            }
            if x < 0 || x - k < 0 {
                continue;
            }
            let pre_x = x;
            while x < n && x - k < m && q[x as usize] == t[(x - k) as usize] {
                x += 1;
            }
            let cf = (x << 1) - k;
            let cell = d * width + off(k);
            history.v[cell] = cf;
            history.trace[cell] = Trace { pred, pre_x };
            if cf > new_besta {
                new_besta = cf;
            }
            if cf > besta {
                besta = cf;
                history.best_d = d;
                history.best_k = k;
                history.best_c = cf;
            }
        }
        lasta = new_besta;
        d += 1;
        if d >= d_cap {
            break;
        }
    }
    if history.best_c < 0 {
        return None;
    }
    Some(history)
}

/// Wave history: V per (wave, diagonal) plus the best point.
struct WaveHistory {
    v: Vec<i64>,
    trace: Vec<Trace>,
    width: usize,
    k_lo: i64,
    best_d: usize,
    best_k: i64,
    best_c: i64,
}

/// Reconstruct the path from the anchor to the best point as operations.
///
/// Returns operations in anchor-first order as `(query index, target index)`
/// pairs with `None` for a gap.
fn traceback_forward(anchor: (i64, i64), history: &WaveHistory) -> Vec<(Option<i64>, Option<i64>)> {
    let width = history.width;
    let k_lo = history.k_lo;
    let off = |k: i64| (k - k_lo) as usize;
    let mut ops: Vec<(Option<i64>, Option<i64>)> = Vec::new();
    let mut d = history.best_d;
    let mut k = history.best_k;
    while d > 0 {
        let tr = history.trace[d * width + off(k)];
        let cf = history.v[d * width + off(k)];
        let x_end = (cf + k) >> 1;
        for x in tr.pre_x..x_end {
            ops.push((Some(x), Some(x - k)));
        }
        let pre_x = tr.pre_x;
        let pre_y = pre_x - k;
        if tr.pred == k - 1 {
            ops.push((Some(pre_x - 1), None));
        } else if tr.pred == k + 1 {
            ops.push((None, Some(pre_y - 1)));
        } else {
            ops.push((Some(pre_x - 1), Some(pre_y - 1)));
        }
        k = tr.pred;
        d -= 1;
    }
    // 0-wave snake from the anchor (on diagonal k0) to the traceback end.
    let k0 = anchor.0 - anchor.1;
    debug_assert_eq!(k, k0);
    let c0 = history.v[off(k0)];
    let x0 = (c0 + k0) >> 1;
    for x in anchor.0..x0 {
        ops.push((Some(x), Some(x - k0)));
    }
    ops
}

/// Bidirectional extension from an anchor point.
///
/// The forward half extends toward the sequence ends; the reverse half runs
/// the same wave on mirrored sequences and is converted back, so both halves
/// meet exactly at the anchor.
pub fn wave_extend(q: &[u8], t: &[u8], band: usize, anchor: (usize, usize)) -> Option<WaveAlign> {
    let n = q.len() as i64;
    let m = t.len() as i64;
    let (aq, at) = (anchor.0 as i64, anchor.1 as i64);
    if aq >= n || at >= m || q[aq as usize] != t[at as usize] {
        return None;
    }

    let fwd = forward_wave(q, t, (aq, at), band)?;
    let fwd_ops = traceback_forward((aq, at), &fwd);

    // Reverse half on mirrored sequences.
    let rq: Vec<u8> = q.iter().rev().copied().collect();
    let rt: Vec<u8> = t.iter().rev().copied().collect();
    let (raq, rat) = (n - 1 - aq, m - 1 - at);
    let rev = forward_wave(&rq, &rt, (raq, rat), band)?;
    let rev_ops = traceback_forward((raq, rat), &rev);

    // Convert the reverse half to forward coordinates (drop the anchor op,
    // which is duplicated by the forward half) and reverse it.
    let mut pre: Vec<(Option<i64>, Option<i64>)> = rev_ops
        .into_iter()
        .skip(1)
        .map(|(qi, ti)| (qi.map(|x| n - 1 - x), ti.map(|y| m - 1 - y)))
        .collect();
    pre.reverse();

    let mut q_aln: Vec<u8> = Vec::new();
    let mut t_aln: Vec<u8> = Vec::new();
    let mut matches = 0usize;
    let mut q_start = aq as usize;
    let mut t_start = at as usize;
    for (qi, ti) in pre.iter().chain(fwd_ops.iter()) {
        match (qi, ti) {
            (Some(x), Some(y)) => {
                q_aln.push(q[*x as usize]);
                t_aln.push(t[*y as usize]);
                if q[*x as usize] == t[*y as usize] {
                    matches += 1;
                }
                q_start = q_start.min(*x as usize);
                t_start = t_start.min(*y as usize);
            }
            (Some(x), None) => {
                q_aln.push(q[*x as usize]);
                t_aln.push(b'-');
                q_start = q_start.min(*x as usize);
            }
            (None, Some(y)) => {
                q_aln.push(b'-');
                t_aln.push(t[*y as usize]);
                t_start = t_start.min(*y as usize);
            }
            (None, None) => unreachable!(),
        }
    }

    Some(WaveAlign {
        matches,
        q_aln,
        t_aln,
        q_start,
        t_start,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_align_fully() {
        let q = b"ACGTACGTACGT";
        let t = b"ACGTACGTACGT";
        let r = wave_extend(q, t, 4, (6, 6)).unwrap();
        assert_eq!(r.q_aln, q);
        assert_eq!(r.t_aln, t);
        assert_eq!(r.matches, q.len());
        assert_eq!((r.q_start, r.t_start), (0, 0));
    }

    #[test]
    fn anchor_mismatch_returns_none() {
        let q = b"ACGTACGTACGT";
        let t = b"ACGTACGTACGT";
        assert!(wave_extend(q, t, 4, (0, 2)).is_none());
    }

    #[test]
    fn one_mismatch_keeps_one_edit() {
        let q = b"ACGTTCGTACGT";
        let t = b"ACGTACGTACGT";
        // Anchor at a matching base near the mismatch.
        let r = wave_extend(q, t, 4, (6, 6)).unwrap();
        assert_eq!(r.q_aln.len(), r.t_aln.len());
        assert_eq!(r.matches + 1, r.q_aln.len(), "exactly one edit");
        // The single edit is a substitution (no gaps).
        assert_eq!(
            r.q_aln
                .iter()
                .zip(&r.t_aln)
                .filter(|(a, b)| a != b && **a != b'-' && **b != b'-')
                .count(),
            1
        );
    }

    #[test]
    fn internal_insertion_produces_valid_path() {
        let mut q = b"ACGTACGTACGT".to_vec();
        q.extend_from_slice(b"TTTT");
        q.extend_from_slice(b"ACGTACGTACGT");
        let t = b"ACGTACGTACGTACGTACGTACGT";
        // Unit-cost edits: a 4-base insertion (4 edits) ties with 4
        // mismatches, so only assert the path is valid and short.
        let r = wave_extend(&q, t, 4, (11, 11)).unwrap();
        assert_eq!(r.q_aln.len(), r.t_aln.len());
        let edits = r.q_aln.len() - r.matches;
        assert!(edits <= 4, "too many edits: {edits}");
        assert!(edits >= 1);
    }
}
