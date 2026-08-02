//! Banded local pairwise alignment for diagonal-band intervals.

use crate::libs::poa::align::AlignmentParams;

/// Result of a banded local alignment: aligned strings plus the offset of
/// the first aligned base within each input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BandedAlign {
    /// Alignment score (sum of match/mismatch/gap penalties).
    pub score: i32,
    /// Query bases with `-` for gaps.
    pub q_aln: Vec<u8>,
    /// Target bases with `-` for gaps (same length as `q_aln`).
    pub t_aln: Vec<u8>,
    /// Offset of the first aligned query base.
    pub q_start: usize,
    /// Offset of the first aligned target base.
    pub t_start: usize,
}

/// Banded Smith-Waterman local alignment restricted to the diagonal band
/// `|(i - j) - diag0| <= band` (i/j are query/target positions).
///
/// Returns `None` when the best local score is non-positive. The band keeps
/// the DP linear in the sequence length (O(n * band)) instead of O(n * m),
/// which is required for whole-chain intervals in large genomes.
pub fn align_banded_local(
    q: &[u8],
    t: &[u8],
    band: usize,
    diag0: i64,
    params: &AlignmentParams,
) -> Option<BandedAlign> {
    let n = q.len();
    let m = t.len();
    if n == 0 || m == 0 {
        return None;
    }
    let width = 2 * band + 1;
    let neg = i32::MIN / 4;
    // Affine-gap local alignment (M/I/D states). M resets to 0 (local);
    // I = insertion in q (q base vs gap), D = deletion (gap vs t base).
    let mut mscore = vec![0i32; (n + 1) * width];
    let mut iscore = vec![neg; (n + 1) * width];
    let mut dscore = vec![neg; (n + 1) * width];
    // Packed traces, 2 bits per state:
    //   M: 0=reset, 1=from M diag, 2=from I diag, 3=from D diag
    //   I: 0=none, 1=open from M, 2=extend from I
    //   D: 0=none, 1=open from M, 2=extend from D
    let mut trace = vec![0u8; (n + 1) * width];
    let sub = |a: u8, b: u8| {
        if a == b {
            params.match_score
        } else {
            params.mismatch_score
        }
    };
    let go = params.gap_open;
    let ge = params.gap_extend;

    let mut best = 0i32;
    let mut best_off = 0usize;
    for i in 1..=n {
        // Only columns inside the diagonal band: |j - i - diag0| <= band.
        let j_lo = ((i as i64 + diag0 - band as i64).max(1).min(m as i64)) as usize;
        let j_hi = ((i as i64 + diag0 + band as i64).clamp(1, m as i64)) as usize;
        for j in j_lo..=j_hi {
            let off = (j as i64 - i as i64 - diag0 + band as i64) as usize;
            let c = i * width + off;
            // I: insertion in q, predecessor (i-1, j) at off+1.
            let mut iv = neg;
            let mut it = 0u8;
            if off + 1 < width {
                let up = (i - 1) * width + off + 1;
                let open = mscore[up].saturating_add(go);
                let ext = iscore[up].saturating_add(ge);
                if open >= ext {
                    iv = open;
                    it = 1;
                } else {
                    iv = ext;
                    it = 2;
                }
            }
            // D: deletion, predecessor (i, j-1) at off-1.
            let mut dv = neg;
            let mut dt = 0u8;
            if off > 0 {
                let left = c - 1;
                let open = mscore[left].saturating_add(go);
                let ext = dscore[left].saturating_add(ge);
                if open >= ext {
                    dv = open;
                    dt = 1;
                } else {
                    dv = ext;
                    dt = 2;
                }
            }
            // M: diagonal predecessor (i-1, j-1) at the same offset; local
            // reset to 0. Tie-break M > I > D for deterministic traceback.
            let diag = (i - 1) * width + off;
            let mut cand = mscore[diag];
            let mut st = 1u8;
            if iscore[diag] > cand {
                cand = iscore[diag];
                st = 2;
            }
            if dscore[diag] > cand {
                cand = dscore[diag];
                st = 3;
            }
            let v = cand.saturating_add(sub(q[i - 1], t[j - 1]));
            let (mv, mt) = if v > 0 { (v, st) } else { (0, 0) };
            mscore[c] = mv;
            iscore[c] = iv;
            dscore[c] = dv;
            trace[c] = mt | (it << 2) | (dt << 4);
            if mv > best {
                best = mv;
                best_off = c;
            }
        }
    }
    if best <= 0 {
        return None;
    }

    // Recover (i, j) from the flat offset: i = offset / width; the stored
    // offset formula gives the j back via the band relation.
    let bi = best_off / width;
    let bo = best_off % width;
    let bj = (bo as i64 + bi as i64 + diag0 - band as i64) as usize;

    let mut q_aln = Vec::with_capacity(n);
    let mut t_aln = Vec::with_capacity(m);
    let (mut i, mut j) = (bi, bj);
    let mut state = 0u8; // 0=M, 1=I, 2=D
    while i > 0 && j > 0 {
        let d = j as i64 - i as i64 - diag0;
        let c = i * width + (d + band as i64) as usize;
        let tr = trace[c];
        match state {
            0 => match tr & 3 {
                0 => break,
                1 => {
                    q_aln.push(q[i - 1]);
                    t_aln.push(t[j - 1]);
                    i -= 1;
                    j -= 1;
                    state = 0;
                }
                2 => {
                    q_aln.push(q[i - 1]);
                    t_aln.push(t[j - 1]);
                    i -= 1;
                    j -= 1;
                    state = 1;
                }
                3 => {
                    q_aln.push(q[i - 1]);
                    t_aln.push(t[j - 1]);
                    i -= 1;
                    j -= 1;
                    state = 2;
                }
                _ => unreachable!(),
            },
            1 => match (tr >> 2) & 3 {
                1 => {
                    q_aln.push(q[i - 1]);
                    t_aln.push(b'-');
                    i -= 1;
                    state = 0;
                }
                2 => {
                    q_aln.push(q[i - 1]);
                    t_aln.push(b'-');
                    i -= 1;
                    state = 1;
                }
                _ => break,
            },
            2 => match (tr >> 4) & 3 {
                1 => {
                    q_aln.push(b'-');
                    t_aln.push(t[j - 1]);
                    j -= 1;
                    state = 0;
                }
                2 => {
                    q_aln.push(b'-');
                    t_aln.push(t[j - 1]);
                    j -= 1;
                    state = 2;
                }
                _ => break,
            },
            _ => unreachable!(),
        }
    }
    q_aln.reverse();
    t_aln.reverse();
    Some(BandedAlign {
        score: best,
        q_start: i,
        t_start: j,
        q_aln,
        t_aln,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(q: &[u8], t: &[u8], band: usize, diag0: i64) -> BandedAlign {
        align_banded_local(q, t, band, diag0, &AlignmentParams::default()).unwrap()
    }

    #[test]
    fn perfect_match() {
        let r = aln(b"ACGTACGT", b"ACGTACGT", 4, 0);
        assert_eq!(r.q_aln, b"ACGTACGT");
        assert_eq!(r.t_aln, b"ACGTACGT");
        assert_eq!((r.q_start, r.t_start), (0, 0));
        assert_eq!(r.score, 8 * 5);
    }

    #[test]
    fn internal_insertion_is_kept() {
        // q has an internal 2-base insertion ("TT") between two 12-base
        // conserved flanks; the affine-gapped solution (106 = 24*5 - 8 - 6)
        // clearly beats the best gapless one (60), so the gap must survive
        // local trimming.
        let mut q = b"ACGTACGTACGT".to_vec();
        q.extend_from_slice(b"TT");
        q.extend_from_slice(b"ACGTACGTACGT");
        let t = b"ACGTACGTACGTACGTACGTACGT";
        let r = aln(&q, t, 4, 0);
        assert_eq!(r.q_aln.len(), r.t_aln.len());
        assert!(
            r.t_aln.contains(&b'-'),
            "insertion gap expected: {:?}",
            r.t_aln
        );
        assert_eq!(r.q_aln, q);
        assert_eq!(r.t_aln.len(), t.len() + 2);
        assert_eq!(r.score, 24 * 5 - 8 - 6);
    }

    #[test]
    fn local_trimming_finds_core() {
        // Non-homologous flanks; only the middle TT pair aligns.
        let q = b"CCCTTT";
        let t = b"AAATT";
        let r = aln(q, t, 4, 0);
        assert_eq!(r.q_aln, b"TT");
        assert_eq!(r.t_aln, b"TT");
        assert_eq!(r.q_start, 3);
        assert_eq!(r.t_start, 3);
    }

    #[test]
    fn band_excludes_off_diagonal() {
        // q and t share no k-mer on the main diagonal; with band 0 nothing
        // aligns, but a shifted band (diag0 = t_pos - q_pos = 2) finds the
        // shifted match q[0..6) vs t[2..8) = "ACGTAC".
        let q = b"ACGTACGT";
        let t = b"TTACGTAC";
        let none = align_banded_local(q, t, 0, 0, &AlignmentParams::default());
        assert!(none.is_none() || none.unwrap().score <= 0);
        let shifted = align_banded_local(q, t, 2, 2, &AlignmentParams::default());
        let r = shifted.expect("shifted band must align");
        assert_eq!(r.q_aln, b"ACGTAC");
        assert_eq!(r.t_aln, b"ACGTAC");
        assert_eq!((r.q_start, r.t_start), (0, 2));
    }

    #[test]
    fn no_homology_returns_none() {
        let q = b"AAAAAAAA";
        let t = b"CCCCCCCC";
        assert!(align_banded_local(q, t, 4, 0, &AlignmentParams::default()).is_none());
    }
}
