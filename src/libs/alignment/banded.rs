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
    let mut score = vec![0i32; (n + 1) * width];
    let mut trace = vec![0u8; (n + 1) * width]; // 0=reset, 1=diag, 2=up, 3=left
    let sub = |a: u8, b: u8| {
        if a == b {
            params.match_score
        } else {
            params.mismatch_score
        }
    };

    let mut best = 0i32;
    let mut best_off = 0usize;
    for i in 1..=n {
        // Only columns inside the diagonal band: |j - i - diag0| <= band.
        let j_lo = ((i as i64 + diag0 - band as i64).max(1).min(m as i64)) as usize;
        let j_hi = ((i as i64 + diag0 + band as i64).clamp(1, m as i64)) as usize;
        for j in j_lo..=j_hi {
            let off = (j as i64 - i as i64 - diag0 + band as i64) as usize;
            let c = i * width + off;
            let mut s = 0i32;
            let mut tr = 0u8;
            // Diagonal predecessor (i-1, j-1) sits on the same band offset.
            let v = score[(i - 1) * width + off] + sub(q[i - 1], t[j - 1]);
            if v > s {
                s = v;
                tr = 1;
            }
            // Up predecessor (i-1, j) shifts the offset by +1.
            if off + 1 < width {
                let v = score[(i - 1) * width + off + 1] + params.gap_open;
                if v > s {
                    s = v;
                    tr = 2;
                }
            }
            // Left predecessor (i, j-1) shifts the offset by -1.
            if off > 0 {
                let v = score[c - 1] + params.gap_open;
                if v > s {
                    s = v;
                    tr = 3;
                }
            }
            score[c] = s;
            trace[c] = tr;
            if s > best {
                best = s;
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
    while i > 0 && j > 0 {
        let d = j as i64 - i as i64 - diag0;
        let c = i * width + (d + band as i64) as usize;
        match trace[c] {
            0 => break,
            1 => {
                q_aln.push(q[i - 1]);
                t_aln.push(t[j - 1]);
                i -= 1;
                j -= 1;
            }
            2 => {
                q_aln.push(q[i - 1]);
                t_aln.push(b'-');
                i -= 1;
            }
            3 => {
                q_aln.push(b'-');
                t_aln.push(t[j - 1]);
                j -= 1;
            }
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
        // conserved flanks; the gapped solution (104) clearly beats the best
        // gapless one (60), so the gap must survive local trimming.
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
        assert_eq!(r.score, 24 * 5 + 2 * (-8));
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
