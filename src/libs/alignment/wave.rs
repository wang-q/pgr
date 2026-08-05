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

/// One edit of the Myers divide-and-conquer edit script.
///
/// Matches are implicit between consecutive ops; `q_pos`/`t_pos` are absolute
/// coordinates in the original sequences.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EditOp {
    /// Delete one query base at `q_pos` (the target has no base there).
    Del { q_pos: usize, t_pos: usize },
    /// Insert one target base at `t_pos` (the query has no base there).
    Ins { q_pos: usize, t_pos: usize },
}

/// Myers `split_nd`: edit distance of `q[0..m)` vs `t[0..n)` and a midpoint
/// `(x, y)` on an optimal path (port of FastGA `align.c:5046`).
///
/// A substitution counts as one edit (FastGA's metric; see the wave update
/// below, where the same-diagonal candidate costs +1).
fn split_nd(q: &[u8], t: &[u8]) -> (usize, usize, usize) {
    let m = q.len() as i64;
    let n = t.len() as i64;
    let mut y = 0i64;
    if n < m {
        while y < n && t[y as usize] == q[y as usize] {
            y += 1;
        }
    } else {
        while y < m && t[y as usize] == q[y as usize] {
            y += 1;
        }
        if y >= m && n == m {
            return (0, m as usize, m as usize);
        }
    }

    // Diagonal k = x - y; VF/VB hold the furthest y on each diagonal. The
    // backward wave reaches k up to 2n-1 for insertion-heavy spans, so the
    // arrays span k in [-(2(m+n)+4), 2(m+n)+3] with pointer-style indexing.
    let off = 2 * (m + n) + 4;
    let size = 2 * off as usize;
    let mut vf = vec![-2i64; size];
    let mut vb = vec![0i64; size];
    vf[off as usize] = y;
    vf[(off - 1) as usize] = -2;

    // Backward wave init: leading match run from the far end.
    let x0 = n - m;
    let mut yb = n - 1;
    if n > m {
        while yb >= x0 && t[yb as usize] == q[(yb - x0) as usize] {
            yb -= 1;
        }
    } else {
        while yb >= 0 && t[yb as usize] == q[(yb - x0) as usize] {
            yb -= 1;
        }
    }
    let mut blow = -x0;
    let mut bhgh = -x0;
    vb[(blow + off) as usize] = yb;
    vb[(blow - 1 + off) as usize] = n + 1;

    let mut flow = 0i64;
    for d in 1i64.. {
        // Forward wave.
        flow -= 1;
        vf[(flow - 1 + off) as usize] = -2;
        let mut am = -2i64;
        let mut ac = -2i64;
        for k in (flow..=d).rev() {
            let ap = ac;
            ac = am + 1;
            am = vf[(k - 1 + off) as usize];
            let cand = if ac < am {
                if ap < am {
                    am
                } else {
                    ap
                }
            } else if ap < ac {
                ac
            } else {
                ap
            };
            if (blow..=bhgh).contains(&k) {
                let r = vb[(k + off) as usize];
                if cand > r {
                    let yy = if ap > r {
                        ap
                    } else if ac > r {
                        ac
                    } else {
                        r + 1
                    };
                    return ((2 * d - 1) as usize, (k + yy) as usize, yy as usize);
                }
            }
            let mut yv = cand;
            let x = m - k;
            // Snake: extend y while both sequences stay in bounds and match.
            while yv < n && yv < x && t[yv as usize] == q[(k + yv) as usize] {
                yv += 1;
            }
            vf[(k + off) as usize] = yv;
        }

        // Backward wave.
        bhgh += 1;
        blow -= 1;
        vb[(blow - 1 + off) as usize] = n + 1;
        let mut am = n + 1;
        let mut ac = n + 1;
        for k in (blow..=bhgh).rev() {
            let ap = ac + 1;
            ac = am;
            am = vb[(k - 1 + off) as usize];
            let cand = if ac > am {
                if ap > am {
                    am
                } else {
                    ap
                }
            } else if ap > ac {
                ac
            } else {
                ap
            };
            if (flow..=d).contains(&k) {
                let r = vf[(k + off) as usize];
                if cand <= r {
                    let yy = if ap <= r {
                        ap
                    } else if ac <= r {
                        ac
                    } else {
                        r
                    };
                    return ((2 * d) as usize, (k + yy) as usize, yy as usize);
                }
            }
            let mut yv = cand - 1;
            let x = -k;
            while yv >= 0
                && yv >= x
                && yv < n
                && (k + yv) < m
                && t[yv as usize] == q[(k + yv) as usize]
            {
                yv -= 1;
            }
            vb[(k + off) as usize] = yv;
        }
    }
    unreachable!("split_nd must meet within the wave loop");
}

/// Myers `dandc_nd`: append the exact edit script of `q` vs `t` to `ops`.
///
/// Coordinates are absolute (`q_abs`/`t_abs` are the slice offsets in the
/// full sequences). Returns the edit distance.
fn dandc_nd(q: &[u8], t: &[u8], q_abs: usize, t_abs: usize, ops: &mut Vec<EditOp>) -> usize {
    let m = q.len();
    let n = t.len();
    if m == 0 {
        for j in 0..n {
            ops.push(EditOp::Ins {
                q_pos: q_abs.saturating_sub(1),
                t_pos: t_abs + j,
            });
        }
        return n;
    }
    if n == 0 {
        for i in 0..m {
            ops.push(EditOp::Del {
                q_pos: q_abs + i,
                t_pos: t_abs,
            });
        }
        return m;
    }
    let (d, x, y) = split_nd(q, t);
    if d > 1 {
        let e1 = dandc_nd(&q[..x], &t[..y], q_abs, t_abs, ops);
        let e2 = dandc_nd(&q[x..], &t[y..], q_abs + x, t_abs + y, ops);
        e1 + e2
    } else if d == 1 {
        if m > n {
            ops.push(EditOp::Del {
                q_pos: q_abs + x - 1,
                t_pos: t_abs + y,
            });
        } else if m < n {
            ops.push(EditOp::Ins {
                q_pos: q_abs.saturating_sub(1),
                t_pos: t_abs + y - 1,
            });
        }
        // m == n: a substitution, implied by the surrounding match run.
        1
    } else {
        0
    }
}

/// Banded unit-cost edit script for a diagonal-restricted path (FastGA's
/// in-box DP): every aligned column keeps `t_pos - q_pos` inside
/// `[k_lo, k_hi]` (absolute coordinates), with `k_lo == k_hi` degenerating
/// to a single diagonal. Used by self mode, where the wave anchors are
/// clipped to one side of diagonal 0 and the exact D&C path must not cross
/// it either. Returns the edit distance.
fn banded_edit_ops(
    q: &[u8],
    t: &[u8],
    q_abs: usize,
    t_abs: usize,
    k_lo: i64,
    k_hi: i64,
    ops: &mut Vec<EditOp>,
) -> usize {
    let m = q.len() as i64;
    let n = t.len() as i64;
    let d0 = t_abs as i64 - q_abs as i64;
    let width = (k_hi - k_lo + 1) as usize;
    let inf = u32::MAX / 4;
    let mut dp = vec![inf; (m as usize + 1) * width];
    let mut tr = vec![0u8; (m as usize + 1) * width];
    dp[(d0 - k_lo) as usize] = 0; // (0, 0): the path starts on the anchor diagonal
    for i in 0..=m {
        let j_lo = (i + k_lo - d0).max(0);
        let j_hi = (i + k_hi - d0).min(n);
        for j in j_lo..=j_hi {
            if i == 0 && j == 0 {
                continue;
            }
            let k = d0 + j - i;
            let off = (k - k_lo) as usize;
            let c = (i as usize) * width + off;
            let mut best = inf;
            let mut best_op = 0u8;
            if i > 0 {
                let kp = d0 + j - (i - 1);
                if (k_lo..=k_hi).contains(&kp) {
                    let pv = dp[((i - 1) as usize) * width + (kp - k_lo) as usize];
                    if pv + 1 < best {
                        best = pv + 1;
                        best_op = 1; // delete q[i-1]
                    }
                }
            }
            if j > 0 {
                let kp = d0 + (j - 1) - i;
                if (k_lo..=k_hi).contains(&kp) {
                    let pv = dp[(i as usize) * width + (kp - k_lo) as usize];
                    if pv + 1 < best {
                        best = pv + 1;
                        best_op = 2; // insert t[j-1]
                    }
                }
            }
            if i > 0 && j > 0 {
                let pv = dp[((i - 1) as usize) * width + off];
                let sub = u32::from(q[(i - 1) as usize] != t[(j - 1) as usize]);
                if pv + sub < best {
                    best = pv + sub;
                    best_op = 3; // match / substitution
                }
            }
            dp[c] = best;
            tr[c] = best_op;
        }
    }
    let d_end = d0 + n - m;
    debug_assert!((k_lo..=k_hi).contains(&d_end), "end diagonal outside band");
    let total = dp[(m as usize) * width + (d_end - k_lo) as usize];
    debug_assert!(total < inf, "no banded path found");
    let mut i = m;
    let mut j = n;
    let mut stack = Vec::new();
    while i > 0 || j > 0 {
        let k = d0 + j - i;
        let off = (k - k_lo) as usize;
        match tr[(i as usize) * width + off] {
            1 => {
                stack.push(EditOp::Del {
                    q_pos: q_abs + i as usize - 1,
                    t_pos: t_abs + j as usize,
                });
                i -= 1;
            }
            2 => {
                stack.push(EditOp::Ins {
                    q_pos: q_abs.saturating_sub(1),
                    t_pos: t_abs + j as usize - 1,
                });
                j -= 1;
            }
            _ => {
                i -= 1;
                j -= 1;
            }
        }
    }
    ops.extend(stack.into_iter().rev());
    total as usize
}

/// Expand an edit script into aligned columns; returns `(q_aln, t_aln, matches)`.
fn ops_to_columns(
    q: &[u8],
    t: &[u8],
    q_start: usize,
    t_start: usize,
    q_end: usize,
    t_end: usize,
    ops: &[EditOp],
) -> (Vec<u8>, Vec<u8>, usize) {
    let mut q_aln = Vec::with_capacity(q.len() + ops.len());
    let mut t_aln = Vec::with_capacity(t.len() + ops.len());
    let mut matches = 0usize;
    let mut pa = 0usize;
    let mut pb = 0usize;
    for op in ops {
        let run = match *op {
            EditOp::Del { q_pos, .. } => q_pos - q_start - pa,
            EditOp::Ins { t_pos, .. } => t_pos - t_start - pb,
        };
        for _ in 0..run {
            q_aln.push(q[q_start + pa]);
            t_aln.push(t[t_start + pb]);
            if q[q_start + pa] == t[t_start + pb] {
                matches += 1;
            }
            pa += 1;
            pb += 1;
        }
        match *op {
            EditOp::Del { .. } => {
                q_aln.push(q[q_start + pa]);
                t_aln.push(b'-');
                pa += 1;
            }
            EditOp::Ins { .. } => {
                q_aln.push(b'-');
                t_aln.push(t[t_start + pb]);
                pb += 1;
            }
        }
    }
    let tail = (q_end - q_start - pa).min(t_end - t_start - pb);
    for _ in 0..tail {
        q_aln.push(q[q_start + pa]);
        t_aln.push(t[t_start + pb]);
        if q[q_start + pa] == t[t_start + pb] {
            matches += 1;
        }
        pa += 1;
        pb += 1;
    }
    debug_assert_eq!(pa, q_end - q_start);
    debug_assert_eq!(pb, t_end - t_start);
    (q_aln, t_aln, matches)
}

/// One wavefront cell: furthest anti (`v`), the last `PATH_LEN` columns as a
/// match bitvector, and the match count in that window.
#[derive(Clone, Copy)]
struct WaveCell {
    v: i64,
    bits: u64,
    m: u32,
}

const PATH_LEN: u32 = 60;
const PATH_TOP: u64 = 1 << PATH_LEN;
const PATH_INT: u64 = PATH_TOP - 1;
const PATH_AVE: u32 = 42; // PATH_LEN * (1 - (1 - 0.7)) for unbiased bases
const WAVE_LAG: i64 = 70;
const D_CAP: usize = 5_000_000;

/// FastGA `forward_wave` from a mid-line: for every diagonal of the band
/// `[k_lo, k_hi]`, the 0-wave starts a match snake at `(mida + k) / 2` and the
/// wavefront expands one diagonal per edit. `minp`/`maxp` hard-clip the
/// diagonal range during the expansion (FastGA's self-mode boundaries).
///
/// Only the tip is kept (no per-wave history): the endpoint is the last wave
/// maximum whose `PATH_LEN`-column window has at least `PATH_AVE` matches
/// (FastGA's trim point). Returns the trim point as `(a, b)` coordinates
/// (exclusive end), or `None` when nothing extends.
fn forward_wave_mid(
    a: &[u8],
    b: &[u8],
    k_lo: i64,
    k_hi: i64,
    mida: i64,
    minp: Option<i64>,
    maxp: Option<i64>,
) -> Option<(i64, i64)> {
    let n = a.len() as i64;
    let m = b.len() as i64;
    let mut low = k_lo.max(-m);
    let mut high = k_hi.min(n);
    if let Some(p) = minp {
        low = low.max(p);
    }
    if let Some(p) = maxp {
        high = high.min(p);
    }
    if low > high {
        return None;
    }
    let dead = WaveCell {
        v: -1,
        bits: 0,
        m: 0,
    };

    // 0-wave: each diagonal starts its snake at the mid-line.
    let mut prev: Vec<WaveCell> = Vec::with_capacity((high - low + 1) as usize);
    let mut besta = mida;
    let mut last_good = mida;
    let mut trim = (mida, (mida + high) >> 1);
    let mut found = false;
    for k in (low..=high).rev() {
        let mut x = (mida + k) >> 1;
        let mut y = x - k;
        if x < 0 || x > n || y < 0 || y > m {
            prev.push(dead);
            continue;
        }
        let mut bits = PATH_INT;
        let mut mcnt = PATH_LEN;
        while x < n && y < m && a[x as usize] == b[y as usize] {
            if (bits & PATH_TOP) == 0 {
                mcnt += 1;
            }
            bits = (bits << 1) | 1;
            x += 1;
            y += 1;
        }
        let c = (x << 1) - k;
        prev.push(WaveCell {
            v: c,
            bits,
            m: mcnt,
        });
        if c > besta {
            besta = c;
            found = true;
            if mcnt >= PATH_AVE {
                last_good = c;
                trim = (c, x);
            }
        }
    }
    if !found {
        return None;
    }

    let mut d = 0usize;
    let mut prev_max = besta;
    let mut last_good_d = 0usize;
    let mut cur: Vec<WaveCell> = Vec::new();
    loop {
        // Expand the band by one diagonal on each side, clipped by the
        // hard boundaries and the sequence length.
        low = (low - 1).max(-m);
        high = (high + 1).min(n);
        if let Some(p) = minp {
            low = low.max(p);
        }
        if let Some(p) = maxp {
            high = high.min(p);
        }
        if low > high {
            break;
        }
        d += 1;
        if d >= D_CAP {
            break;
        }
        let width = (high - low + 1) as usize;
        cur.resize(width, dead);
        cur[..width].fill(dead);
        let poff = |k: i64| (k - (low + 1)) as usize;
        let coff = |k: i64| (k - low) as usize;
        let mut new_max = -1i64;
        for k in (low..=high).rev() {
            let am = if k > low + 1 { prev[poff(k - 1)] } else { dead };
            let ac = if (low + 1..=high - 1).contains(&k) {
                prev[poff(k)]
            } else {
                dead
            };
            let ap = if k < high - 1 {
                prev[poff(k + 1)]
            } else {
                dead
            };
            // Three-branch update (FastGA order); the chosen predecessor also
            // provides the tip bitvector and match count.
            let (cand, mut bits, mut mcnt) = if ac.v < am.v {
                if am.v < ap.v {
                    (ap.v + 1, ap.bits, ap.m)
                } else {
                    (am.v + 1, am.bits, am.m)
                }
            } else if ac.v < ap.v {
                (ap.v + 1, ap.bits, ap.m)
            } else {
                (ac.v + 2, ac.bits, ac.m)
            };
            if cand < 0 {
                continue;
            }
            // The edit step shifts a 0 into the tip.
            if (bits & PATH_TOP) != 0 {
                mcnt -= 1;
            }
            bits <<= 1;
            let mut x = (cand + k) >> 1;
            let mut y = x - k;
            if x > n {
                x = n;
                y = x - k;
            }
            if y > m {
                y = m;
                x = y + k;
            }
            if x < 0 || x > n || y < 0 || y > m {
                continue;
            }
            while x < n && y < m && a[x as usize] == b[y as usize] {
                if (bits & PATH_TOP) == 0 {
                    mcnt += 1;
                }
                bits = (bits << 1) | 1;
                x += 1;
                y += 1;
            }
            let c = (x << 1) - k;
            cur[coff(k)] = WaveCell {
                v: c,
                bits,
                m: mcnt,
            };
            if c > new_max {
                new_max = c;
            }
            if c > besta {
                besta = c;
                if mcnt >= PATH_AVE {
                    last_good = c;
                    trim = (c, x);
                    last_good_d = d;
                }
            }
        }
        if new_max <= prev_max {
            break; // every active diagonal hit a sequence end
        }
        prev_max = new_max;
        // WAVE_LAG pruning: keep only cells near the best point.
        let gate = besta - WAVE_LAG;
        let mut alow = low;
        let mut ahigh = high;
        while alow <= ahigh && cur[coff(alow)].v < gate {
            alow += 1;
        }
        while ahigh >= alow && cur[coff(ahigh)].v < gate {
            ahigh -= 1;
        }
        if alow > ahigh {
            break;
        }
        prev.clear();
        prev.extend_from_slice(&cur[coff(alow)..=coff(ahigh)]);
        low = alow;
        high = ahigh;
        if besta - last_good >= TRIM_MLAG {
            break; // the tip has outrun the last good point
        }
        // FastGA lets a frozen trim crawl TRIM_MLAG (250 anti) past the last
        // good point before stopping; calls whose trim has been frozen for
        // this long without recovering return the same endpoint, so stop the
        // crawl early and save most of the wasted waves.
        if d - last_good_d > 60 {
            break;
        }
    }

    let (trim_anti, trim_x) = trim;
    Some((trim_x, trim_anti - trim_x))
}

/// Result of a FastGA-style `Local_Alignment`: the aligned span of the full
/// contig sequences plus its exact edit script.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalAlign {
    /// Query bases with `-` for gaps.
    pub q_aln: Vec<u8>,
    /// Target bases with `-` for gaps (same length as `q_aln`).
    pub t_aln: Vec<u8>,
    /// First aligned query base (orientation space).
    pub q_start: usize,
    /// First aligned target base.
    pub t_start: usize,
    /// End of the aligned query span (exclusive).
    pub q_end: usize,
    /// End of the aligned target span (exclusive).
    pub t_end: usize,
    /// Number of equal columns.
    pub matches: usize,
    /// Edit distance of the span (substitution = 1).
    pub diffs: usize,
}

/// FastGA `Local_Alignment`: bidirectional mid-line waves over the full
/// sequences within the diagonal band, plus an exact Myers D&C edit script.
///
/// `q` is the query (orientation space), `t` the target; `rt`/`rq` are the
/// reversed sequences (reused across calls of one tube); `amid` is the
/// mid-line anti-diagonal and `[dgmin, dgmax]` the tube's diagonal band.
/// `selfie` applies FastGA's self-mode diagonal boundaries: a same-contig
/// forward self-alignment must not cross diagonal 0 (the exact self-identity
/// line), so tubes entirely on one side are clipped there and tubes
/// straddling 0 are skipped.
#[allow(clippy::too_many_arguments)]
pub fn local_alignment(
    q: &[u8],
    t: &[u8],
    rt: &[u8],
    rq: &[u8],
    dgmin: i64,
    dgmax: i64,
    amid: i64,
    selfie: bool,
) -> Option<LocalAlign> {
    const DUB_TRIM: i64 = 45;
    let n = q.len() as i64;
    let m = t.len() as i64;
    if n == 0 || m == 0 {
        return None;
    }
    let (minp, maxp) = if selfie {
        if dgmin > 0 {
            (Some(1), None)
        } else if dgmax < 0 {
            (None, Some(-1))
        } else {
            return None;
        }
    } else {
        (None, None)
    };
    // The reverse wave runs on mirrored sequences (k' = m - n - k), so its
    // hard boundaries are the mirror of the forward ones.
    let (r_minp, r_maxp) = (maxp.map(|p| (m - n) - p), minp.map(|p| (m - n) - p));
    // Forward wave from the mid-line up (a = target, b = query).
    let (at, bt) = forward_wave_mid(t, q, dgmin, dgmax, amid, minp, maxp)?;
    // Reverse wave on mirrored sequences (extends downward in original space).
    let mida_rev = n + m - 2 - amid;
    let k_lo = (m - n) - dgmax;
    let k_hi = (m - n) - dgmin;
    let (ar, br) = forward_wave_mid(rt, rq, k_lo, k_hi, mida_rev, r_minp, r_maxp)?;
    // The mirrored trim point is an exclusive end; the original path start
    // (inclusive) is the mirror of `end - 1`.
    let (ab, bb) = (m - ar, n - br);
    let fshort = at + bt - amid < DUB_TRIM;
    let rshort = amid - (ab + bb) < DUB_TRIM;
    let (at, bt, ab, bb) = match (fshort, rshort) {
        (true, true) => return None,
        (true, false) => {
            let (at, bt) = forward_wave_mid(t, q, ab - bb, ab - bb, ab + bb, minp, maxp)?;
            (at, bt, ab, bb)
        }
        (false, true) => {
            let mida_rev = n + m - 2 - (at + bt);
            let k = (m - n) - (at - bt);
            let (ar, br) = forward_wave_mid(rt, rq, k, k, mida_rev, r_minp, r_maxp)?;
            (at, bt, m - ar, n - br)
        }
        (false, false) => (at, bt, ab, bb),
    };
    if at <= ab || bt <= bb {
        return None;
    }
    let q_span = &q[bb as usize..bt as usize];
    let t_span = &t[ab as usize..at as usize];
    let mut ops = Vec::new();
    let diffs = if selfie {
        // FastGA's in-box DP keeps the path near the wave (its band is the
        // trace-point diagonal spread); pgr has only the anchors, whose
        // diagonals can drift from the tube band at copy boundaries, so the
        // band is the union of the tube band and the anchor diagonals,
        // clipped to the non-zero side in self mode.
        let dg0 = ab - bb;
        let dg1 = at - bt;
        let lo = dg0.min(dg1).min(dgmin).max(minp.unwrap_or(i64::MIN));
        let hi = dg0.max(dg1).max(dgmax).min(maxp.unwrap_or(i64::MAX));
        banded_edit_ops(q_span, t_span, bb as usize, ab as usize, lo, hi, &mut ops)
    } else {
        dandc_nd(q_span, t_span, bb as usize, ab as usize, &mut ops)
    };
    let (q_aln, t_aln, matches) = ops_to_columns(
        q,
        t,
        bb as usize,
        ab as usize,
        bt as usize,
        at as usize,
        &ops,
    );
    Some(LocalAlign {
        q_aln,
        t_aln,
        q_start: bb as usize,
        t_start: ab as usize,
        q_end: bt as usize,
        t_end: at as usize,
        matches,
        diffs,
    })
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
    // Safety cap on stored waves (memory bound); realistic tubes never reach
    // it because the trim stops ~TRIM_MLAG waves after the best point.
    const D_CAP: usize = 500_000;

    let mut history = WaveHistory {
        // Reserve ~256 waves of cells (the trim stops ~TRIM_MLAG waves after
        // the best point); a 4096-wave reservation wasted most of its memory
        // (24 B/cell) on every tube call.
        v: Vec::with_capacity((256 * width).min(D_CAP)),
        trace: Vec::with_capacity((256 * width).min(D_CAP)),
        width,
        k_lo,
        best_d: 0,
        best_k: k0,
        best_c: -1,
    };
    let off = |k: i64| (k - k_lo) as usize;

    // 0-wave: snake from the anchor point on its diagonal.
    history.v.extend(std::iter::repeat_n(-1, width));
    history
        .trace
        .extend(std::iter::repeat_n(Trace::default(), width));
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
    let mut prev_max = c;
    let mut d = 1usize;
    while besta - prev_max < TRIM_MLAG {
        let lo = (k0 - d as i64).max(k_lo);
        let hi = (k0 + d as i64).min(k_hi);
        if lo > hi {
            break;
        }
        history.v.extend(std::iter::repeat_n(-1, width));
        history
            .trace
            .extend(std::iter::repeat_n(Trace::default(), width));
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
        if new_besta <= prev_max {
            break; // no progress: every diagonal hit a sequence end
        }
        prev_max = new_besta;
        d += 1;
        if d >= D_CAP {
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

    /// Edit distance with substitution = 1 (FastGA `split_nd` metric).
    fn dp_edit(a: &[u8], b: &[u8]) -> usize {
        let (n, m) = (a.len(), b.len());
        let mut prev: Vec<usize> = (0..=m).collect();
        let mut cur = vec![0usize; m + 1];
        for i in 1..=n {
            cur[0] = i;
            for j in 1..=m {
                let sub = if a[i - 1] == b[j - 1] { 0 } else { 1 };
                cur[j] = (prev[j - 1] + sub).min(prev[j] + 1).min(cur[j - 1] + 1);
            }
            std::mem::swap(&mut prev, &mut cur);
        }
        prev[m]
    }

    fn rand_seq(len: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..len)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize & 3]
            })
            .collect()
    }

    fn all_seqs(len: usize) -> Vec<Vec<u8>> {
        if len == 0 {
            return vec![Vec::new()];
        }
        let mut out = Vec::new();
        for p in 0..(1usize << len) {
            out.push(
                (0..len)
                    .map(|i| if (p >> i) & 1 == 0 { b'A' } else { b'C' })
                    .collect(),
            );
        }
        out
    }

    #[test]
    fn dandc_exhaustive_small_cases() {
        for m in 0..=5usize {
            for n in 0..=5usize {
                for q in all_seqs(m) {
                    for t in all_seqs(n) {
                        let mut ops = Vec::new();
                        let d = dandc_nd(&q, &t, 0, 0, &mut ops);
                        assert_eq!(d, dp_edit(&q, &t), "D mismatch q={q:?} t={t:?} ops={ops:?}");
                        let res = std::panic::catch_unwind(|| {
                            ops_to_columns(&q, &t, 0, 0, q.len(), t.len(), &ops)
                        });
                        let (qa, ta, matches) = match res {
                            Ok(v) => v,
                            Err(_) => {
                                panic!("ops_to_columns panic q={q:?} t={t:?} ops={ops:?} d={d}")
                            }
                        };
                        assert_eq!(
                            qa.len(),
                            ta.len(),
                            "length mismatch q={q:?} t={t:?} ops={ops:?}"
                        );
                        assert_eq!(
                            qa.len(),
                            matches + d,
                            "edit count mismatch q={q:?} t={t:?} ops={ops:?}"
                        );
                        let qq: Vec<u8> = qa.iter().copied().filter(|&c| c != b'-').collect();
                        let tt: Vec<u8> = ta.iter().copied().filter(|&c| c != b'-').collect();
                        assert_eq!(qq, q);
                        assert_eq!(tt, t);
                    }
                }
            }
        }
    }

    #[test]
    fn forward_wave_mid_extends_identical_sequences_to_end() {
        let s = b"ACGTACGTACGT";
        // Mid-line at anti 11 (odd parity): the snake still reaches the end.
        let (x, y) = forward_wave_mid(s, s, -4, 4, 11, None, None).unwrap();
        assert_eq!((x, y), (12, 12), "identical sequences must align fully");
    }

    #[test]
    fn local_alignment_covers_conserved_region() {
        let cons = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 100 bp
        let mut q = b"TTTTTTT".to_vec();
        q.extend_from_slice(cons);
        q.extend_from_slice(b"GGGGGGG");
        let mut t = b"CCCCCCC".to_vec();
        t.extend_from_slice(cons);
        t.extend_from_slice(b"AAAAAAA");
        // Mid-line anti through the middle of the conserved block.
        let rt: Vec<u8> = t.iter().rev().copied().collect();
        let rq: Vec<u8> = q.iter().rev().copied().collect();
        let aln = local_alignment(&q, &t, &rt, &rq, -2, 2, 114, false).unwrap();
        assert!(
            aln.q_start <= 7 && aln.q_end >= 7 + cons.len(),
            "query span {:?}..{:?}",
            aln.q_start,
            aln.q_end
        );
        assert!(
            aln.t_start <= 7 && aln.t_end >= 7 + cons.len(),
            "target span {:?}..{:?}",
            aln.t_start,
            aln.t_end
        );
        assert_eq!(aln.q_aln.len(), aln.t_aln.len());
        let qq: Vec<u8> = aln.q_aln.iter().copied().filter(|&c| c != b'-').collect();
        let tt: Vec<u8> = aln.t_aln.iter().copied().filter(|&c| c != b'-').collect();
        assert_eq!(qq, &q[aln.q_start..aln.q_end]);
        assert_eq!(tt, &t[aln.t_start..aln.t_end]);
    }

    #[test]
    fn local_alignment_self_clips_diagonal_zero() {
        // q/t are the same length; the conserved block sits on a positive
        // diagonal (+3) in the first case and a negative one (-3) in the
        // second. Self mode must never let the path cross diagonal 0.
        let cons = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 100 bp
        let pad = |n: usize, seed: u64| {
            let mut x = seed;
            (0..n)
                .map(|_| {
                    x = x
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    b"ACGT"[(x >> 33) as usize & 3]
                })
                .collect::<Vec<u8>>()
        };

        // Positive diagonal: q's block is 3 bp left of t's.
        let t = {
            let mut v = pad(20, 1);
            v.extend_from_slice(cons);
            v.extend_from_slice(&pad(20, 2));
            v
        };
        let q = {
            let mut v = pad(17, 3);
            v.extend_from_slice(cons);
            v.extend_from_slice(&pad(23, 4));
            v
        };
        let rt: Vec<u8> = t.iter().rev().copied().collect();
        let rq: Vec<u8> = q.iter().rev().copied().collect();
        let aln = local_alignment(&q, &t, &rt, &rq, 3, 3, 137, true).unwrap();
        let mut t_i = aln.t_start as i64;
        let mut q_i = aln.q_start as i64;
        let mut min_dg = i64::MAX;
        for (qc, tc) in aln.q_aln.iter().zip(&aln.t_aln) {
            if *qc != b'-' && *tc != b'-' {
                min_dg = min_dg.min(t_i - q_i);
            }
            if *qc != b'-' {
                q_i += 1;
            }
            if *tc != b'-' {
                t_i += 1;
            }
        }
        assert!(
            min_dg >= 1,
            "positive-diagonal self path crossed 0: {min_dg}"
        );

        // Negative diagonal: q's block is 3 bp right of t's.
        let q = {
            let mut v = pad(23, 5);
            v.extend_from_slice(cons);
            v.extend_from_slice(&pad(17, 6));
            v
        };
        let rq: Vec<u8> = q.iter().rev().copied().collect();
        let aln = local_alignment(&q, &t, &rt, &rq, -3, -3, 137, true).unwrap();
        let mut t_i = aln.t_start as i64;
        let mut q_i = aln.q_start as i64;
        let mut max_dg = i64::MIN;
        for (qc, tc) in aln.q_aln.iter().zip(&aln.t_aln) {
            if *qc != b'-' && *tc != b'-' {
                max_dg = max_dg.max(t_i - q_i);
            }
            if *qc != b'-' {
                q_i += 1;
            }
            if *tc != b'-' {
                t_i += 1;
            }
        }
        assert!(
            max_dg <= -1,
            "negative-diagonal self path crossed 0: {max_dg}"
        );

        // A tube straddling diagonal 0 is skipped entirely.
        assert!(
            local_alignment(&q, &t, &rt, &rq, -2, 2, 137, true).is_none(),
            "straddling tube must be skipped in self mode"
        );
    }

    #[test]
    fn split_nd_matches_dp_edit_distance() {
        for seed in 0..40u64 {
            let q = rand_seq((seed % 12) as usize + 1, seed);
            let t = rand_seq(((seed * 7) % 12) as usize + 1, seed + 999);
            let (d, x, y) = split_nd(&q, &t);
            assert_eq!(d, dp_edit(&q, &t), "q={q:?} t={t:?}");
            assert!(x <= q.len() && y <= t.len(), "split out of range");
        }
        // Identical case.
        let s = b"ACGTACGT";
        assert_eq!(split_nd(s, s), (0, 8, 8));
    }

    #[test]
    fn split_nd_random_larger_cases() {
        for seed in 0..2000u64 {
            let m = (seed % 9) as usize + 3;
            let n = ((seed * 13) % 9) as usize + 3;
            let q = rand_seq(m, seed);
            let t = rand_seq(n, seed + 77);
            let (d, x, y) = split_nd(&q, &t);
            assert_eq!(d, dp_edit(&q, &t), "q={q:?} t={t:?}");
            assert!(x <= q.len() && y <= t.len(), "split out of range");
        }
    }

    #[test]
    fn split_nd_exhaustive_binary_6() {
        for m in 1..=6usize {
            for n in 1..=6usize {
                for q in all_seqs(m) {
                    for t in all_seqs(n) {
                        let (d, x, y) = split_nd(&q, &t);
                        assert_eq!(
                            d,
                            dp_edit(&q, &t),
                            "D mismatch q={q:?} t={t:?} split=({x},{y})"
                        );
                        assert!(
                            x <= q.len() && y <= t.len(),
                            "split out of range q={q:?} t={t:?} ({x},{y})"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn split_nd_random_catch_panics() {
        for seed in 0..20000u64 {
            let m = (seed % 8) as usize + 2;
            let n = 12 - m;
            let q = rand_seq(m, seed);
            let t = rand_seq(n, seed + 555);
            let res = std::panic::catch_unwind(|| split_nd(&q, &t));
            if res.is_err() {
                panic!("split_nd panic q={q:?} t={t:?}");
            }
        }
    }

    #[test]
    fn dandc_script_reconstructs_exact_alignment() {
        for seed in 0..40u64 {
            let q = rand_seq((seed % 12) as usize + 1, seed);
            let t = rand_seq(((seed * 7) % 12) as usize + 1, seed + 999);
            let mut ops = Vec::new();
            let d = dandc_nd(&q, &t, 0, 0, &mut ops);
            assert_eq!(d, dp_edit(&q, &t));
            let (qa, ta, matches) = ops_to_columns(&q, &t, 0, 0, q.len(), t.len(), &ops);
            assert_eq!(qa.len(), ta.len());
            assert_eq!(qa.len(), matches + d, "edits = len - matches");
            // De-gapped columns must equal the original sequences.
            let qq: Vec<u8> = qa.iter().copied().filter(|&c| c != b'-').collect();
            let tt: Vec<u8> = ta.iter().copied().filter(|&c| c != b'-').collect();
            assert_eq!(qq, q);
            assert_eq!(tt, t);
        }
        // Empty-adjacent cases (handled before split_nd).
        let s = b"ACGT";
        let mut ops = Vec::new();
        assert_eq!(dandc_nd(b"", s, 0, 0, &mut ops), 4);
        assert_eq!(ops.len(), 4);
        assert!(ops.iter().all(|o| matches!(o, EditOp::Ins { .. })));
        let mut ops = Vec::new();
        assert_eq!(dandc_nd(s, b"", 0, 0, &mut ops), 4);
        assert!(ops.iter().all(|o| matches!(o, EditOp::Del { .. })));
    }

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
