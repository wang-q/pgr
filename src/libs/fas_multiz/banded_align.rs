//! Banded DP alignment of two FasBlock reference sequences.
//!
//! [`banded_align_refs`] computes column-to-column alignments between the
//! reference entries of two [`FasBlock`]s using a banded dynamic programming
//! algorithm. The DP is a direct port of multiz's `mz_yama` engine: three
//! states (C substitution / D deletion / I insertion), quasi-natural gap
//! costs that depend on the last two edge types, and sum-of-pairs scoring
//! over all species columns of the two blocks. End-gaps (overhangs) are not
//! charged a gap-open penalty. The result is consumed by [`super::merge`] to
//! merge blocks into a unified alignment.

use super::{find_ref_entry, FasMultizConfig, FasMultizGapModel};
use crate::libs::chain::sub_matrix::SubMatrix;
use crate::libs::chain::GapCalc;
use crate::libs::fmt::fas::{FasBlock, FasEntry};
use std::collections::BTreeMap;

#[allow(clippy::type_complexity)]
pub(super) fn banded_align_refs(
    blocks: [&FasBlock; 2],
    ref_name: &str,
    cfg: &FasMultizConfig,
) -> anyhow::Result<Option<(Vec<Option<usize>>, Vec<Option<usize>>)>> {
    let submat = match &cfg.score_matrix {
        Some(name) => SubMatrix::from_name(name)?,
        None => SubMatrix::hoxd55(),
    };
    Ok(banded_align_refs_inner(blocks, ref_name, cfg, &submat))
}

#[allow(clippy::type_complexity)]
fn banded_align_refs_inner(
    blocks: [&FasBlock; 2],
    ref_name: &str,
    cfg: &FasMultizConfig,
    submat: &SubMatrix,
) -> Option<(Vec<Option<usize>>, Vec<Option<usize>>)> {
    let ref_a = find_ref_entry(blocks[0], ref_name)?;
    let ref_b = find_ref_entry(blocks[1], ref_name)?;

    let sa = ref_a.seq();
    let sb = ref_b.seq();

    let n = sa.len();
    let m = sb.len();

    if n == 0 || m == 0 {
        return None;
    }

    // Reference-anchored band constraints (multiz `pre_yama` LB/RB): pair the
    // k-th base of each reference (both sequences are the same locus), which
    // pins the DP grid row for ref_a column i to ref_b column j. `smooth`
    // then makes the bounds monotone and expands them into a radius-wide
    // sausage around the reference diagonal. Without this the free end-gaps
    // would let the DP dump the column-count difference at the block ends.
    let radius = cfg.radius;
    let mut lb = vec![0usize; n + 1];
    let mut rb = vec![m; n + 1];
    {
        let mut ia = 0usize;
        let mut ib = 0usize;
        while ia < n && ib < m {
            while ia < n && sa[ia] == b'-' {
                ia += 1;
            }
            while ib < m && sb[ib] == b'-' {
                ib += 1;
            }
            if ia >= n || ib >= m {
                break;
            }
            // Each ref_a base column is visited once, so a plain assignment
            // pins the grid row to the matching ref_b base column (multiz
            // `pre_yama` uses the same one-to-one walk with sentinel checks).
            lb[ia] = ib;
            rb[ia] = ib;
            ia += 1;
            ib += 1;
        }
    }

    // Smooth the bounds exactly as multiz `smooth`: monotone first, then the
    // radius-wide sausage.
    let radi = radius.min(n);
    let mut run = 0usize;
    for v in lb.iter_mut() {
        run = run.max(*v);
        *v = run;
    }
    run = m;
    for v in rb.iter_mut().rev() {
        run = run.min(*v);
        *v = run;
    }
    for i in (radi + 1..=n).rev() {
        lb[i] = lb[i].saturating_sub(radi).min(lb[i - radi]);
    }
    for v in lb[..=radi.min(n)].iter_mut() {
        *v = 0;
    }
    for i in 0..n.saturating_sub(radi) {
        rb[i] = (rb[i].saturating_add(radi)).min(m).max(rb[i + radi]);
    }
    for v in rb[n.saturating_sub(radi)..].iter_mut() {
        *v = m;
    }
    for i in 0..=n {
        if rb[i] < lb[i] {
            return None;
        }
    }
    // Per-row storage: row i occupies `row_len[i]` cells for j in lb[i]..=rb[i].
    let mut row_off = vec![0usize; n + 1];
    let mut total = 0usize;
    for i in 0..=n {
        row_off[i] = total;
        total += rb[i] + 1 - lb[i];
    }
    let cell = |i: usize, j: usize| -> Option<usize> {
        if j < lb[i] || j > rb[i] {
            None
        } else {
            Some(row_off[i] + (j - lb[i]))
        }
    };

    let (gap_open_pen, gap_extend_pen) =
        if let (Some(open), Some(extend)) = (cfg.gap_open, cfg.gap_extend) {
            let scale = cfg.match_score as f64 / 100.0;
            let open_scaled = (open as f64 * scale).round() as i32;
            let extend_scaled = (extend as f64 * scale).round() as i32;
            (-open_scaled, -extend_scaled)
        } else {
            match cfg.gap_model {
                FasMultizGapModel::Constant => (cfg.gap_score, cfg.gap_score),
                FasMultizGapModel::Medium | FasMultizGapModel::Loose => {
                    let gap_calc = match cfg.gap_model {
                        FasMultizGapModel::Medium => GapCalc::medium(),
                        FasMultizGapModel::Loose => GapCalc::loose(),
                        FasMultizGapModel::Constant => {
                            unreachable!("Constant gap model already handled in outer branch")
                        }
                    };
                    let c1 = gap_calc.calc(1, 0).max(1);
                    let c2 = gap_calc.calc(2, 0).max(c1 + 1);
                    let open_raw = 2 * c1 - c2;
                    let extend_raw = c2 - c1;
                    let scale = cfg.match_score as f64 / 100.0;
                    let open_scaled = (open_raw as f64 * scale).round() as i32;
                    let extend_scaled = (extend_raw as f64 * scale).round() as i32;
                    (-open_scaled, -extend_scaled)
                }
            }
        };

    // All-species column profiles. `col_a[i]` holds the bases of every species
    // in blocks[0] at column i (sorted by species name), `col_b[j]` for
    // blocks[1]. A substitution edge scores sum-of-pairs over all K x L
    // species pairs, matching multiz `SS` (base-base matrix score; base
    // against a gap costs one gap-extension; gap-gap is 0).
    let mut seqs_a: Vec<(&str, &[u8])> = Vec::new();
    let mut map_a: BTreeMap<&str, &FasEntry> = BTreeMap::new();
    for (entry, name) in blocks[0].entries.iter().zip(blocks[0].names.iter()) {
        map_a.insert(name.as_str(), entry);
    }
    let mut names_a: Vec<&str> = map_a.keys().copied().collect();
    names_a.sort_unstable();
    for name in names_a {
        seqs_a.push((name, map_a[name].seq()));
    }

    let mut seqs_b: Vec<(&str, &[u8])> = Vec::new();
    let mut map_b: BTreeMap<&str, &FasEntry> = BTreeMap::new();
    for (entry, name) in blocks[1].entries.iter().zip(blocks[1].names.iter()) {
        map_b.insert(name.as_str(), entry);
    }
    let mut names_b: Vec<&str> = map_b.keys().copied().collect();
    names_b.sort_unstable();
    for name in names_b {
        seqs_b.push((name, map_b[name].seq()));
    }

    let col_a: Vec<Vec<u8>> = (0..n)
        .map(|i| {
            seqs_a
                .iter()
                .map(|(_, s)| s.get(i).copied().unwrap_or(b'-'))
                .collect()
        })
        .collect();
    let col_b: Vec<Vec<u8>> = (0..m)
        .map(|j| {
            seqs_b
                .iter()
                .map(|(_, s)| s.get(j).copied().unwrap_or(b'-'))
                .collect()
        })
        .collect();

    let k = col_a[0].len();
    let l = col_b[0].len();
    let na_col: Vec<usize> = col_a
        .iter()
        .map(|c| c.iter().filter(|b| **b != b'-').count())
        .collect();
    let nb_col: Vec<usize> = col_b
        .iter()
        .map(|c| c.iter().filter(|b| **b != b'-').count())
        .collect();

    // Substitution score for one column pair (multiz `SS`): base-base uses
    // the matrix; a base against a gap costs one gap-extension; gap-gap is 0.
    let ss = |ba: u8, bb: u8| -> i32 {
        if ba == b'-' && bb == b'-' {
            0
        } else if ba == b'-' || bb == b'-' {
            gap_extend_pen
        } else {
            submat.get_score(ba as char, bb as char) / 50
        }
    };
    // Quasi-natural gap-open lookup (multiz `GAP`): of the 16 configurations
    // of the last two edge types, six charge one gap-open penalty.
    let gop = |s: bool, t: bool, u: bool, v: bool| -> i32 {
        match (s, t, u, v) {
            (false, false, false, true)
            | (false, false, true, false)
            | (false, true, true, false)
            | (true, false, false, true)
            | (true, true, false, true)
            | (true, true, true, false) => gap_open_pen,
            _ => 0,
        }
    };

    // C (substitution), D (deletion: A column, B all gaps) and I (insertion:
    // B column, A all gaps) states, plus a trace byte encoding the previous
    // node type per state (flag_c | flag_d<<2 | flag_i<<4, values 0=C,1=I,2=D).
    let mut c = vec![i32::MIN; total];
    let mut d = vec![i32::MIN; total];
    let mut ins = vec![i32::MIN; total];
    let mut trace = vec![0u8; total];

    let k0 = cell(0, 0)?;
    c[k0] = 0;
    d[k0] = 0;
    ins[k0] = 0;

    // First row: only insertions (B columns against a fully dashed A).
    for j in 1..=rb[0] {
        let kj = cell(0, j)?;
        let prev = cell(0, j - 1)?;
        let ext = (nb_col[j - 1] as i32) * (k as i32) * gap_extend_pen;
        ins[kj] = ins[prev].saturating_add(ext);
        // The row-0 insertion chain always continues from the I state, so the
        // traceback keeps walking left (flag_i = 1, the I code).
        trace[kj] = 1 << 4;
    }

    for i in 1..=n {
        for j in lb[i]..=rb[i] {
            let k = cell(i, j)?;

            // ---- I state: B column j-1 inserted, A all dashes ----
            let mut i_best = i32::MIN;
            let mut i_flag = 0u8;
            if j > 0 {
                if let Some(pk) = cell(i, j - 1) {
                    let mut x = c[pk];
                    let mut y = d[pk];
                    let mut z = ins[pk];
                    // No gap-open penalty on the last row (end-gap).
                    if i < n {
                        let ca = &col_a[i - 1];
                        let cb = &col_b[j - 1];
                        let pb: &[u8] = if j > 1 { &col_b[j - 2] } else { &[] };
                        for &ba in ca.iter() {
                            let s = ba == b'-';
                            for (ib, &bb) in cb.iter().enumerate() {
                                let v = bb == b'-';
                                let t = j > 1 && pb[ib] == b'-';
                                if j > lb[i - 1] + 1 {
                                    x = x.saturating_add(gop(s, t, true, v));
                                }
                                y = y.saturating_add(gop(s, true, true, v));
                                if j > lb[i] + 1 {
                                    z = z.saturating_add(gop(true, t, true, v));
                                }
                            }
                        }
                    }
                    let ext = (nb_col[j - 1] as i32) * (k as i32) * gap_extend_pen;
                    if x >= y && x >= z {
                        i_best = x;
                        i_flag = 0;
                    } else if y > z {
                        i_best = y;
                        i_flag = 2;
                    } else {
                        i_best = z;
                        i_flag = 1;
                    }
                    i_best = i_best.saturating_add(ext);
                }
            }

            // ---- C state: A column i-1 aligned to B column j-1 ----
            let mut c_best = i32::MIN;
            let mut c_flag = 0u8;
            if i > 0 && j > 0 {
                if let Some(pk) = cell(i - 1, j - 1) {
                    let mut x = c[pk];
                    let mut y = d[pk];
                    let mut z = ins[pk];
                    // No gap-open penalty at the start (first column).
                    if j > 1 {
                        let ca = &col_a[i - 1];
                        let pa: &[u8] = if i > 1 { &col_a[i - 2] } else { &[] };
                        let cb = &col_b[j - 1];
                        let pb = &col_b[j - 2];
                        for (ia, &ba) in ca.iter().enumerate() {
                            let u = ba == b'-';
                            let s = i > 1 && pa[ia] == b'-';
                            for (ib, &bb) in cb.iter().enumerate() {
                                let v = bb == b'-';
                                let t = pb[ib] == b'-';
                                if i > 1 && j > lb[i - 2] + 1 {
                                    x = x.saturating_add(gop(s, t, u, v));
                                }
                                if i > 1 {
                                    y = y.saturating_add(gop(s, true, u, v));
                                }
                                if j > lb[i - 1] + 1 {
                                    z = z.saturating_add(gop(true, t, u, v));
                                }
                            }
                        }
                    }
                    let mut s = 0;
                    for &ba in &col_a[i - 1] {
                        for &bb in &col_b[j - 1] {
                            s += ss(ba, bb);
                        }
                    }
                    if x >= y && x >= z {
                        c_best = x;
                        c_flag = 0;
                    } else if y > z {
                        c_best = y;
                        c_flag = 2;
                    } else {
                        c_best = z;
                        c_flag = 1;
                    }
                    c_best = c_best.saturating_add(s);
                }
            }

            // ---- D state: A column i-1 inserted, B all dashes ----
            let mut d_best = i32::MIN;
            let mut d_flag = 0u8;
            if i > 0 {
                if let Some(pk) = cell(i - 1, j) {
                    let mut x = c[pk];
                    let mut y = d[pk];
                    let mut z = ins[pk];
                    // No gap-open penalty at the start or end (0<col<N).
                    if j > 0 && j < m {
                        let ca = &col_a[i - 1];
                        let pa: &[u8] = if i > 1 { &col_a[i - 2] } else { &[] };
                        let cb = &col_b[j - 1];
                        for (ia, &ba) in ca.iter().enumerate() {
                            let u = ba == b'-';
                            let s = i > 1 && pa[ia] == b'-';
                            for &bb in cb.iter() {
                                let t = bb == b'-';
                                if i > 1 && j > lb[i - 2] {
                                    x = x.saturating_add(gop(s, t, u, true));
                                }
                                if i > 1 {
                                    y = y.saturating_add(gop(s, true, u, true));
                                }
                                if j > lb[i - 1] {
                                    z = z.saturating_add(gop(true, t, u, true));
                                }
                            }
                        }
                    }
                    let ext = (na_col[i - 1] as i32) * (l as i32) * gap_extend_pen;
                    if x >= y && x >= z {
                        d_best = x;
                        d_flag = 0;
                    } else if y > z {
                        d_best = y;
                        d_flag = 2;
                    } else {
                        d_best = z;
                        d_flag = 1;
                    }
                    d_best = d_best.saturating_add(ext);
                }
            }

            c[k] = c_best;
            d[k] = d_best;
            ins[k] = i_best;
            trace[k] = c_flag | (d_flag << 2) | (i_flag << 4);
        }
    }

    let mut i = n;
    let mut j = m;
    let k = cell(i, j)?;
    let mut node = if c[k] >= d[k] && c[k] >= ins[k] {
        0
    } else if d[k] >= ins[k] {
        2
    } else {
        1
    };

    let mut map_a = Vec::new();
    let mut map_b = Vec::new();

    while i > 0 || j > 0 {
        let k = cell(i, j)?;
        let st = trace[k];
        match node {
            1 => {
                map_a.push(None);
                map_b.push(Some(j - 1));
                j -= 1;
                node = st >> 4;
            }
            2 => {
                map_a.push(Some(i - 1));
                map_b.push(None);
                i -= 1;
                node = (st >> 2) & 3;
            }
            _ => {
                map_a.push(Some(i - 1));
                map_b.push(Some(j - 1));
                i -= 1;
                j -= 1;
                node = st & 3;
            }
        }
    }
    map_a.reverse();
    map_b.reverse();

    if map_a.len() != map_b.len() || map_a.is_empty() {
        return None;
    }

    // Trim only leading/trailing columns that are all gaps across every
    // species of both blocks. Single-sided edge columns (overhangs) are
    // legitimate under the free end-gap model and keep their content.
    let all_gap = |ma: Option<usize>, mb: Option<usize>| -> bool {
        let a = match ma {
            Some(i) => col_a[i].iter().all(|b| *b == b'-'),
            None => true,
        };
        let b = match mb {
            Some(j) => col_b[j].iter().all(|b| *b == b'-'),
            None => true,
        };
        a && b
    };
    let mut left = 0usize;
    let mut right = map_a.len();
    while left < right && all_gap(map_a[left], map_b[left]) {
        left += 1;
    }
    while right > left && all_gap(map_a[right - 1], map_b[right - 1]) {
        right -= 1;
    }

    if left >= right {
        return None;
    }

    let map_a = map_a[left..right].to_vec();
    let map_b = map_b[left..right].to_vec();

    Some((map_a, map_b))
}
