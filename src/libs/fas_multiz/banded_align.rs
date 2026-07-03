//! Banded DP alignment of two FasBlock reference sequences.
//!
//! [`banded_align_refs`] computes column-to-column alignments between the
//! reference entries of two [`FasBlock`]s using a banded dynamic programming
//! algorithm with affine gaps. The result is consumed by [`super::merge`] to
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
    use std::cmp::min;

    let ref_a = find_ref_entry(blocks[0], ref_name)?;
    let ref_b = find_ref_entry(blocks[1], ref_name)?;

    let sa = ref_a.seq();
    let sb = ref_b.seq();

    let n = sa.len();
    let m = sb.len();

    if n == 0 || m == 0 {
        return None;
    }

    let band = cfg
        .radius
        .max(((n as isize - m as isize).unsigned_abs()) as usize);

    let width = 2 * band + 1;
    let mut score = vec![i32::MIN; (n + 1) * width];
    let mut gap_i = vec![i32::MIN; (n + 1) * width];
    let mut gap_j = vec![i32::MIN; (n + 1) * width];
    let mut trace = vec![0i8; (n + 1) * width];

    let idx = |i: usize, j: usize| -> Option<usize> {
        let band_start = i.saturating_sub(band);
        let band_end = min(m, i + band);
        if j < band_start || j > band_end {
            None
        } else {
            let offset = j + band - i;
            Some(i * width + offset)
        }
    };

    if let Some(k) = idx(0, 0) {
        score[k] = 0;
        gap_i[k] = i32::MIN;
        gap_j[k] = i32::MIN;
        trace[k] = 0;
    } else {
        return None;
    }

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
                        FasMultizGapModel::Constant => return None, // unreachable in this branch
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

    let mut profiles: Vec<(&[u8], &[u8])> = Vec::new();
    let mut map_a: BTreeMap<&str, &FasEntry> = BTreeMap::new();
    for (entry, name) in blocks[0].entries.iter().zip(blocks[0].names.iter()) {
        map_a.insert(name.as_str(), entry);
    }
    for (entry, name) in blocks[1].entries.iter().zip(blocks[1].names.iter()) {
        if let Some(ea) = map_a.get(name.as_str()) {
            profiles.push((ea.seq(), entry.seq()));
        }
    }

    for i in 0..=n {
        let band_start = i.saturating_sub(band);
        let band_end = min(m, i + band);
        for j in band_start..=band_end {
            let k = match idx(i, j) {
                Some(v) => v,
                None => continue,
            };
            if i == 0 && j == 0 {
                continue;
            }

            let mut best = i32::MIN;
            let mut bt = 0i8;

            if i > 0 && j > 0 {
                if let Some(pk) = idx(i - 1, j - 1) {
                    let mut s = 0;
                    for (pa, pb) in &profiles {
                        let ba = pa[i - 1];
                        let bb = pb[j - 1];
                        if ba == b'-' && bb == b'-' {
                            continue;
                        } else if ba == b'-' || bb == b'-' {
                            s += gap_open_pen + gap_extend_pen;
                        } else {
                            let raw = submat.get_score(ba as char, bb as char);
                            s += raw / 50;
                        }
                    }
                    let cand = score[pk].saturating_add(s);
                    best = cand;
                    bt = 1;
                }
            }

            let mut gi_val = i32::MIN;
            if i > 0 {
                if let Some(pk_score) = idx(i - 1, j) {
                    let from_match = score[pk_score]
                        .saturating_add(gap_open_pen)
                        .saturating_add(gap_extend_pen);
                    let from_gap = gap_i[pk_score].saturating_add(gap_extend_pen);
                    gi_val = from_match.max(from_gap);
                    if gi_val > best {
                        best = gi_val;
                        bt = 2;
                    }
                }
            }

            let mut gj_val = i32::MIN;
            if j > 0 {
                if let Some(pk_score) = idx(i, j - 1) {
                    let from_match = score[pk_score]
                        .saturating_add(gap_open_pen)
                        .saturating_add(gap_extend_pen);
                    let from_gap = gap_j[pk_score].saturating_add(gap_extend_pen);
                    gj_val = from_match.max(from_gap);
                    if gj_val > best {
                        best = gj_val;
                        bt = 3;
                    }
                }
            }

            score[k] = best;
            gap_i[k] = gi_val;
            gap_j[k] = gj_val;
            trace[k] = bt;
        }
    }

    let mut i = n;
    let mut j = m;

    idx(i, j)?;

    let mut map_a = Vec::new();
    let mut map_b = Vec::new();

    while i > 0 || j > 0 {
        let k = match idx(i, j) {
            Some(v) => v,
            None => break,
        };
        let bt = trace[k];
        if bt == 1 {
            if i == 0 || j == 0 {
                break;
            }
            let pi = i - 1;
            let pj = j - 1;
            map_a.push(Some(pi));
            map_b.push(Some(pj));
            i -= 1;
            j -= 1;
        } else if bt == 2 {
            if i == 0 {
                break;
            }
            let pi = i - 1;
            map_a.push(Some(pi));
            map_b.push(None);
            i -= 1;
        } else if bt == 3 {
            if j == 0 {
                break;
            }
            let pj = j - 1;
            map_a.push(None);
            map_b.push(Some(pj));
            j -= 1;
        } else {
            break;
        }
    }

    map_a.reverse();
    map_b.reverse();

    if map_a.len() != map_b.len() || map_a.is_empty() {
        return None;
    }

    let mut left = 0usize;
    let mut right = map_a.len();

    while left < right {
        let a = map_a[left];
        let b = map_b[left];
        if a.is_some() && b.is_some() {
            break;
        }
        left += 1;
    }

    while right > left {
        let a = map_a[right - 1];
        let b = map_b[right - 1];
        if a.is_some() && b.is_some() {
            break;
        }
        right -= 1;
    }

    if left >= right {
        return None;
    }

    let map_a = map_a[left..right].to_vec();
    let map_b = map_b[left..right].to_vec();

    Some((map_a, map_b))
}
