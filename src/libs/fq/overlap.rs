//! Overlap detection between paired reads (BBMerge-compatible).

/// Result of ratio-based overlap detection.
///
/// Mirrors BBMergeOverlapper's `rvector` outputs: `bad` is the integer
/// mismatch count of the winning overlap (`rvector[2]`) and `ambig`
/// corresponds to `rvector[4]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Overlap {
    /// Best insert size (merged-read length), or -1 when none was found.
    pub insert: i32,
    /// Integer mismatch count of the best overlap.
    pub bad: i32,
    /// Whether the best overlap was ambiguous.
    pub ambig: bool,
    /// Best overlap length (number of overlapping bases).
    pub best_overlap: i32,
    /// Best float mismatch score (`bad` in BBMerge).
    pub best_bad: f32,
    /// Best float match score (`good`).
    pub best_good: f32,
    /// Best ratio.
    pub best_ratio: f32,
    /// Second-best insert size.
    pub second_insert: i32,
    /// Second-best overlap length.
    pub second_overlap: i32,
    /// Second-best float mismatch score.
    pub second_bad: f32,
    /// Second-best float match score.
    pub second_good: f32,
    /// Second-best ratio.
    pub second_ratio: f32,
    /// Second-best integer mismatch count.
    pub second_bad_int: i32,
}

/// BBMergeOverlapper.mateByOverlapRatioJava (no-quality path).
///
/// `b` must already be reverse-complemented. `g_incr`/`b_incr` are the
/// increments for matching/mismatching bases (bbmerge uses 0.95/0.95).
#[allow(clippy::too_many_arguments)]
pub fn mate_by_overlap_ratio(
    a: &[u8],
    b: &[u8],
    min_overlap0: usize,
    min_overlap: usize,
    min_insert0: usize,
    min_insert: usize,
    max_ratio: f32,
    min_second_ratio: f32,
    margin: f32,
    offset: f32,
    g_incr: f32,
    b_incr: f32,
) -> Overlap {
    // Java: minOverlap=max(4, minOverlap0, minOverlap);
    //       minOverlap0=mid(4, minOverlap0, minOverlap);
    let min_overlap = min_overlap.max(min_overlap0).max(4);
    let min_overlap0 = min_overlap0.clamp(4, min_overlap);

    let alen = a.len() as i32;
    let blen = b.len() as i32;
    let min_length = alen.min(blen) as usize;
    let n = b'N';
    let mut max_ratio = max_ratio;
    let x = find_best_ratio(
        a,
        b,
        min_overlap0,
        min_overlap,
        min_insert,
        max_ratio,
        offset,
        g_incr,
        b_incr,
    );
    if x > max_ratio {
        // rvector[2]=minLength, rvector[4]=0
        return Overlap {
            insert: -1,
            bad: min_length as i32,
            ambig: false,
            best_overlap: -1,
            best_bad: min_length as f32,
            best_good: 0.0,
            best_ratio: 1.0,
            second_insert: 0,
            second_overlap: 0,
            second_bad: 0.0,
            second_good: 0.0,
            second_ratio: 1.0,
            second_bad_int: -1,
        };
    }
    max_ratio = max_ratio.min(x);
    let margin2 = (margin + offset) / min_length as f32;
    let mut best_insert = -1i32;
    let mut best_overlap = -1i32;
    let mut best_bad = min_length as f32;
    let mut best_ratio = 1f32;
    let mut best_ambig = false;
    let mut second_insert = 0i32;
    let mut second_overlap = 0i32;
    let mut second_bad = 0f32;
    let mut second_best_ratio = 1f32;
    let mut best_bad_int = -1i32;
    let mut second_bad_int = -1i32;
    let extra_mult = 1.2f32;
    let extra_badlimit = 20f32;
    let largest = alen + blen - min_overlap0 as i32;
    let smallest = min_insert0 as i32;
    let mut insert = largest;
    while insert >= smallest {
        let istart = if insert <= blen { 0 } else { insert - blen };
        let jstart = if insert >= blen { 0 } else { blen - insert };
        let overlap_length = (alen - istart).min(blen - jstart).min(insert);
        let badlimit = extra_mult * (best_ratio.min(max_ratio) * margin * overlap_length as f32)
            + 1.0
            + extra_badlimit;
        let mut good = 0f32;
        let mut bad = 0f32;
        let mut bad_int = 0i32;
        let imax = istart + overlap_length;
        let mut i = istart;
        let mut j = jstart;
        while i < imax && bad <= badlimit {
            let ca = a[i as usize];
            let cb = b[j as usize];
            if ca == cb {
                if ca != n {
                    good += g_incr;
                }
            } else {
                bad += b_incr;
                bad_int += 1;
            }
            i += 1;
            j += 1;
        }
        if bad <= badlimit {
            if bad == 0.0 && good > min_overlap0 as f32 && good < min_overlap as f32 {
                // rvector[2]=bestBadInt (=-1), rvector[4]=1
                return Overlap {
                    insert: -1,
                    bad: best_bad_int,
                    ambig: true,
                    best_overlap,
                    best_bad,
                    best_good: 0.0,
                    best_ratio,
                    second_insert,
                    second_overlap,
                    second_bad,
                    second_good: 0.0,
                    second_ratio: second_best_ratio,
                    second_bad_int,
                };
            }
            let ratio = (bad + offset) / overlap_length as f32;
            if ratio < best_ratio * margin {
                let this_ambig = ratio * margin >= best_ratio || good < min_overlap as f32;
                if ratio < best_ratio {
                    second_best_ratio = best_ratio;
                    second_insert = best_insert;
                    second_overlap = best_overlap;
                    second_bad = best_bad;
                    second_bad_int = best_bad_int;
                    best_insert = insert;
                    best_overlap = overlap_length;
                    best_bad = bad;
                    best_ratio = ratio;
                    best_bad_int = bad_int;
                } else if ratio < second_best_ratio {
                    second_best_ratio = ratio;
                    second_insert = insert;
                    second_overlap = overlap_length;
                    second_bad = bad;
                    second_bad_int = bad_int;
                }
                best_ambig = this_ambig;
                if (best_ambig && best_ratio < margin2) || second_best_ratio < min_second_ratio {
                    // rvector[2]=bestBadInt, rvector[4]=1
                    return Overlap {
                        insert: -1,
                        bad: best_bad_int,
                        ambig: true,
                        best_overlap,
                        best_bad,
                        best_good: 0.0,
                        best_ratio,
                        second_insert,
                        second_overlap,
                        second_bad,
                        second_good: 0.0,
                        second_ratio: second_best_ratio,
                        second_bad_int,
                    };
                }
            }
        }
        insert -= 1;
    }
    if second_best_ratio < min_second_ratio {
        best_ambig = true;
    }
    if !best_ambig && best_ratio > max_ratio {
        best_insert = -1;
    }
    Overlap {
        insert: best_insert,
        bad: best_bad_int,
        ambig: best_ambig,
        best_overlap,
        best_bad,
        best_good: 0.0,
        best_ratio,
        second_insert,
        second_overlap,
        second_bad,
        second_good: 0.0,
        second_ratio: second_best_ratio,
        second_bad_int,
    }
}

/// findBestRatio: fast pre-screen for overlap detection.
#[allow(clippy::too_many_arguments)]
fn find_best_ratio(
    a: &[u8],
    b: &[u8],
    min_overlap0: usize,
    min_overlap: usize,
    min_insert: usize,
    max_ratio: f32,
    offset: f32,
    g_incr: f32,
    b_incr: f32,
) -> f32 {
    let alen = a.len() as i32;
    let blen = b.len() as i32;
    let n = b'N';
    let mut best_ratio = max_ratio + 0.0001;
    let halfmax = max_ratio * 0.5;
    let largest = alen + blen - min_overlap as i32;
    let smallest = min_insert as i32;
    let mut insert = largest;
    while insert >= smallest {
        let istart = if insert <= blen { 0 } else { insert - blen };
        let jstart = if insert >= blen { 0 } else { blen - insert };
        let overlap_length = (alen - istart).min(blen - jstart).min(insert);
        let badlimit = best_ratio * overlap_length as f32 + 20.0;
        let mut good = 0f32;
        let mut bad = 0f32;
        let imax = istart + overlap_length;
        let mut i = istart;
        let mut j = jstart;
        while i < imax && bad <= badlimit {
            let ca = a[i as usize];
            let cb = b[j as usize];
            if ca == cb {
                if ca != n {
                    good += g_incr;
                }
            } else {
                bad += b_incr;
            }
            i += 1;
            j += 1;
        }
        if bad <= badlimit {
            if bad == 0.0 && good > min_overlap0 as f32 && good < min_overlap as f32 {
                return 100.0;
            }
            let ratio = (bad + offset) / overlap_length as f32;
            if ratio < best_ratio {
                best_ratio = ratio;
                if good >= min_overlap as f32 && ratio < halfmax {
                    return best_ratio;
                }
            }
        }
        insert -= 1;
    }
    best_ratio
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_overlap_finds_insert() {
        // r1 = 20 nt, b = reverse-complemented read2 (r1's last 12 nt).
        let r1 = b"GCTAAAGACAATTACATAAC";
        let b = b"ACAATTACATAAC";
        let o = mate_by_overlap_ratio(r1, b, 4, 8, 5, 20, 0.09, 0.1, 5.5, 0.55, 0.95, 0.95);
        assert!(o.insert > 0);
        assert!(!o.ambig);
        assert_eq!(o.bad, 0);
        assert_eq!(o.insert, 20);
    }

    #[test]
    fn clamps_min_overlap_like_java() {
        // min_overlap0=0 must be raised to 4 before the scan, so the largest
        // tested insert is alen+blen-4 rather than alen+blen.
        let a = b"ATACACGTCAGCACGAAACTTGTT";
        let b = b"CACGAAACTTGTT"; // reverse-complemented read2 tail
        let o = mate_by_overlap_ratio(a, b, 0, 8, 5, 24, 0.09, 0.1, 5.5, 0.55, 0.95, 0.95);
        assert!(o.insert > 0);
        assert_eq!(o.bad, 0);
    }
}
