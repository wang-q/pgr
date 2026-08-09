//! Quality-aware read error detection (quorum anchor+extend signals).

use super::quality::QualityTable;

/// Why a read was flagged as bad (quorum sub/trunc events, no correction).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadError {
    /// No high-quality anchor k-mer found.
    NoAnchor,
    /// A position with no continuation (quorum truncation).
    Truncation,
    /// A position whose base would be corrected (quorum substitution).
    Substitution,
}

/// Thresholds mirroring quorum's error-correction parameters.
#[derive(Debug, Clone)]
pub struct CheckParams {
    /// K-mer length.
    pub k: usize,
    /// Bases to skip before searching for an anchor.
    pub skip: usize,
    /// Consecutive anchor k-mers required (quorum `--good`).
    pub good: usize,
    /// Minimum count for a high-quality anchor k-mer (quorum `--anchor-count`).
    pub anchor_count: usize,
    /// Count above which a base is trusted before the cutoff check.
    pub min_count: u64,
    /// Trusted count for keeping the current base (quorum `--cutoff`).
    pub cutoff: u64,
    /// Prior error rate / 3 used in the Poisson collision test.
    pub collision_prob: f64,
    /// Poisson probability threshold below which the base is kept.
    pub poisson_threshold: f64,
}

impl Default for CheckParams {
    fn default() -> Self {
        Self {
            k: 17,
            skip: 0,
            good: 1,
            anchor_count: 1,
            min_count: 1,
            cutoff: 4,
            collision_prob: 0.01 / 3.0,
            poisson_threshold: 1.0e-06,
        }
    }
}

/// Check one read; `Ok(())` keeps it, `Err` gives the discard reason.
///
/// Mirrors quorum: find a run of `good` high-quality anchor k-mers, then
/// extend in both directions and flag any position that quorum would
/// substitute or truncate. Only the presence of an error is reported; no
/// corrected sequence is produced.
pub fn check_read(table: &QualityTable, seq: &[u8], p: &CheckParams) -> Result<(), ReadError> {
    let Some((start, kx, kxr)) = find_anchor(table, seq, p) else {
        return Err(ReadError::NoAnchor);
    };
    extend(table, seq, start + p.k, kx, kxr, p)?;

    // Mirror the read to extend leftwards with the same forward logic.
    let rev: Vec<u8> = seq.iter().rev().copied().collect();
    let rev_start = seq.len() - start - p.k;
    let Some((rkx, rkxr)) = roll_window(&rev, rev_start, p.k) else {
        return Err(ReadError::Truncation);
    };
    extend(table, &rev, rev_start + p.k, rkx, rkxr, p)
}

/// First position whose window completes `good` consecutive high-quality
/// k-mers with count >= `anchor_count`; returns the window start and its
/// rolling keys.
fn find_anchor(table: &QualityTable, seq: &[u8], p: &CheckParams) -> Option<(usize, u128, u128)> {
    if seq.len() < p.k {
        return None;
    }
    let (kmask, rc_top) = masks(p.k);
    let codes = super::base_codes();
    let mut kx: u128 = 0;
    let mut kxr: u128 = 0;
    let mut valid = 0usize;
    let mut found = 0usize;
    for (i, &b) in seq.iter().enumerate().skip(p.skip) {
        let code = codes[b as usize];
        if code == 4 {
            kx = 0;
            kxr = 0;
            valid = 0;
            found = 0;
            continue;
        }
        kx = ((kx << 2) | code as u128) & kmask;
        kxr = (kxr >> 2) | (((3 - code) as u128) << rc_top);
        valid += 1;
        if valid >= p.k {
            let start = i + 1 - p.k;
            // quorum get_val(): only high-quality k-mers count as anchors.
            let anchor_ok =
                matches!(table.get(kx.min(kxr)), Some((c, 1)) if c as usize >= p.anchor_count);
            found = if anchor_ok { found + 1 } else { 0 };
            if found >= p.good {
                return Some((start, kx, kxr));
            }
        }
    }
    None
}

/// Rolling 2-bit keys of `seq[start..start+k]`; `None` if it contains N.
fn roll_window(seq: &[u8], start: usize, k: usize) -> Option<(u128, u128)> {
    let (kmask, rc_top) = masks(k);
    let codes = super::base_codes();
    let mut kx: u128 = 0;
    let mut kxr: u128 = 0;
    for &b in &seq[start..start + k] {
        let code = codes[b as usize];
        if code == 4 {
            return None;
        }
        kx = ((kx << 2) | code as u128) & kmask;
        kxr = (kxr >> 2) | (((3 - code) as u128) << rc_top);
    }
    Some((kx, kxr))
}

fn masks(k: usize) -> (u128, u32) {
    let kmask = if 2 * k >= 128 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    };
    (kmask, (2 * k - 2) as u32)
}

/// Extend rightwards from `start` (just past the anchor window), flagging
/// positions quorum would correct. `kx`/`kxr` are the window's rolling keys.
fn extend(
    table: &QualityTable,
    seq: &[u8],
    start: usize,
    mut kx: u128,
    mut kxr: u128,
    p: &CheckParams,
) -> Result<(), ReadError> {
    let (kmask, rc_top) = masks(p.k);
    let codes = super::base_codes();
    for &b in &seq[start..] {
        let code = codes[b as usize];
        if code == 4 {
            // quorum would substitute or truncate an N; flag the read.
            return Err(ReadError::Substitution);
        }
        kx = ((kx << 2) | code as u128) & kmask;
        kxr = (kxr >> 2) | (((3 - code) as u128) << rc_top);
        let (counts, ucode, level, count) = best_alternatives(table, kx, kxr, rc_top);
        if count == 0 {
            return Err(ReadError::Truncation);
        }
        if count == 1 {
            if ucode != code as usize {
                return Err(ReadError::Substitution);
            }
            continue;
        }
        let ori = counts[code as usize];
        if ori > p.min_count {
            if ori >= p.cutoff {
                continue;
            }
            let p_err = (counts[0] + counts[1] + counts[2] + counts[3]) as f64 * p.collision_prob;
            if poisson_term(p_err, ori) < p.poisson_threshold {
                continue;
            }
            return Err(ReadError::Substitution);
        }
        if level == 0 && ori == 0 {
            return Err(ReadError::Truncation);
        }
        return Err(ReadError::Substitution);
    }
    Ok(())
}

/// Per-base alternative counts at the highest quality level (quorum
/// `get_best_alternatives`): replaces the current (newest) base with A/C/G/T
/// and queries each canonical k-mer.
fn best_alternatives(
    table: &QualityTable,
    kx: u128,
    kxr: u128,
    rc_top: u32,
) -> ([u64; 4], usize, u8, usize) {
    let mut counts = [0u64; 4];
    let mut level = 0u8;
    let mut ucode = 0usize;
    let mut count = 0usize;
    for i in 0..4u128 {
        let kx2 = (kx & !3) | i;
        let kxr2 = (kxr & !(3u128 << rc_top)) | ((3 - i) << rc_top);
        if let Some((c, q)) = table.get(kx2.min(kxr2)) {
            if q >= level {
                if q > level && count > 0 {
                    counts = [0u64; 4];
                    count = 0;
                }
                counts[i as usize] = c as u64;
                ucode = i as usize;
                level = q;
                count += 1;
            }
        }
    }
    (counts, ucode, level, count)
}

/// Poisson term used by quorum's collision test (exact for i<11, Stirling
/// approximation beyond).
fn poisson_term(lambda: f64, i: u64) -> f64 {
    const FACTS: [f64; 11] = [
        1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0,
    ];
    const TAU: f64 = 6.283185307179583;
    if i < 11 {
        (-lambda).exp() * lambda.powf(i as f64) / FACTS[i as usize]
    } else {
        (-lambda + i as f64).exp() * (lambda / i as f64).powf(i as f64) / (TAU * i as f64).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::kmer::quality::build_table;

    fn params(k: usize) -> CheckParams {
        CheckParams {
            k,
            ..Default::default()
        }
    }

    #[test]
    fn clean_read_passes() {
        let read = b"ACGTACGTACGTACGTACGT".to_vec();
        let quals = vec![b"I".repeat(20); 50];
        let seqs = vec![read.clone(); 50];
        let table = build_table(&seqs, &quals, 8, 38, 127);
        assert_eq!(check_read(&table, &read, &params(8)), Ok(()));
    }

    #[test]
    fn mutated_read_is_flagged() {
        let clean = b"ACGTACGTACGTACGTACGT".to_vec();
        let quals = vec![b"I".repeat(20); 50];
        let seqs = vec![clean.clone(); 50];
        let table = build_table(&seqs, &quals, 8, 38, 127);

        // Position 10 G->C: windows covering it are absent from the table.
        let mut mutated = clean.clone();
        mutated[10] = b'C';
        assert_eq!(
            check_read(&table, &mutated, &params(8)),
            Err(ReadError::Substitution)
        );
    }

    #[test]
    fn all_low_quality_reads_have_no_anchor() {
        let read = b"ACGTACGTACGTACGTACGT".to_vec();
        let low_quals = vec![b"#".repeat(20); 50];
        let seqs = vec![read.clone(); 50];
        let table = build_table(&seqs, &low_quals, 8, 38, 127);
        assert_eq!(
            check_read(&table, &read, &params(8)),
            Err(ReadError::NoAnchor)
        );
    }

    #[test]
    fn uncovered_region_truncates() {
        // Training covers a 24 bp mixed sequence; a read extending past it
        // with TTTT hits k-mers absent from the table (no continuation).
        let covered = b"ACGTACGTACGTGGGGCCCCAAAA".to_vec();
        let seqs = vec![covered.clone(); 50];
        let quals = vec![b"I".repeat(24); 50];
        let table = build_table(&seqs, &quals, 8, 38, 127);

        let mut long = covered.clone();
        long.extend_from_slice(b"TTTT");
        assert_eq!(
            check_read(&table, &long, &params(8)),
            Err(ReadError::Truncation)
        );
    }

    #[test]
    fn short_read_has_no_anchor() {
        // Reads shorter than k cannot contain any k-mer: no anchor, no panic.
        let seqs = vec![b"ACGTACGTACGT".to_vec(); 10];
        let quals = vec![b"I".repeat(12); 10];
        let table = build_table(&seqs, &quals, 8, 38, 127);
        assert_eq!(
            check_read(&table, b"ACG".as_slice(), &params(8)),
            Err(ReadError::NoAnchor)
        );
    }
}
