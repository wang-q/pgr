//! Quality-weighted k-mer histogram (quorum `histo_mer_database` equivalent).

use crate::libs::ds::radix_sort::radix_sort_u128_par;
use rayon::prelude::*;
use std::io::Write;

/// Histogram length and count cap of quorum's `histo_mer_database` (0..=1000).
pub const QUAL_HLEN: usize = 1001;

/// Quality-biased canonical k-mer count table (quorum `hash_with_quality`).
///
/// A k-mer counts as high quality iff all `k` bases of its window score at
/// least `thresh` (quality chars, e.g. Phred+33 ASCII). Per k-mer, the final
/// count is the number of high-quality occurrences when any exist, otherwise
/// the number of low-quality occurrences; the quality flag mirrors that
/// choice. The aggregation is order-independent, matching `hash_with_quality`
/// (a low-quality occurrence never raises the count once a high-quality one
/// was seen, and a high-quality occurrence resets a low-quality count).
#[derive(Debug, Clone)]
pub struct QualityTable {
    /// K-mer length (bp).
    pub k: usize,
    /// Sorted canonical 2-bit k-mers.
    pub keys: Vec<u128>,
    /// Quorum-biased counts (capped by the build `count_cap`).
    pub counts: Vec<u32>,
    /// Quality flag: 1 if any high-quality occurrence was seen.
    pub qualities: Vec<u8>,
}

impl QualityTable {
    /// Count and quality of `key`, or `None` when absent from the table.
    pub fn get(&self, key: u128) -> Option<(u32, u8)> {
        let idx = self.keys.partition_point(|&x| x < key);
        if idx < self.keys.len() && self.keys[idx] == key {
            Some((self.counts[idx], self.qualities[idx]))
        } else {
            None
        }
    }
}

/// Build a quality-biased count table from reads.
pub fn build_table(
    seqs: &[Vec<u8>],
    quals: &[Vec<u8>],
    k: usize,
    thresh: u8,
    count_cap: u64,
) -> QualityTable {
    let per_read: Vec<Vec<(u128, u8)>> = seqs
        .par_iter()
        .zip(quals)
        .map(|(seq, qual)| {
            let mut v = Vec::new();
            quality_keys(seq, qual, k, thresh, |key, high| v.push((key, high)));
            v
        })
        .collect();
    let n: usize = per_read.iter().map(Vec::len).sum();
    let mut keys: Vec<u128> = Vec::with_capacity(n);
    let mut flags: Vec<u8> = Vec::with_capacity(n);
    for v in per_read {
        for (key, high) in v {
            keys.push(key);
            flags.push(high);
        }
    }
    radix_sort_u128_par(&mut keys, &mut flags, 2 * k as u32);

    let mut out_keys = Vec::with_capacity(keys.len());
    let mut counts = Vec::with_capacity(keys.len());
    let mut qualities = Vec::with_capacity(keys.len());
    let mut i = 0usize;
    while i < keys.len() {
        let key = keys[i];
        let mut n_high = 0u64;
        let mut n_low = 0u64;
        while i < keys.len() && keys[i] == key {
            if flags[i] != 0 {
                n_high += 1;
            } else {
                n_low += 1;
            }
            i += 1;
        }
        let (count, quality) = if n_high > 0 {
            (n_high.min(count_cap), 1u8)
        } else {
            (n_low.min(count_cap), 0u8)
        };
        out_keys.push(key);
        counts.push(count as u32);
        qualities.push(quality);
    }
    QualityTable {
        k,
        keys: out_keys,
        counts,
        qualities,
    }
}

/// Histogram bins from a quality table (quorum `histo_mer_database` format).
pub fn histogram(table: &QualityTable) -> Vec<[u64; 2]> {
    let mut hist = vec![[0u64; 2]; QUAL_HLEN];
    for (&count, &quality) in table.counts.iter().zip(&table.qualities) {
        hist[(count as usize).min(QUAL_HLEN - 1)][quality as usize] += 1;
    }
    hist
}

/// Emit `(canonical key, all-k-bases-high-quality)` for every N-free window.
///
/// Mirrors quorum's `quality_mer_counter`: non-ACGT bases reset both the
/// low-quality and high-quality stretches, a base below `thresh` resets only
/// the high-quality stretch, and a window is high quality iff `high_len`
/// reached `k` at its last base.
fn quality_keys(seq: &[u8], qual: &[u8], k: usize, thresh: u8, mut emit: impl FnMut(u128, u8)) {
    if seq.len() < k {
        return;
    }
    let codes = super::base_codes();
    let kmask = if 2 * k >= 128 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    };
    let rc_top = (2 * k - 2) as u32;
    let mut kx: u128 = 0;
    let mut kxr: u128 = 0;
    let mut low_len = 0usize;
    let mut high_len = 0usize;
    for (&b, &q) in seq.iter().zip(qual) {
        let code = codes[b as usize];
        if code == 4 {
            kx = 0;
            kxr = 0;
            low_len = 0;
            high_len = 0;
            continue;
        }
        kx = ((kx << 2) | code as u128) & kmask;
        kxr = (kxr >> 2) | (((3 - code) as u128) << rc_top);
        low_len += 1;
        if q >= thresh {
            high_len += 1;
        } else {
            high_len = 0;
        }
        if low_len >= k {
            emit(kx.min(kxr), u8::from(high_len >= k));
        }
    }
}

/// Write the histogram in quorum's `histo_mer_database` format:
/// `count n_lowq n_highq` per non-empty count bin.
pub fn write_hist(w: &mut impl Write, hist: &[[u64; 2]]) -> anyhow::Result<()> {
    for (i, [lo, hi]) in hist.iter().enumerate() {
        if *lo > 0 || *hi > 0 {
            writeln!(w, "{i} {lo} {hi}")?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::StdRng, Rng, SeedableRng};
    use std::collections::HashMap;

    /// Deterministic pseudo-random DNA block.
    fn random_dna(len: usize, rng: &mut StdRng) -> Vec<u8> {
        (0..len)
            .map(|_| *b"ACGT".get(rng.random_range(0..4)).unwrap())
            .collect()
    }

    /// Random quality chars spanning both sides of the default threshold 38.
    fn random_quals(len: usize, rng: &mut StdRng) -> Vec<u8> {
        (0..len).map(|_| rng.random_range(33..74) as u8).collect()
    }

    /// Strict sequential simulation of quorum `hash_with_quality::add`.
    ///
    /// Applies the nval update rule in event order with a count cap
    /// (`max_val = 2^bits - 1`), mirroring the C implementation line by line.
    fn quorum_sequential(events: &[(u128, u8)], count_cap: u64) -> HashMap<u128, (u64, u8)> {
        let mut vals: HashMap<u128, u64> = HashMap::new();
        for &(key, quality) in events {
            let nval = vals.entry(key).or_insert(0);
            *nval = if (*nval & 1) < u64::from(quality) {
                3
            } else if (*nval >> 1) == count_cap || (*nval & 1) > u64::from(quality) {
                *nval
            } else {
                *nval + 2
            };
        }
        vals.into_iter()
            .map(|(k, v)| (k, (v >> 1, (v & 1) as u8)))
            .collect()
    }

    #[test]
    fn high_quality_occurrences_dominate_count() {
        // k=4, thresh 38 (Phred+33 '!'+5): 'I' (73) is high, '#' (35) is low.
        let seqs = vec![
            b"ACGTACGT".to_vec(), // 5 windows, all high
            b"ACGTACGT".to_vec(), // same k-mers, all high -> count 2
        ];
        let quals = vec![b"IIIIIIII".to_vec(), b"IIIIIIII".to_vec()];
        let table = build_table(&seqs, &quals, 4, 38, u64::MAX);
        let hist = histogram(&table);
        // Canonical set is {ACGT (palindrome), CGTA (rc of TACG), GTAC
        // (palindrome)}: ACGT and CGTA occur 4x, GTAC 2x, all high.
        assert_eq!(hist[0], [0, 0]);
        assert_eq!(hist[4][1], 2);
        assert_eq!(hist[2][1], 1);
        assert_eq!(hist[4][0], 0);
    }

    #[test]
    fn low_quality_never_pollutes_high_count() {
        // Same k-mer seen 3x low then 1x high: count must be 1 (high), not 4.
        let seqs = vec![
            b"ACGTACGT".to_vec(), // 5 low windows
            b"ACGTACGT".to_vec(), // 5 high windows
        ];
        let quals = vec![b"########".to_vec(), b"IIIIIIII".to_vec()];
        let table = build_table(&seqs, &quals, 4, 38, u64::MAX);
        let hist = histogram(&table);
        // Each canonical k-mer has both low and high occurrences; the count
        // must equal the high-occurrence count only.
        assert_eq!(hist[2][1], 2); // ACGT, CGTA: 2 high occurrences each
        assert_eq!(hist[1][1], 1); // GTAC: 1 high occurrence
        assert_eq!(hist[5][0], 0, "low occurrences must not add to the count");
    }

    #[test]
    fn n_splits_stretches() {
        // N splits the window; low-quality base splits only the high stretch.
        let seqs = vec![b"ACGTNNACGT".to_vec()];
        let quals = vec![b"##########".to_vec()];
        let table = build_table(&seqs, &quals, 4, 38, u64::MAX);
        let hist = histogram(&table);
        // Only the two N-free 4-base windows (left and right of the N run)
        // emit; both are the same canonical k-mer (ACGT) at low quality.
        assert_eq!(hist[2][0], 1);
        assert_eq!(hist[2][1], 0);
    }

    #[test]
    fn write_hist_matches_quorum_format() {
        let mut hist = vec![[0u64; 2]; QUAL_HLEN];
        hist[1][0] = 3;
        hist[4][1] = 7;
        let mut out = Vec::new();
        write_hist(&mut out, &hist).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "1 3 0\n4 0 7\n");
    }

    #[test]
    fn count_caps_at_bits_limit() {
        // quorum create_database default `-b 7` caps counts at 127.
        let seqs: Vec<Vec<u8>> = (0..130).map(|_| b"ACGTACGTACGTACGT".to_vec()).collect();
        let quals: Vec<Vec<u8>> = (0..130).map(|_| b"IIIIIIIIIIIIIIII".to_vec()).collect();
        let table = build_table(&seqs, &quals, 4, 38, 127);
        let hist = histogram(&table);
        assert_eq!(hist[127][1], 3); // ACGT/CGTA capped, GTAC at 130->127
        assert_eq!(hist[130][1], 0);
        // Direct query agrees with the histogram bins (ACGT key = 0b00011011).
        assert_eq!(table.get(0b00011011), Some((127, 1)));
    }

    #[test]
    fn aggregation_matches_sequential_quorum_semantics() {
        // Random reads with random qualities: the order-independent aggregate
        // must agree with a strict sequential simulation of quorum's add().
        let mut rng = StdRng::seed_from_u64(20260809);
        let cap = 127;
        for k in [4usize, 8, 17] {
            for _ in 0..20 {
                let n_reads = rng.random_range(1..30);
                let mut seqs = Vec::new();
                let mut quals = Vec::new();
                for _ in 0..n_reads {
                    let len = rng.random_range(k..80);
                    seqs.push(random_dna(len, &mut rng));
                    quals.push(random_quals(len, &mut rng));
                }
                let mut events = Vec::new();
                for (seq, qual) in seqs.iter().zip(&quals) {
                    quality_keys(seq, qual, k, 38, |key, high| {
                        events.push((key, high));
                    });
                }
                let sim = quorum_sequential(&events, cap);
                let table = build_table(&seqs, &quals, k, 38, cap);
                let hist = histogram(&table);
                let mut expect = vec![[0u64; 2]; QUAL_HLEN];
                for &(count, quality) in sim.values() {
                    expect[(count as usize).min(QUAL_HLEN - 1)][quality as usize] += 1;
                }
                assert_eq!(hist, expect, "k={k} reads={n_reads}");
            }
        }
    }
}
