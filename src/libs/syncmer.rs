//! Closed syncmer sampling for DNA and protein sequences.
//!
//! Ported from Richard Durbin's syng `seqhash.c`, implementing closed
//! syncmers per Edgar (2021, <https://peerj.com/articles/10805/>). A closed
//! syncmer is a window of `window` s-mers whose minimal s-mer hash falls at
//! the first or last position; this guarantees a sparse but complete cover
//! of any sequence (average depth ~2). See `notes/references/syng.md`.

use std::collections::VecDeque;

use crate::libs::hash::{Hasher, RapidHash};

/// Parameters for syncmer sampling.
///
/// `smer` is the small k-mer length used for hashing (syng's `k`).
/// `window` is the number of s-mers per syncmer window (syng's `w`).
/// A syncmer spans `smer + window - 1` bases and contains `window` s-mers;
/// it is emitted iff the minimal s-mer hash in the window is at the first
/// or last s-mer. `seed` seeds the DNA rolling hash (ignored for protein).
#[derive(Debug, Clone, Copy)]
pub struct SyncmerParams {
    pub smer: usize,
    pub window: usize,
    pub seed: u64,
}

impl SyncmerParams {
    /// Default DNA parameters matching syng (smer=8, window=55, seed=7).
    pub fn default_dna() -> Self {
        Self {
            smer: 8,
            window: 55,
            seed: 7,
        }
    }

    /// Validate parameters; smer must fit 2-bit-per-base in a u64.
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.smer == 0 {
            anyhow::bail!("syncmer smer must be positive");
        }
        if self.smer >= 32 {
            anyhow::bail!("syncmer smer must be < 32, got {}", self.smer);
        }
        if self.window == 0 {
            anyhow::bail!("syncmer window must be positive");
        }
        // Cap the window so the pgi build path's `(window + 2).next_power_of_two()`
        // ring buffer cannot overflow (an extreme `--window` like `usize::MAX`
        // would otherwise panic on `window + 2`) or force a multi-GB allocation.
        // Real windows are tiny (default 5); 1,000,000 is a generous upper bound.
        if self.window > 1_000_000 {
            anyhow::bail!("syncmer window must be <= 1000000, got {}", self.window);
        }
        Ok(())
    }
}

/// Alphabet-agnostic core: given the s-mer hash stream of a sequence,
/// return the endpoint index of each closed syncmer window.
///
/// A window `[i, i+window)` is a closed syncmer iff the window minimum
/// hash appears at index `i` (first) or `i+window-1` (last). This returns
/// that endpoint index — preferring `i` (first) when both endpoints are
/// minimal — so that a sequence and its reverse complement yield the same
/// hash set (required for Mash/Jaccard). syng emits the window-start s-mer
/// instead, which is fine for graph paths but not strand-symmetric.
#[cfg(test)]
fn closed_syncmers_from_hashes(hashes: &[u64], window: usize) -> Vec<usize> {
    closed_syncmers_stream(window, hashes.iter().copied().map(|h| (h, ())))
        .into_iter()
        .map(|(pos, _, _)| pos)
        .collect()
}

/// Streaming closed-syncmer detector over an on-the-fly s-mer hash stream.
///
/// O(window) memory: a monotonic deque of (position, hash) for the window
/// minimum plus a ring buffer of the last `window` raw entries for endpoint
/// checks. Emits `(endpoint_index, hash, extra)` per closed syncmer window;
/// the extra (e.g. the canonical-strand flag) travels with the endpoint
/// s-mer. Produces exactly the same output as [`closed_syncmers_from_hashes`]
/// (which is now a thin wrapper over this).
fn closed_syncmers_stream<E: Copy>(
    window: usize,
    hashes: impl IntoIterator<Item = (u64, E)>,
) -> Vec<(usize, u64, E)> {
    let mut out = Vec::new();
    if window == 0 {
        return out;
    }
    // Monotonic deque of positions with non-decreasing hashes; front is the
    // current window minimum. Ring buffer holds the last `window` raw entries
    // so the window-start hash/extra can be recovered for endpoint checks.
    let mut dq: VecDeque<(usize, u64)> = VecDeque::new();
    let mut ring: VecDeque<(u64, E)> = VecDeque::new();
    for (j, (h, e)) in hashes.into_iter().enumerate() {
        while let Some(&(_, bh)) = dq.back() {
            if bh <= h {
                break;
            }
            dq.pop_back();
        }
        dq.push_back((j, h));
        ring.push_back((h, e));
        if ring.len() > window {
            ring.pop_front();
        }

        if j >= window - 1 {
            let start = j + 1 - window;
            while let Some(&(fp, _)) = dq.front() {
                if fp < start {
                    dq.pop_front();
                } else {
                    break;
                }
            }
            let min_val = dq.front().expect("deque non-empty within full window").1;
            // Closed syncmer: the minimum appears at the first or last s-mer.
            // Check value (not argmin position) so ties are symmetric under
            // reversal; prefer the window start on ties.
            let (sh, se) = *ring.front().expect("ring holds a full window");
            if sh == min_val {
                out.push((start, sh, se));
            } else if h == min_val {
                out.push((j, h, e));
            }
        }
    }
    out
}

/// 2-bit encode (A/C/G/T/U->0..3; N/IUPAC/invalid->0=a) via `nt::NT_VAL`.
fn encode_base(b: u8) -> u64 {
    let v = crate::libs::nt::NT_VAL[b as usize];
    if v <= 3 {
        v as u64
    } else {
        0
    }
}

/// Splitmix64-style deterministic pseudo-random odd factor from a seed.
/// syng uses libc `random()`; we use a portable equivalent. Exact factor
/// values differ from syng but are deterministic and uniformly distributed.
pub(crate) fn hash_factor(seed: u64) -> u64 {
    let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    (z ^ (z >> 31)) | 1 // odd, like syng's | 0x01
}

/// Compute closed syncmers over a DNA sequence with O(window) memory.
///
/// Returns `(canonical_hash, pos, is_forward)` per syncmer, where the
/// hash/strand are those of the window's minimal s-mer. syng's
/// `syncmerNext` emits the window-start s-mer for graph path construction;
/// pgr emits the min s-mer instead so a sequence and its reverse complement
/// yield the same hash set (required for Mash/Jaccard). The canonical hash
/// is `min(kHash(h), kHash(hRC))`.
pub fn syncmer_dna(seq: &[u8], params: &SyncmerParams) -> anyhow::Result<Vec<(u64, usize, bool)>> {
    let w = params.window;
    let syncs = closed_syncmers_stream(w, dna_canonical_hashes(seq, params)?);
    Ok(syncs
        .into_iter()
        .map(|(pos, h, fwd)| (h, pos, fwd))
        .collect())
}

/// Streaming canonical hash iterator over the s-mers of `seq`.
///
/// Yields `(canonical_hash, is_forward)` per s-mer without materializing the
/// hash array; `syncmer_dna` feeds it into the O(window) syncmer detector.
fn dna_canonical_hashes<'a>(
    seq: &'a [u8],
    params: &SyncmerParams,
) -> anyhow::Result<impl Iterator<Item = (u64, bool)> + 'a> {
    params.validate()?;
    let k = params.smer;
    let n = seq.len();
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let shift: u32 = (64 - 2 * k) as u32;
    let factor = hash_factor(params.seed);
    let pattern_rc: [u64; 4] = std::array::from_fn(|i| ((3 - i) as u64) << (2 * (k - 1)));
    let k_hash = move |x: u64| x.wrapping_mul(factor) >> shift;

    // First s-mer (positions 0..k).
    let mut h: u64 = 0;
    let mut h_rc: u64 = 0;
    for &byte in seq.iter().take(k) {
        let b = encode_base(byte);
        h = (h << 2) | b;
        h_rc = (h_rc >> 2) | pattern_rc[b as usize];
    }

    let mut next_pos = k;
    let mut first = true;
    let exhausted = n < k;
    let hasher = move || {
        if exhausted {
            return None;
        }
        if first {
            first = false;
        } else {
            let byte = *seq.get(next_pos)?;
            next_pos += 1;
            let b = encode_base(byte);
            h = ((h << 2) & mask) | b;
            h_rc = (h_rc >> 2) | pattern_rc[b as usize];
        }
        let (hf, hr) = (k_hash(h), k_hash(h_rc));
        Some((if hf < hr { hf } else { hr }, hf < hr))
    };
    Ok(std::iter::from_fn(hasher))
}

/// Compute closed syncmers over a protein sequence with O(window) memory.
///
/// Uses the provided byte-string hasher on each s-mer; no canonical strand
/// (proteins have no reverse complement). Returns `(hash, pos)` per syncmer.
pub fn syncmer_protein<H: Hasher>(
    seq: &[u8],
    params: &SyncmerParams,
    hasher: H,
) -> anyhow::Result<Vec<(u64, usize)>> {
    params.validate()?;
    let k = params.smer;
    let w = params.window;
    let n = seq.len();
    if n < k {
        return Ok(Vec::new());
    }
    let syncs = closed_syncmers_stream(w, seq.windows(k).map(|smer| (hasher.hash(smer), ())));
    Ok(syncs.into_iter().map(|(pos, h, _)| (h, pos)).collect())
}

/// Build a syncmer hash set, dispatching on sequence type.
///
/// Protein uses `RapidHash` on s-mer bytes (no canonical); DNA uses the
/// 2-bit canonical rolling hash. Drop-in parallel to `libs::hash::seq_mins`.
pub fn seq_syncmer_set(
    seq: &[u8],
    params: &SyncmerParams,
    is_protein: bool,
) -> anyhow::Result<rapidhash::RapidHashSet<u64>> {
    if is_protein {
        let v = syncmer_protein(seq, params, RapidHash)?;
        Ok(v.into_iter().map(|(h, _)| h).collect())
    } else {
        let v = syncmer_dna(seq, params)?;
        Ok(v.into_iter().map(|(h, _, _)| h).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_core_min_in_middle_not_syncmer() {
        // hashes = [5,3,8,3,6], window=3
        // [5,3,8] min@1 (middle) -> no; [3,8,3] min@1 (start) -> yes; [8,3,6] min@3 (middle) -> no
        let hashes = vec![5, 3, 8, 3, 6];
        assert_eq!(closed_syncmers_from_hashes(&hashes, 3), vec![1]);
    }

    #[test]
    fn test_core_min_at_first_and_last() {
        // [1,5,5] min@0 (first) -> yes, argmin=0; [5,5,1] min@3 (last) -> yes, argmin=3
        let hashes = vec![1, 5, 5, 1];
        assert_eq!(closed_syncmers_from_hashes(&hashes, 3), vec![0, 3]);
    }

    #[test]
    fn test_core_too_short() {
        assert!(closed_syncmers_from_hashes(&[1, 2], 3).is_empty());
        assert!(closed_syncmers_from_hashes(&[], 3).is_empty());
        assert!(closed_syncmers_from_hashes(&[1, 2, 3], 0).is_empty());
    }

    #[test]
    fn test_core_window_one_emits_all() {
        // window=1: every position is both first and last -> all are syncmers
        let hashes = vec![7, 3, 9];
        assert_eq!(closed_syncmers_from_hashes(&hashes, 1), vec![0, 1, 2]);
    }

    #[test]
    fn test_dna_basic() {
        let params = SyncmerParams {
            smer: 4,
            window: 3,
            seed: 7,
        };
        let seq = b"ACGTACGTACGTACGT";
        let syncs = syncmer_dna(seq, &params).unwrap();
        assert!(!syncs.is_empty());
        for (_, pos, _) in &syncs {
            assert!(*pos <= seq.len() - params.smer);
        }
    }

    #[test]
    fn test_dna_too_short() {
        let params = SyncmerParams {
            smer: 8,
            window: 4,
            seed: 7,
        };
        // length 10 < smer+window-1 = 11
        assert!(syncmer_dna(b"ACGTACGTAC", &params).unwrap().is_empty());
    }

    #[test]
    fn test_dna_invalid_params() {
        let p = SyncmerParams {
            smer: 0,
            window: 3,
            seed: 7,
        };
        assert!(syncmer_dna(b"ACGT", &p).is_err());
        let p = SyncmerParams {
            smer: 32,
            window: 3,
            seed: 7,
        };
        assert!(syncmer_dna(b"ACGT", &p).is_err());
    }

    #[test]
    fn test_validate_rejects_extreme_window() {
        // An extreme window must be rejected (the pgi build sizes a ring
        // buffer from `window + 2`; capping prevents overflow / OOM).
        let p = SyncmerParams {
            smer: 8,
            window: usize::MAX,
            seed: 7,
        };
        assert!(p.validate().is_err());
        let p = SyncmerParams {
            smer: 8,
            window: 1_000_001,
            seed: 7,
        };
        assert!(p.validate().is_err());
        let p = SyncmerParams {
            smer: 8,
            window: 1_000_000,
            seed: 7,
        };
        assert!(p.validate().is_ok());
    }

    #[test]
    fn test_dna_canonical_revcomp() {
        // A sequence and its reverse complement yield the same syncmer hash set.
        let params = SyncmerParams {
            smer: 5,
            window: 4,
            seed: 7,
        };
        let seq = b"ACGTACGTACGTACGTGGCGCGCATATATACGTACGT";
        let rc = revcomp_dna(seq);
        let set1: std::collections::HashSet<u64> = syncmer_dna(seq, &params)
            .unwrap()
            .into_iter()
            .map(|(h, _, _)| h)
            .collect();
        let set2: std::collections::HashSet<u64> = syncmer_dna(&rc, &params)
            .unwrap()
            .into_iter()
            .map(|(h, _, _)| h)
            .collect();
        assert_eq!(set1, set2, "DNA syncmer set must be strand-symmetric");
    }

    #[test]
    fn test_dna_density_reasonable() {
        // Closed syncmers have average depth ~2/(window+1). Verify the syncmer
        // count is in a reasonable band around that expectation. Closed
        // syncmers do NOT guarantee every position is covered (sequence ends
        // may be uncovered; syng patches them with X/Y ends), so we check
        // density rather than full cover.
        let params = SyncmerParams {
            smer: 6,
            window: 8,
            seed: 7,
        };
        let seq = b"ACGTACGTACGTGGGGCCCCACGTACGTACGTGGGGCCCCttttacgtACGTACGTGGGGCCCC";
        let syncs = syncmer_dna(seq, &params).unwrap();
        assert!(!syncs.is_empty());
        let n_smers = seq.len() - params.smer + 1;
        let expected = 2.0 * n_smers as f64 / (params.window + 1) as f64;
        let actual = syncs.len() as f64;
        assert!(
            actual > expected * 0.3 && actual < expected * 3.0,
            "syncmer density {} far from expected ~{}",
            actual,
            expected
        );
    }

    #[test]
    fn test_protein_basic() {
        let params = SyncmerParams {
            smer: 3,
            window: 5,
            seed: 7,
        };
        let seq = b"ACDEFGHIKLMNPQRSTVWY"; // 20 aa
        let syncs = syncmer_protein(seq, &params, RapidHash).unwrap();
        assert!(!syncs.is_empty());
        for (_, pos) in &syncs {
            assert!(*pos <= seq.len() - params.smer);
        }
    }

    #[test]
    fn test_dispatch_dna_and_protein() {
        let dna_params = SyncmerParams::default_dna();
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let set = seq_syncmer_set(seq, &dna_params, false).unwrap();
        assert!(!set.is_empty());

        let prot_params = SyncmerParams {
            smer: 3,
            window: 5,
            seed: 7,
        };
        let pseq = b"ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY";
        let pset = seq_syncmer_set(pseq, &prot_params, true).unwrap();
        assert!(!pset.is_empty());
    }

    fn revcomp_dna(seq: &[u8]) -> Vec<u8> {
        crate::libs::nt::rev_comp(seq).collect()
    }

    /// Closed syncmers guarantee a bounded gap between consecutive syncmer
    /// endpoints: for `w >= 2`, consecutive endpoint positions differ by at
    /// most `2*(w-1)` (tight bound, verified exhaustively over all-distinct
    /// permutations and 4-letter alphabets for w in {3,4,5}, and randomly for
    /// w up to 32). This bounded-gap property (no large unsampled region) is
    /// the coverage guarantee minimizers lack. Sequence ends are NOT
    /// guaranteed (only interior gaps). Proof sketch: take the minimum-hash
    /// position `p*` of a non-endpoint run; since `p*` is not an endpoint,
    /// each of its two w-windows must hold a smaller value, and those smaller
    /// values must lie outside the run (p* is the run min), forcing the run
    /// endpoints `a < b` to satisfy `b-w+1 <= p* <= a+w-1`, i.e. `b-a <=
    /// 2*(w-1)`.
    #[test]
    fn test_core_bounded_gap_property() {
        // Deterministic LCG for reproducible random hash streams.
        fn lcg(state: &mut u64) -> u64 {
            *state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *state
        }

        for &w in &[2usize, 3, 4, 5, 8, 16] {
            for &n in &[w, w + 1, w + 5, 20, 50, 100, 200, 500] {
                for seed in 0..8u64 {
                    let mut state = 0xA5A5_5A5A_A5A5_5A5A_u64
                        .wrapping_mul(w as u64)
                        .wrapping_add(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15));
                    let hashes: Vec<u64> = (0..n).map(|_| lcg(&mut state)).collect();
                    let positions = closed_syncmers_from_hashes(&hashes, w);
                    let unique: Vec<usize> = positions
                        .into_iter()
                        .collect::<std::collections::BTreeSet<_>>()
                        .into_iter()
                        .collect();
                    for pair in unique.windows(2) {
                        assert!(
                            pair[1] - pair[0] <= 2 * w - 2,
                            "gap {} exceeds 2(w-1)={} (w={}, n={}, seed={})",
                            pair[1] - pair[0],
                            2 * w - 2,
                            w,
                            n,
                            seed
                        );
                    }
                }
            }
        }
    }
}
