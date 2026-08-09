//! SIMD-accelerated per-byte nucleotide statistics.
//!
//! Platform policy (same as `libs::hv` and `libs::poa::simd`): a hand-written
//! AVX2 path is the primary x86-64 implementation, dispatched at runtime with
//! `is_x86_feature_detected!`; other targets fall through to portable
//! implementations. All paths are bit-identical to the scalar reference.

use wide::u8x32;

/// SIMD implementation selector, mainly for benchmarking and diagnostics.
pub enum SimdPath {
    /// Runtime detection: AVX2 on capable x86-64, portable elsewhere.
    Auto,
    /// Force the portable `wide`/scalar path.
    Wide,
    /// Force the AVX2 path (falls back to portable on CPUs without AVX2).
    Avx2,
}

/// Lowercased ASCII values of A/C/G/T/U (the `NT_VAL` 0..=3 set).
const VALID_LOWER: [u8; 5] = [0x61, 0x63, 0x67, 0x74, 0x75];

/// Lowercased ASCII values mapped to `Nt::N` (IUPAC ambiguous + N + X).
const N_LOWER: [u8; 12] = [
    0x62, 0x64, 0x68, 0x6B, 0x6D, 0x6E, 0x72, 0x73, 0x76, 0x77, 0x78, 0x79,
];

/// True if `b` maps to one of A/C/G/T/U (case-insensitive), matching
/// `nt::to_nt(b) != Nt::N | Nt::Invalid`.
#[inline]
fn is_valid_base(b: u8) -> bool {
    matches!(b | 0x20, 0x61 | 0x63 | 0x67 | 0x74 | 0x75)
}

/// True if `b` maps to `Nt::N` (IUPAC ambiguous codes + N + X), matching
/// `nt::NT_VAL[b] == Nt::N`.
#[inline]
fn is_n_base(b: u8) -> bool {
    matches!(
        b | 0x20,
        0x62 | 0x64 | 0x68 | 0x6B | 0x6D | 0x6E | 0x72 | 0x73 | 0x76 | 0x77 | 0x78 | 0x79
    )
}

/// True if `b` is a lowercase ASCII letter.
#[inline]
fn is_lower_base(b: u8) -> bool {
    b.is_ascii_lowercase()
}

/// Counts A/C/G/T/U bases (case-insensitive); equivalent to
/// `seq.iter().filter(|b| !matches!(nt::to_nt(b), Nt::N | Nt::Invalid)).count()`.
pub fn count_valid(seq: &[u8]) -> usize {
    count_valid_with(SimdPath::Auto, seq)
}

/// [`count_valid`] with an explicit implementation path.
pub fn count_valid_with(path: SimdPath, seq: &[u8]) -> usize {
    match path {
        SimdPath::Wide => count_valid_wide(seq),
        SimdPath::Auto | SimdPath::Avx2 => {
            #[cfg(target_arch = "x86_64")]
            if is_x86_feature_detected!("avx2") {
                // SAFETY: gated on runtime AVX2 support.
                return unsafe { avx2::count_valid_avx2(seq) };
            }
            count_valid_wide(seq)
        }
    }
}

/// Counts bases mapped to `Nt::N` (IUPAC ambiguous codes + N + X);
/// equivalent to `nt::count_n`.
pub fn count_n(seq: &[u8]) -> usize {
    count_n_with(SimdPath::Auto, seq)
}

/// [`count_n`] with an explicit implementation path.
pub fn count_n_with(path: SimdPath, seq: &[u8]) -> usize {
    let scalar = |s: &[u8]| s.iter().filter(|&&b| is_n_base(b)).count();
    match path {
        SimdPath::Wide => scalar(seq),
        SimdPath::Auto | SimdPath::Avx2 => {
            #[cfg(target_arch = "x86_64")]
            if is_x86_feature_detected!("avx2") {
                // SAFETY: gated on runtime AVX2 support.
                return unsafe { avx2::count_n_avx2(seq) };
            }
            scalar(seq)
        }
    }
}

/// Builds a per-base mask bitmap: word `w` bit `k` is set iff `seq[w*32+k]`
/// is masked (`gap_only`: N-family; otherwise lowercase letters or N-family).
pub fn masked_bitmap(seq: &[u8], gap_only: bool) -> Vec<u32> {
    masked_bitmap_with(SimdPath::Auto, seq, gap_only)
}

/// [`masked_bitmap`] with an explicit implementation path.
pub fn masked_bitmap_with(path: SimdPath, seq: &[u8], gap_only: bool) -> Vec<u32> {
    let scalar = |s: &[u8]| {
        s.chunks(32)
            .map(|chunk| {
                chunk.iter().enumerate().fold(0u32, |w, (k, &b)| {
                    let masked = if gap_only {
                        is_n_base(b)
                    } else {
                        is_n_base(b) || is_lower_base(b)
                    };
                    w | ((masked as u32) << k)
                })
            })
            .collect()
    };
    match path {
        SimdPath::Wide => scalar(seq),
        SimdPath::Auto | SimdPath::Avx2 => {
            #[cfg(target_arch = "x86_64")]
            if is_x86_feature_detected!("avx2") {
                // SAFETY: gated on runtime AVX2 support.
                return unsafe { avx2::masked_bitmap_avx2(seq, gap_only) };
            }
            scalar(seq)
        }
    }
}

/// Portable `wide` path for [`count_valid`]: 32-byte lanes, equality via
/// saturating-subtraction masks (bit 7 set when unequal), popcount per word.
fn count_valid_wide(seq: &[u8]) -> usize {
    let or20 = u8x32::splat(0x20);
    let mut count = 0usize;
    let (chunks, remainder) = seq.as_chunks::<32>();
    for chunk in chunks {
        let t = u8x32::from(*chunk) | or20;
        let mut bits = 0u32;
        for &val in &VALID_LOWER {
            bits |= eq_bits(t, u8x32::splat(val));
        }
        count += bits.count_ones() as usize;
    }
    count + remainder.iter().filter(|&&b| is_valid_base(b)).count()
}

/// Lane-wise `a == x` as a 32-bit mask (bit set per matching lane).
#[inline]
fn eq_bits(a: u8x32, x: u8x32) -> u32 {
    let d = a.saturating_sub(x) | x.saturating_sub(a);
    // d != 0  <=>  bit 7 of (d | -d) is set; -d wraps as (!d + 1).
    let ne = d | (!d + u8x32::splat(1));
    !ne.to_bitmask()
}

#[cfg(target_arch = "x86_64")]
mod avx2 {
    use super::*;
    use std::arch::x86_64::*;

    #[inline]
    unsafe fn load(ptr: *const u8) -> __m256i {
        _mm256_loadu_si256(ptr as *const __m256i)
    }

    #[inline]
    unsafe fn set1(b: u8) -> __m256i {
        _mm256_set1_epi8(b as i8)
    }

    #[inline]
    unsafe fn popcount(m: __m256i) -> usize {
        _mm256_movemask_epi8(m).count_ones() as usize
    }

    /// `(v | 0x20)` membership in `vals` as a byte mask.
    #[inline]
    unsafe fn eq_any_lower(v: __m256i, vals: &[u8]) -> __m256i {
        let or20 = set1(0x20);
        let t = _mm256_or_si256(v, or20);
        let mut m = _mm256_setzero_si256();
        for &val in vals {
            m = _mm256_or_si256(m, _mm256_cmpeq_epi8(t, set1(val)));
        }
        m
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn count_valid_avx2(seq: &[u8]) -> usize {
        let mut count = 0usize;
        let (chunks, remainder) = seq.as_chunks::<32>();
        for chunk in chunks {
            let m = eq_any_lower(load(chunk.as_ptr()), &VALID_LOWER);
            count += popcount(m);
        }
        count + remainder.iter().filter(|&&b| is_valid_base(b)).count()
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn count_n_avx2(seq: &[u8]) -> usize {
        let mut count = 0usize;
        let (chunks, remainder) = seq.as_chunks::<32>();
        for chunk in chunks {
            let m = eq_any_lower(load(chunk.as_ptr()), &N_LOWER);
            count += popcount(m);
        }
        count + remainder.iter().filter(|&&b| is_n_base(b)).count()
    }

    /// One 32-byte word of the mask bitmap.
    #[inline]
    unsafe fn mask_word(v: __m256i, gap_only: bool) -> u32 {
        if gap_only {
            _mm256_movemask_epi8(eq_any_lower(v, &N_LOWER)) as u32
        } else {
            // Lowercase letters, or uppercase letters whose lowercased form is
            // in the N-family (IUPAC ambiguous codes + N + X).
            let lower = _mm256_cmpeq_epi8(_mm256_max_epu8(v, set1(0x61)), v);
            let upper = _mm256_cmpeq_epi8(_mm256_min_epu8(v, set1(0x5A)), v);
            let t = _mm256_or_si256(v, set1(0x20));
            let mut n_lower = _mm256_setzero_si256();
            for &val in &N_LOWER {
                n_lower = _mm256_or_si256(n_lower, _mm256_cmpeq_epi8(t, set1(val)));
            }
            let iupac_upper = _mm256_and_si256(upper, n_lower);
            _mm256_movemask_epi8(_mm256_or_si256(lower, iupac_upper)) as u32
        }
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn masked_bitmap_avx2(seq: &[u8], gap_only: bool) -> Vec<u32> {
        let n_words = seq.len().div_ceil(32);
        let mut out = vec![0u32; n_words];
        let (chunks, remainder) = seq.as_chunks::<32>();
        for (wi, chunk) in chunks.iter().enumerate() {
            out[wi] = mask_word(load(chunk.as_ptr()), gap_only);
        }
        let start = chunks.len();
        for (k, &b) in remainder.iter().enumerate() {
            let masked = if gap_only {
                is_n_base(b)
            } else {
                is_n_base(b) || is_lower_base(b)
            };
            if masked {
                out[start] |= 1 << k;
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::StdRng, Rng, SeedableRng};

    fn random_seq(rng: &mut StdRng, len: usize) -> Vec<u8> {
        let alphabet = *b"ACGTacgtNnMRWSYKVHDXBXmrwsykvhdbx0123456789-*._ ZzUu";
        (0..len)
            .map(|_| alphabet[rng.random_range(0..alphabet.len())])
            .collect()
    }

    fn scalar_masked(seq: &[u8], gap_only: bool) -> Vec<u32> {
        seq.chunks(32)
            .map(|chunk| {
                chunk.iter().enumerate().fold(0u32, |w, (k, &b)| {
                    let masked = if gap_only {
                        is_n_base(b)
                    } else {
                        is_n_base(b) || is_lower_base(b)
                    };
                    w | ((masked as u32) << k)
                })
            })
            .collect()
    }

    #[test]
    fn counts_match_scalar_random() {
        let mut rng = StdRng::seed_from_u64(2026);
        for len in [0usize, 1, 31, 32, 33, 63, 64, 65, 1000, 10000] {
            let seq = random_seq(&mut rng, len);
            let expected_valid = seq.iter().filter(|&&b| is_valid_base(b)).count();
            let expected_n = seq.iter().filter(|&&b| is_n_base(b)).count();
            assert_eq!(count_valid(&seq), expected_valid, "valid len={len}");
            assert_eq!(count_n(&seq), expected_n, "n len={len}");
        }
    }

    #[test]
    fn counts_match_nt_reference() {
        // Cross-check against the existing nt::NT_VAL semantics.
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..20 {
            let seq = random_seq(&mut rng, 500);
            let n_ref = crate::libs::nt::count_n(&seq);
            assert_eq!(count_n(&seq), n_ref);
            let valid_ref = seq
                .iter()
                .filter(|&&b| {
                    !matches!(
                        crate::libs::nt::to_nt(b),
                        crate::libs::nt::Nt::N | crate::libs::nt::Nt::Invalid
                    )
                })
                .count();
            assert_eq!(count_valid(&seq), valid_ref);
        }
    }

    #[test]
    fn bitmap_matches_scalar() {
        let mut rng = StdRng::seed_from_u64(99);
        for len in [0usize, 1, 31, 32, 33, 64, 65, 1000] {
            for gap_only in [false, true] {
                let seq = random_seq(&mut rng, len);
                assert_eq!(
                    masked_bitmap(&seq, gap_only),
                    scalar_masked(&seq, gap_only),
                    "bitmap len={len} gap={gap_only}"
                );
            }
        }
    }

    #[test]
    fn find_masked_regions_matches_scalar() {
        let mut rng = StdRng::seed_from_u64(123);
        for len in [0usize, 1, 31, 32, 33, 64, 65, 500] {
            for gap_only in [false, true] {
                let seq = random_seq(&mut rng, len);
                let bitmap = masked_bitmap(&seq, gap_only);
                let regions = super::super::fmt::fa::find_masked_regions(&seq, gap_only);
                // Reconstruct expected regions from the scalar mask.
                let mut expected = Vec::new();
                let mut begin: Option<usize> = None;
                let mut end: Option<usize> = None;
                for (i, &b) in seq.iter().enumerate() {
                    let masked = if gap_only {
                        is_n_base(b)
                    } else {
                        is_n_base(b) || is_lower_base(b)
                    };
                    if masked {
                        if begin.is_none() {
                            begin = Some(i);
                        }
                        end = Some(i);
                    } else if let (Some(b), Some(e)) = (begin, end) {
                        expected.push((b, e));
                        begin = None;
                        end = None;
                    }
                }
                if let (Some(b), Some(e)) = (begin, end) {
                    expected.push((b, e));
                }
                assert_eq!(regions, expected, "regions len={len} gap={gap_only}");
                let _ = bitmap;
            }
        }
    }
}
