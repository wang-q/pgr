use rand::{RngCore, SeedableRng};
use rapidhash::RapidRng;
use wide::{bytemuck, i32x8, u16x8, u32x8, u8x16};

#[cfg(target_arch = "x86_64")]
mod rng_jump {
    // RapidRng's state is a counter advanced by a constant step, so output j
    // of a seed (1-based) is mix(seed + j*SECRET0, ...): cheap random access
    // enables a block-major loop that keeps a chunk of the HV in registers
    // while sweeping all seeds over it, and lets independent seed streams
    // overlap in the CPU.
    pub(super) const RAPID_SECRET0: u64 = 0x2d358dccaa6c78a5;
    pub(super) const RAPID_SECRET1: u64 = 0x8bb84b93962eacc9;

    #[inline(always)]
    pub(super) fn rapid_mix(a: u64, b: u64) -> u64 {
        let r = (a as u128) * (b as u128);
        (r as u64) ^ ((r >> 64) as u64)
    }

    #[inline(always)]
    pub(super) fn rnd_at(seed: u64, j: u64) -> u64 {
        let s = seed.wrapping_add(j.wrapping_mul(RAPID_SECRET0));
        rapid_mix(s, s ^ RAPID_SECRET1)
    }
}
#[cfg(target_arch = "x86_64")]
use rng_jump::rnd_at;

/// AVX2 bit encoding, block-major. Bit-identical to the scalar `hash_hv_bit`
/// (each 64-dim chunk consumes all 64 bits of one u64 of the RapidRng stream:
/// low 32 bits → dims 0..32, high 32 bits → dims 32..64; ±1 values accumulated
/// in eight 8-lane registers per chunk). The −N offset is deferred: per seed
/// only `2·bit` is accumulated (values balance around 0), and N is subtracted
/// once per chunk when storing — one vector op less per seed per group.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn hash_hv_bit_avx2(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    use std::arch::x86_64::*;
    let n = seed_vec.len() as i32;
    let mut hv = vec![0i32; hv_d];
    let num_chunk = hv_d / 64;
    let one = _mm256_set1_epi32(1);
    let s0 = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    let s1 = _mm256_set_epi32(15, 14, 13, 12, 11, 10, 9, 8);
    let s2 = _mm256_set_epi32(23, 22, 21, 20, 19, 18, 17, 16);
    let s3 = _mm256_set_epi32(31, 30, 29, 28, 27, 26, 25, 24);
    let nv = _mm256_set1_epi32(n);
    for b in 0..num_chunk {
        let j = b as u64 + 1;
        let mut a = [_mm256_setzero_si256(); 8];
        for &seed in seed_vec {
            let r = rnd_at(seed, j);
            let vlo = _mm256_set1_epi32((r as u32) as i32);
            let vhi = _mm256_set1_epi32(((r >> 32) as u32) as i32);
            for k in 0..4 {
                let shift = [s0, s1, s2, s3][k];
                let bl = _mm256_and_si256(_mm256_srlv_epi32(vlo, shift), one);
                a[k] = _mm256_add_epi32(a[k], _mm256_slli_epi32(bl, 1));
                let bh = _mm256_and_si256(_mm256_srlv_epi32(vhi, shift), one);
                a[k + 4] = _mm256_add_epi32(a[k + 4], _mm256_slli_epi32(bh, 1));
            }
        }
        let base = b * 64;
        for (k, acc) in a.iter().enumerate() {
            _mm256_storeu_si256(
                hv[base + k * 8..base + (k + 1) * 8].as_mut_ptr() as *mut _,
                _mm256_sub_epi32(*acc, nv),
            );
        }
    }
    // Tail dims (hv_d not a multiple of 64) keep the −N offset, matching the
    // portable fallback exactly.
    for v in &mut hv[num_chunk * 64..] {
        *v = -n;
    }
    hv
}

/// AVX2 i8 encoding, block-major. Bit-identical to the scalar `hash_hv_i8`:
/// one u64 (8 dims) is sign-extended to 32-bit lanes with a single
/// `vpmovsxbd` per seed.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn hash_hv_i8_avx2(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    use std::arch::x86_64::*;
    let mut hv = vec![0i32; hv_d];
    let num_chunk = hv_d / 8;
    for b in 0..num_chunk {
        let j = b as u64 + 1;
        let mut acc = _mm256_setzero_si256();
        for &seed in seed_vec {
            let bytes = rnd_at(seed, j).to_ne_bytes();
            let v = _mm256_cvtepi8_epi32(_mm_loadl_epi64(bytes.as_ptr() as *const __m128i));
            acc = _mm256_add_epi32(acc, v);
        }
        _mm256_storeu_si256(hv[b * 8..b * 8 + 8].as_mut_ptr() as *mut _, acc);
    }
    hv
}

/// Generates a hypervector (HV) from a set of k-mer hash values using a SIMD-optimized implementation.
///
/// # Arguments
/// * `kmer_hash_set` - A set of k-mer hash values.
/// * `hv_d` - The dimension of the hypervector.
///
/// # Returns
/// A hypervector of dimension `hv_d` represented as a `Vec<i32>`.
///
/// # Formula
/// The hypervector is generated as:
/// \[
/// \mathbf{H} = \sum_{i=1}^{N} (hv^{i} \times 2 - 1)
/// \]
/// where \(N\) is the number of k-mer hash values, and \(hv^{i}\) is a binary hypervector derived from the k-mer hash.
///
/// # Notes
/// This function uses SIMD instructions to process 8 bits at a time, improving performance over the serial implementation.
pub fn hash_hv_bit(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    // Platform policy: AVX2 (256-bit) is the primary x86-64 path; other
    // targets/CPUs (aarch64 NEON, scalar, ...) fall through to the portable
    // wide implementation below. All paths are bit-identical.
    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx2") {
        return unsafe { hash_hv_bit_avx2(seed_vec, hv_d) };
    }
    let num_seed = seed_vec.len();
    let num_chunk = hv_d / 64;
    let mut hv = vec![-(num_seed as i32); hv_d];

    // Loop through all seeds
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);

        // Consume all 64 bits of each u64: low 32 bits → dims 0..32,
        // high 32 bits → dims 32..64 (half the RNG calls of the 32-bit path).
        for i in 0..num_chunk {
            let rnd_bits = rng.next_u64();
            let halves = [(rnd_bits as u32), (rnd_bits >> 32) as u32];

            for (k, half) in halves.iter().enumerate() {
                // Use SIMD to process 8 bits at a time
                for j in (0..32).step_by(8) {
                    let bit_mask = u32x8::splat(1);
                    let shift = u32x8::from([
                        j as u32,
                        (j + 1) as u32,
                        (j + 2) as u32,
                        (j + 3) as u32,
                        (j + 4) as u32,
                        (j + 5) as u32,
                        (j + 6) as u32,
                        (j + 7) as u32,
                    ]);
                    let bits = (u32x8::splat(*half) >> shift) & bit_mask;

                    // Convert bits to i32 and shift left by 1 (0/1 bit pattern
                    // is identical in u32 and i32, so a reinterpret cast suffices)
                    let bits_i32: i32x8 = bytemuck::cast(bits);
                    let bits_i32 = bits_i32 << i32x8::splat(1);

                    let base = i * 64 + k * 32 + j;
                    let mut hv_simd = i32x8::from(&hv[base..base + 8]);
                    hv_simd += bits_i32;
                    hv[base..base + 8].copy_from_slice(&hv_simd.to_array());
                }
            }
        }
    }

    hv
}

/// Generates a hypervector (HV) using i8 values (-128..=127) as the basic unit.
///
/// # Arguments
/// * `seed_vec` - A set of k-mer hash values (seeds).
/// * `hv_d` - The dimension of the hypervector.
///
/// # Returns
/// A hypervector of dimension `hv_d` represented as a `Vec<i32>`.
///
/// # Notes
/// This implementation avoids bit manipulation overhead by using `i8` directly,
/// but requires more RNG calls (1 u64 per 8 dimensions) compared to the bit-based approach.
/// It uses SIMD to process 8 dimensions at a time.
pub fn hash_hv_i8(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    // Platform policy: same AVX2 dispatch as `hash_hv_bit`, portable wide
    // fallback elsewhere; results are bit-identical.
    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx2") {
        return unsafe { hash_hv_i8_avx2(seed_vec, hv_d) };
    }
    // Initialize HV with 0.
    // We accumulate random i8 values (-128..=127) directly.
    let mut hv = vec![0i32; hv_d];

    // Loop through all seeds
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);

        // Process 8 dimensions per chunk (1 u64 = 8 bytes)
        let num_chunk = hv_d / 8;

        for i in 0..num_chunk {
            let rnd_bits = rng.next_u64();
            let bytes = rnd_bits.to_ne_bytes();

            // Sign-extend each byte as i8, then to i32: (b << 24) >> 24
            // with an arithmetic right shift equals b as i8 as i32.
            let mut arr = [0u8; 16];
            arr[..8].copy_from_slice(&bytes);
            let vec_u8 = u8x16::from(arr);
            let vec_u16 = u16x8::from_u8x16_low(vec_u8);
            let vec_i32 = i32x8::from_u16x8(vec_u16);
            let vec_vals = (vec_i32 << i32x8::splat(24)) >> i32x8::splat(24);

            // Load current HV values
            let mut hv_simd = i32x8::from(&hv[i * 8..(i + 1) * 8]);

            // Accumulate
            hv_simd += vec_vals;

            // Store back
            hv[i * 8..(i + 1) * 8].copy_from_slice(&hv_simd.to_array());
        }
    }

    hv
}

/// Splitmix64 step (deterministic, well-distributed for sparse projection).
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

/// Sparse hypervector projection: each seed updates `s` random dimensions
/// with ±1 (locality-sensitive hashing style).
///
/// Dense encodings saturate for large seed sets (every dimension accumulates
/// contributions from all seeds); the sparse projection keeps the shared-seed
/// signal dominant, so cosine similarity on the result approximates the k-mer
/// set overlap (see notes/benchmarks/dist-cohort-validation.md).
pub fn hash_hv_sparse(seed_vec: &[u64], hv_d: usize, s: usize) -> Vec<i32> {
    let mut hv = vec![0i32; hv_d];
    for &seed in seed_vec {
        let mut x = splitmix64(seed);
        for _ in 0..s {
            x = splitmix64(x);
            let idx = (x % hv_d as u64) as usize;
            if ((x >> 32) & 1) == 1 {
                hv[idx] += 1;
            } else {
                hv[idx] -= 1;
            }
        }
    }
    hv
}

#[allow(dead_code)]
fn hash_hv_i8_serial(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    // Initialize HV with 0.
    let mut hv = vec![0i32; hv_d];

    // Loop through all seeds
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);

        // Process dimensions in chunks of 8 (matching SIMD behavior for RNG alignment)
        let num_chunk = hv_d / 8;

        for i in 0..num_chunk {
            let rnd_bits = rng.next_u64();
            let bytes = rnd_bits.to_ne_bytes();

            // Iterate over each byte in the u64
            for j in 0..8 {
                let val_i8 = bytes[j] as i8;
                hv[i * 8 + j] += val_i8 as i32;
            }
        }
    }

    hv
}

#[allow(dead_code)]
fn hv_norm_l2_sq_serial(hv: &[i32]) -> f32 {
    let norm_sq = hv
        .iter()
        .fold(0., |sum: f32, &num| sum + (num as f32 * num as f32));
    norm_sq
}

/// Computes the squared L2 norm of a hypervector using a SIMD-optimized implementation.
///
/// # Arguments
/// * `hv` - The hypervector represented as a slice of `i32`.
///
/// # Returns
/// The squared L2 norm of the hypervector as an `f32`.
pub fn hv_norm_l2_sq(hv: &[i32]) -> f32 {
    let a_f32: Vec<f32> = hv.iter().map(|&x| x as f32).collect();
    crate::libs::linalg::norm_l2_sq(&a_f32)
}

/// Computes the cardinality of a set represented by a hypervector.
///
/// # Arguments
/// * `hv` - The hypervector represented as a slice of `i32`.
/// * `hv_d` - The dimension of the hypervector.
///
/// # Returns
/// The cardinality of the set as a `usize`.
///
/// # Formula
/// The cardinality is computed as:
/// \[
/// |\mathcal{S}_k(A)| = \frac{\|\mathbf{H}_A\|_2^2}{D}
/// \]
/// where \(\|\mathbf{H}_A\|_2^2\) is the squared L2 norm of the hypervector, and \(D\) is the dimension of the hypervector.
pub fn hv_cardinality(hv: &[i32]) -> usize {
    let norm_sq = hv_norm_l2_sq(hv);
    (norm_sq / hv.len() as f32) as usize
}

/// Computes the dot product of two hypervectors.
///
/// # Arguments
/// * `a` - The first hypervector represented as a slice of `i32`.
/// * `b` - The second hypervector represented as a slice of `i32`.
///
/// # Returns
/// The dot product of the two hypervectors as an `f32`.
pub fn hv_dot(a: &[i32], b: &[i32]) -> f32 {
    let hv_d_sqrt = (a.len() as f32).sqrt();
    let a_f32: Vec<_> = a.iter().map(|&x| x as f32 / hv_d_sqrt).collect();
    let b_f32: Vec<_> = b.iter().map(|&x| x as f32 / hv_d_sqrt).collect();

    crate::libs::linalg::dot_product(&a_f32, &b_f32)
}

/// A hypervector entry with its source name and the resulting HV set.
#[derive(Debug, Default, Clone)]
pub struct HvEntry {
    pub name: String,
    pub set: Vec<i32>,
}

/// Pairwise distance metrics between two hypervectors.
#[derive(Debug, Clone)]
pub struct HvDistances {
    pub card1: usize,
    pub card2: usize,
    pub inter: usize,
    pub union: usize,
    pub mash: f32,
    pub jaccard: f32,
    pub containment: f32,
}

/// Calculate Jaccard, Containment, and Mash distance between two hypervector sets.
pub fn calc_distances(s1: &[i32], s2: &[i32], kmer: usize) -> HvDistances {
    let card1 = hv_cardinality(s1);
    let card2 = hv_cardinality(s2);

    let inter = hv_dot(s1, s2).min(card1 as f32).min(card2 as f32);
    let union = card1 as f32 + card2 as f32 - inter;

    let jaccard = inter / union;
    let containment = inter / card1 as f32;
    let mash = crate::libs::hash::mash_distance(jaccard as f64, kmer) as f32;

    HvDistances {
        card1,
        card2,
        inter: inter as usize,
        union: union as usize,
        mash,
        jaccard,
        containment,
    }
}

/// Load a single FASTA file into one `HvEntry` by merging all sequences' minimizers.
pub fn load_hv_from_fasta(
    infile: &str,
    hasher: &str,
    kmer: usize,
    window: usize,
    dim: usize,
) -> anyhow::Result<HvEntry> {
    let mut fa_in = crate::libs::fmt::fa::reader(infile)?;

    let mut file_set = rapidhash::RapidHashSet::default();

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;
        let seq = record.sequence();

        let set: rapidhash::RapidHashSet<u64> =
            crate::libs::hash::seq_mins(&seq[..], hasher, kmer, window)?;

        file_set.extend(set);
    }

    let seed_vec: Vec<u64> = file_set.into_iter().collect();
    let hv: Vec<i32> = hash_hv_i8(&seed_vec, dim);
    let entry = HvEntry {
        name: infile.to_string(),
        set: hv,
    };

    Ok(entry)
}

/// Load a single FASTA file into one `HvEntry` by merging all sequences' syncmers.
///
/// Drop-in parallel to [`load_hv_from_fasta`]; `smer` is the s-mer length and
/// `window` the number of s-mers per syncmer window. `is_protein` dispatches
/// to the protein byte-hash path (DNA uses the 2-bit canonical rolling hash).
pub fn load_hv_from_fasta_syncmer(
    infile: &str,
    smer: usize,
    window: usize,
    is_protein: bool,
    dim: usize,
) -> anyhow::Result<HvEntry> {
    let params = crate::libs::syncmer::SyncmerParams {
        smer,
        window,
        seed: 7,
    };
    params.validate()?;

    let mut fa_in = crate::libs::fmt::fa::reader(infile)?;

    let mut file_set = rapidhash::RapidHashSet::default();

    for result in fa_in.records() {
        let record = result?;
        let seq = record.sequence();
        let set = crate::libs::syncmer::seq_syncmer_set(&seq[..], &params, is_protein)?;
        file_set.extend(set);
    }

    let seed_vec: Vec<u64> = file_set.into_iter().collect();
    let hv: Vec<i32> = hash_hv_i8(&seed_vec, dim);
    Ok(HvEntry {
        name: infile.to_string(),
        set: hv,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use rapidhash::RapidHashSet;

    /// Serial ±1 reference for `hash_hv_bit` (AVX-512 path parity check).
    fn hash_hv_bit_serial(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
        let num_seed = seed_vec.len();
        let mut hv = vec![-(num_seed as i32); hv_d];
        for hash in seed_vec {
            let mut rng = RapidRng::seed_from_u64(*hash);
            for i in 0..(hv_d / 64) {
                let rnd_bits = rng.next_u64();
                let halves = [(rnd_bits as u32), (rnd_bits >> 32) as u32];
                for (k, half) in halves.iter().enumerate() {
                    for j in 0..32 {
                        let idx = i * 64 + k * 32 + j;
                        hv[idx] += (((half >> j) & 1) << 1) as i32;
                    }
                }
            }
        }
        hv
    }

    #[test]
    fn test_hash_hv() {
        // Generate random input data
        let mut rng = rand::rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..1000).map(|_| rng.random::<u64>()).collect();
        let seed_vec: Vec<u64> = kmer_hash_set.into_iter().collect();
        let hv_d = 4096;

        // Run the SIMD version
        let hv = hash_hv_bit(&seed_vec, hv_d);

        // Check the dimension of the hypervector
        assert_eq!(hv.len(), hv_d, "Hypervector dimension mismatch!");

        // Check that the hypervector is not all zeros
        assert!(
            hv.iter().any(|&x| x != 0),
            "Hypervector should not be all zeros!"
        );
    }

    #[test]
    fn test_hash_hv_i8() {
        // Generate random input data
        let mut rng = rand::rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..1000).map(|_| rng.random::<u64>()).collect();
        let seed_vec: Vec<u64> = kmer_hash_set.into_iter().collect();
        let hv_d = 4096;

        // Run the i8 SIMD version
        let hv = hash_hv_i8(&seed_vec, hv_d);

        // Check the dimension of the hypervector
        assert_eq!(hv.len(), hv_d, "Hypervector dimension mismatch!");

        // Check that the hypervector is not all zeros
        assert!(
            hv.iter().any(|&x| x != 0),
            "Hypervector should not be all zeros!"
        );
    }

    #[test]
    fn test_hash_hv_i8_jaccard_dc_bias() {
        // Regression for the FASTA `dist hv` dimension mismatch (hv.md §3.4):
        // i8 bytes have mean -0.5, so each dimension accumulates a DC bias
        // (~ -N/2) that dominates the dot product and inflates the Jaccard
        // estimate towards a set-size-dependent constant instead of the true
        // shared-seed fraction.
        let n = 3000usize;
        let shared = 500usize;
        let hv_d = 4096;

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let pool: Vec<u64> = (0..(2 * n - shared)).map(|_| rng.random()).collect();

        let set_a: Vec<u64> = pool[..n].to_vec();
        let set_b: Vec<u64> = pool[..shared]
            .iter()
            .chain(pool[n..].iter())
            .copied()
            .collect();
        assert_eq!(set_b.len(), n);

        let d = calc_distances(&hash_hv_i8(&set_a, hv_d), &hash_hv_i8(&set_b, hv_d), 21);

        let true_j = shared as f32 / (2 * n - shared) as f32; // 500/5500 ≈ 0.0909
                                                              // The DC bias lifts the estimate far above the truth (~0.154, hv.md §3.4).
        assert!(
            (d.jaccard - 0.154).abs() < 0.02,
            "i8 Jaccard estimate {} should reproduce the DC-bias inflation (truth {})",
            d.jaccard,
            true_j
        );
    }

    #[test]
    fn test_hash_hv_i8_serial_vs_simd() {
        // Generate random input data
        let mut rng = rand::rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..1000).map(|_| rng.random::<u64>()).collect();
        let seed_vec: Vec<u64> = kmer_hash_set.into_iter().collect();
        let hv_d = 4096;

        // Run normal version
        let result_serial = hash_hv_i8_serial(&seed_vec, hv_d);

        // Run SIMD version
        let result_simd = hash_hv_i8(&seed_vec, hv_d);

        // Compare results
        assert_eq!(
            result_serial, result_simd,
            "SIMD version does not match serial version for i8 implementation!"
        );
    }

    #[test]
    fn test_hash_hv_bit_serial_vs_simd() {
        let mut rng = rand::rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..2000).map(|_| rng.random::<u64>()).collect();
        let seed_vec: Vec<u64> = kmer_hash_set.into_iter().collect();
        for hv_d in [1024usize, 4096, 16384] {
            assert_eq!(
                hash_hv_bit_serial(&seed_vec, hv_d),
                hash_hv_bit(&seed_vec, hv_d),
                "AVX-512 bit encoding mismatch at dim {hv_d}"
            );
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn test_hash_hv_bit_avx2_serial_vs_simd() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let mut rng = rand::rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..2000).map(|_| rng.random::<u64>()).collect();
        let seed_vec: Vec<u64> = kmer_hash_set.into_iter().collect();
        // 1056 = 33×32 + 0 tail; 1064 = 33×32 + 8 tail (tail keeps −N).
        for hv_d in [1024usize, 4096, 16384, 1056, 1064] {
            let simd = unsafe { hash_hv_bit_avx2(&seed_vec, hv_d) };
            assert_eq!(
                hash_hv_bit_serial(&seed_vec, hv_d),
                simd,
                "AVX2 bit encoding mismatch at dim {hv_d}"
            );
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn test_hash_hv_i8_avx2_serial_vs_simd() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let mut rng = rand::rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..2000).map(|_| rng.random::<u64>()).collect();
        let seed_vec: Vec<u64> = kmer_hash_set.into_iter().collect();
        for hv_d in [1024usize, 4096, 16384] {
            let simd = unsafe { hash_hv_i8_avx2(&seed_vec, hv_d) };
            assert_eq!(
                hash_hv_i8_serial(&seed_vec, hv_d),
                simd,
                "AVX2 i8 encoding mismatch at dim {hv_d}"
            );
        }
    }

    #[test]
    fn test_hv_norm_l2_sq() {
        // Create a simple hypervector
        let hv = vec![1, 2, 3, 4, 5];

        // Compute the squared L2 norm
        let norm_sq = hv_norm_l2_sq(&hv);

        // Expected result: 1^2 + 2^2 + 3^2 + 4^2 + 5^2 = 55
        assert_eq!(norm_sq, 55.0, "Squared L2 norm calculation is incorrect!");
    }

    #[test]
    fn test_hv_norm_l2_sq_serial_vs_simd() {
        let hv: Vec<_> = (1..=32).collect();

        let result_scalar = hv_norm_l2_sq_serial(&hv);
        let result_simd = hv_norm_l2_sq(&hv);

        println!("Scalar result: {}", result_scalar);
        println!("SIMD result: {}", result_simd);

        assert_eq!(result_scalar, result_simd, "Results do not match!");
    }

    #[test]
    fn test_hv_cardinality() {
        // Create a simple hypervector
        let hv = vec![1, 2, 3, 4, 5];

        // Compute the cardinality
        let cardinality = hv_cardinality(&hv);

        // Expected result: (1^2 + 2^2 + 3^2 + 4^2 + 5^2) / 5 = 55 / 5 = 11
        assert_eq!(cardinality, 11, "Cardinality calculation is incorrect!");
    }

    #[test]
    fn test_hv_dot() {
        // Create two simple hypervectors
        let a = vec![1, 2, 3, 4, 5];
        let b = vec![2, 3, 4, 5, 6];

        // Compute the dot product
        let dot = hv_dot(&a, &b);

        // Expected result: (1*2 + 2*3 + 3*4 + 4*5 + 5*6) / 5 = 14
        assert_eq!(dot, 14.0, "Dot product calculation is incorrect!");
    }

    #[test]
    fn test_hv_dot_orthogonal() {
        // Create two orthogonal hypervectors
        let a = vec![1, 0, 0];
        let b = vec![0, 1, 0];

        // Compute the dot product
        let dot = hv_dot(&a, &b);

        // Expected result: (1*0 + 0*1 + 0*0) / 3 = 0
        assert_eq!(
            dot, 0.0,
            "Dot product of orthogonal vectors should be zero!"
        );
    }
}
