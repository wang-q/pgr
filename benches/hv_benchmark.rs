use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::hv::{hash_hv_bit, hash_hv_i8, hash_hv_sparse};
use rand::rngs::{SmallRng, StdRng};
use rand::{Rng, RngCore, SeedableRng};
use rapidhash::{RapidHashSet, RapidRng};
use wide::{bytemuck, i32x8, u16x8, u32x8, u8x16};

// Note: `hash_hv_bit` / `hash_hv_i8` dispatch to AVX2 (block-major jump-ahead
// RapidRng, see src/libs/hv.rs) on x86-64, else the portable wide path.
// AVX-512 implementations are kept below as reference only (not part of the
// runtime dispatch); the `*_avx512_ref` benchmarks compare them side by side
// on AVX-512-capable CPUs.

// ---------------------------------------------------------------------------
// AVX-512 reference implementations (same RNG stream as the AVX2 path).
// ---------------------------------------------------------------------------

const RAPID_S0: u64 = 0x2d358dccaa6c78a5;
const RAPID_S1: u64 = 0x8bb84b93962eacc9;

#[inline(always)]
fn rapid_mix(a: u64, b: u64) -> u64 {
    let r = (a as u128) * (b as u128);
    (r as u64) ^ ((r >> 64) as u64)
}

#[inline(always)]
fn rnd_at(seed: u64, j: u64) -> u64 {
    let s = seed.wrapping_add(j.wrapping_mul(RAPID_S0));
    rapid_mix(s, s ^ RAPID_S1)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl")]
unsafe fn bit_avx512_impl(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    use std::arch::x86_64::*;
    let n = seed_vec.len() as i32;
    let mut hv = vec![-n; hv_d];
    let num_chunk = hv_d / 64;
    let one = _mm512_set1_epi32(1);
    let shifts_lo = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    let shifts_hi = _mm512_set_epi32(
        31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16,
    );
    for b in 0..num_chunk {
        let j = b as u64 + 1;
        let mut acc = [_mm512_setzero_si512(); 4];
        for &seed in seed_vec {
            let r = rnd_at(seed, j);
            let vlo = _mm512_set1_epi32((r as u32) as i32);
            let vhi = _mm512_set1_epi32(((r >> 32) as u32) as i32);
            for k in 0..2 {
                let shift = if k == 0 { shifts_lo } else { shifts_hi };
                let bl = _mm512_and_si512(_mm512_srlv_epi32(vlo, shift), one);
                acc[k] = _mm512_add_epi32(acc[k], _mm512_sub_epi32(_mm512_slli_epi32(bl, 1), one));
                let bh = _mm512_and_si512(_mm512_srlv_epi32(vhi, shift), one);
                acc[k + 2] =
                    _mm512_add_epi32(acc[k + 2], _mm512_sub_epi32(_mm512_slli_epi32(bh, 1), one));
            }
        }
        let base = b * 64;
        for k in 0..4 {
            _mm512_storeu_si512(
                hv[base + k * 16..base + k * 16 + 16].as_mut_ptr() as *mut _,
                acc[k],
            );
        }
    }
    hv
}

pub fn bit_avx512_ref(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx512f")
        && is_x86_feature_detected!("avx512bw")
        && is_x86_feature_detected!("avx512vl")
    {
        return unsafe { bit_avx512_impl(seed_vec, hv_d) };
    }
    hash_hv_bit(seed_vec, hv_d)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl")]
unsafe fn i8_avx512_impl(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    use std::arch::x86_64::*;
    let mut hv = vec![0i32; hv_d];
    let num_chunk = hv_d / 16;
    for b in 0..num_chunk {
        let j1 = (2 * b) as u64 + 1;
        let j2 = (2 * b) as u64 + 2;
        let mut acc = _mm512_setzero_si512();
        for &seed in seed_vec {
            let r1 = rnd_at(seed, j1).to_ne_bytes();
            let r2 = rnd_at(seed, j2).to_ne_bytes();
            let mut arr = [0u8; 16];
            arr[..8].copy_from_slice(&r1);
            arr[8..].copy_from_slice(&r2);
            let bytes = _mm_loadu_si128(arr.as_ptr() as *const __m128i);
            acc = _mm512_add_epi32(acc, _mm512_cvtepi8_epi32(bytes));
        }
        _mm512_storeu_si512(hv[b * 16..b * 16 + 16].as_mut_ptr() as *mut _, acc);
    }
    hv
}

pub fn i8_avx512_ref(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx512f")
        && is_x86_feature_detected!("avx512bw")
        && is_x86_feature_detected!("avx512vl")
    {
        return unsafe { i8_avx512_impl(seed_vec, hv_d) };
    }
    hash_hv_i8(seed_vec, hv_d)
}

// ---------------------------------------------------------------------------
// RNG candidate comparison (2026-08-08): block-major AVX2 bit encoding with
// alternative jump-ahead RNGs vs the RapidRng baseline (`hash_hv_bit`). The
// three candidates share the counter + mix structure required for O(1) jump
// ahead: constant (zero-cost upper bound), splitmix64, wyrand. Bodies are
// instruction-identical to `hash_hv_bit_avx2` except for the RNG line.
// ---------------------------------------------------------------------------

#[inline(always)]
fn rnd_const(_seed: u64, _j: u64) -> u64 {
    0x9E37_79B9_7F4A_7C15
}

#[inline(always)]
fn rnd_splitmix(seed: u64, j: u64) -> u64 {
    let mut x = seed.wrapping_add(j.wrapping_mul(0x9E37_79B9_7F4A_7C15));
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

#[inline(always)]
fn rnd_wyrand(seed: u64, j: u64) -> u64 {
    let mut s = seed.wrapping_add(j.wrapping_mul(0xA076_1D64_78BD_642F));
    s ^= s >> 32;
    s.wrapping_mul(0xE703_7ED1_A0B4_28DB)
}

#[inline(always)]
fn rnd_raw_seed(seed: u64, _j: u64) -> u64 {
    // No mix: isolates whether the RNG mix itself or the per-seed broadcast
    // data-flow is the dominant cost (invalid as a real RNG, measurement only).
    seed
}

/// MINSTD LCG (a=16807, m=2^31-1) with O(log j) jump-ahead via modular
/// exponentiation: out(j) = seed·a^j mod m.
#[inline(always)]
fn rnd_lcg(seed: u64, j: u64) -> u64 {
    const A: u64 = 16807;
    const M: u64 = 2147483647;
    let mut result = 1u64;
    let mut base = A;
    let mut e = j;
    while e > 0 {
        if e & 1 == 1 {
            result = (result as u128 * base as u128 % M as u128) as u64;
        }
        base = (base as u128 * base as u128 % M as u128) as u64;
        e >>= 1;
    }
    (seed as u128 * result as u128 % M as u128) as u64
}

/// PCG-XSH-RR 64-bit with O(log j) jump-ahead over the underlying LCG
/// (a=6364136223846793005, c=1442695040888963407, modulus 2^64).
#[inline(always)]
fn rnd_pcg(seed: u64, j: u64) -> u64 {
    const A: u64 = 6364136223846793005;
    const C: u64 = 1442695040888963407;
    // Jump-ahead: state_j = seed·A^j + C·Σ_{i< j}A^i (mod 2^64).
    let mut a_pow = 1u64;
    let mut sum = 0u64;
    let mut base = A;
    let mut base_sum = 1u64;
    let mut e = j;
    while e > 0 {
        if e & 1 == 1 {
            a_pow = a_pow.wrapping_mul(base);
            sum = sum.wrapping_mul(base).wrapping_add(base_sum);
        }
        base_sum = base_sum.wrapping_mul(base).wrapping_add(base_sum);
        base = base.wrapping_mul(base);
        e >>= 1;
    }
    let state = seed.wrapping_mul(a_pow).wrapping_add(C.wrapping_mul(sum));
    let xorshifted = (((state >> 18) ^ state) >> 27) as u32;
    let rot = (state >> 59) as u32;
    xorshifted.rotate_right(rot) as u64
}

macro_rules! bit_avx2_rng_variant {
    ($name:ident, $rnd:expr) => {
        #[cfg(target_arch = "x86_64")]
        #[target_feature(enable = "avx2")]
        unsafe fn $name(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
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
                    // black_box prevents the compiler from hoisting the
                    // broadcast for constant RNGs (LICM), so every variant
                    // pays the per-seed `set1` cost like a real RNG stream.
                    let r = black_box($rnd(seed, j));
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
            for v in &mut hv[num_chunk * 64..] {
                *v = -n;
            }
            hv
        }
    };
}

bit_avx2_rng_variant!(bit_avx2_rnd_const, rnd_const);
bit_avx2_rng_variant!(bit_avx2_rnd_splitmix, rnd_splitmix);
bit_avx2_rng_variant!(bit_avx2_rnd_wyrand, rnd_wyrand);
bit_avx2_rng_variant!(bit_avx2_rnd_raw_seed, rnd_raw_seed);
bit_avx2_rng_variant!(bit_avx2_rnd_rapid, rnd_at);
bit_avx2_rng_variant!(bit_avx2_rnd_lcg, rnd_lcg);
bit_avx2_rng_variant!(bit_avx2_rnd_pcg, rnd_pcg);

/// Statistical sanity check for the 64-bit variant: both halves of `rnd_at`
/// output must behave like independent balanced bit sources (±1 encoding
/// needs ~50% ones; independence gives ~25% lo&hi overlap).
fn verify_64bit_halves() {
    let mut rng = StdRng::seed_from_u64(7);
    let n = 100_000u64;
    let mut lo_ones = 0u64;
    let mut hi_ones = 0u64;
    let mut cross = 0u64;
    for _ in 0..n {
        let r = rnd_at(rng.random(), rng.random::<u64>() % 1000 + 1);
        let lo = r as u32;
        let hi = (r >> 32) as u32;
        lo_ones += lo.count_ones() as u64;
        hi_ones += hi.count_ones() as u64;
        cross += (lo & hi).count_ones() as u64;
    }
    let total = n as f64 * 32.0;
    let (p_lo, p_hi, p_cross) = (
        lo_ones as f64 / total,
        hi_ones as f64 / total,
        cross as f64 / total,
    );
    assert!((p_lo - 0.5).abs() < 0.01, "lo half imbalanced: {p_lo}");
    assert!((p_hi - 0.5).abs() < 0.01, "hi half imbalanced: {p_hi}");
    assert!(
        (p_cross - 0.25).abs() < 0.01,
        "lo/hi halves correlated: {p_cross}"
    );
}

pub fn encode_hash_hd_rapid(seed_vec: &[u64], hv_d: usize) -> Vec<i16> {
    let mut hv = vec![-(seed_vec.len() as i16); hv_d];

    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);

        for i in 0..(hv_d / 64) {
            let rnd_bits = rng.next_u64();

            for j in 0..64 {
                hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
            }
        }
    }

    hv
}

pub fn encode_hash_hd_std(seed_vec: &[u64], hv_d: usize) -> Vec<i16> {
    let mut hv = vec![-(seed_vec.len() as i16); hv_d];

    for hash in seed_vec {
        let mut rng = StdRng::seed_from_u64(*hash);

        for i in 0..(hv_d / 64) {
            let rnd_bits = rng.next_u64();

            for j in 0..64 {
                hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
            }
        }
    }

    hv
}

pub fn encode_hash_hd_small(seed_vec: &[u64], hv_d: usize) -> Vec<i16> {
    let mut hv = vec![-(seed_vec.len() as i16); hv_d];

    for hash in seed_vec {
        let mut rng = SmallRng::seed_from_u64(*hash);

        for i in 0..(hv_d / 64) {
            let rnd_bits = rng.next_u64();

            for j in 0..64 {
                hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
            }
        }
    }

    hv
}

// ---------------------------------------------------------------------------
// Scalar RNG baselines (2026-08-08): classic / alternative generators in the
// same per-seed streaming pattern as encode_hash_hd_rapid/std/small. Each
// generator is seeded once per k-mer hash and produces the D random bits.
// ---------------------------------------------------------------------------

/// MT19937 (Mersenne Twister), standard 32-bit implementation.
struct Mt19937 {
    mt: [u32; 624],
    index: usize,
}

impl Mt19937 {
    #[inline(always)]
    fn seed_from_u64(seed: u64) -> Self {
        let mut mt = [0u32; 624];
        mt[0] = seed as u32;
        for i in 1..624 {
            mt[i] = 1812433253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        let mut rng = Mt19937 { mt, index: 624 };
        rng.next_u32();
        rng
    }

    #[inline(always)]
    fn next_u32(&mut self) -> u32 {
        if self.index >= 624 {
            for i in 0..624 {
                let y = (self.mt[i] & 0x8000_0000) | (self.mt[(i + 1) % 624] & 0x7fff_ffff);
                self.mt[i] = self.mt[(i + 397) % 624]
                    ^ (y >> 1)
                    ^ (if y & 1 != 0 { 0x9908_b0df } else { 0 });
            }
            self.index = 0;
        }
        let mut y = self.mt[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    #[inline(always)]
    fn next_u64(&mut self) -> u64 {
        let lo = self.next_u32() as u64;
        let hi = self.next_u32() as u64;
        (hi << 32) | lo
    }
}

/// PCG-XSH-RR 32-bit (classic LCG + output function).
struct Pcg32 {
    state: u64,
}

impl Pcg32 {
    #[inline(always)]
    fn seed_from_u64(seed: u64) -> Self {
        const A: u64 = 6364136223846793005;
        const C: u64 = 1442695040888963407;
        let mut rng = Pcg32 { state: 0 };
        rng.state = rng.state.wrapping_mul(A).wrapping_add(C);
        rng.state = rng.state.wrapping_add(seed);
        rng.state = rng.state.wrapping_mul(A).wrapping_add(C);
        rng
    }

    #[inline(always)]
    fn next_u32(&mut self) -> u32 {
        const A: u64 = 6364136223846793005;
        const C: u64 = 1442695040888963407;
        let old = self.state;
        self.state = old.wrapping_mul(A).wrapping_add(C);
        let xorshifted = (((old >> 18) ^ old) >> 27) as u32;
        let rot = (old >> 59) as u32;
        xorshifted.rotate_right(rot)
    }

    #[inline(always)]
    fn next_u64(&mut self) -> u64 {
        let lo = self.next_u32() as u64;
        let hi = self.next_u32() as u64;
        (hi << 32) | lo
    }
}

/// xorshift64* (Marsaglia).
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    #[inline(always)]
    fn seed_from_u64(seed: u64) -> Self {
        Xorshift64 {
            state: if seed == 0 {
                0x9E37_79B9_7F4A_7C15
            } else {
                seed
            },
        }
    }

    #[inline(always)]
    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
}

/// wyrand (HyperGen's HDC encoding RNG).
struct Wyrand {
    state: u64,
}

impl Wyrand {
    #[inline(always)]
    fn seed_from_u64(seed: u64) -> Self {
        Wyrand { state: seed }
    }

    #[inline(always)]
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0xA076_1D64_78BD_642F);
        let mut s = self.state;
        s ^= s >> 32;
        s.wrapping_mul(0xE703_7ED1_A0B4_28DB)
    }
}

/// splitmix64.
struct Splitmix64 {
    state: u64,
}

impl Splitmix64 {
    #[inline(always)]
    fn seed_from_u64(seed: u64) -> Self {
        Splitmix64 { state: seed }
    }

    #[inline(always)]
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut x = self.state;
        x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        x ^ (x >> 31)
    }
}

macro_rules! encode_hash_hd_scalar {
    ($name:ident, $rng:ty) => {
        pub fn $name(seed_vec: &[u64], hv_d: usize) -> Vec<i16> {
            let mut hv = vec![-(seed_vec.len() as i16); hv_d];
            for hash in seed_vec {
                let mut rng = <$rng>::seed_from_u64(*hash);
                for i in 0..(hv_d / 64) {
                    let rnd_bits = rng.next_u64();
                    for j in 0..64 {
                        hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
                    }
                }
            }
            hv
        }
    };
}

encode_hash_hd_scalar!(encode_hash_hd_mt19937, Mt19937);
encode_hash_hd_scalar!(encode_hash_hd_pcg32, Pcg32);
encode_hash_hd_scalar!(encode_hash_hd_xorshift, Xorshift64);
encode_hash_hd_scalar!(encode_hash_hd_wyrand, Wyrand);
encode_hash_hd_scalar!(encode_hash_hd_splitmix, Splitmix64);

/// Scalar i8 baseline: per-seed serial RapidRng stream, one byte per 8 dims
/// (mirrors `hash_hv_i8_serial` semantics; SIMD speedup reference for hv.md
/// §2.1 — the RNG stream is scalar and dominates, so SIMD gains are modest).
pub fn encode_hash_hd_i8_serial(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    let mut hv = vec![0i32; hv_d];
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);
        for i in 0..(hv_d / 8) {
            let rnd_bits = rng.next_u64();
            let bytes = rnd_bits.to_ne_bytes();
            for j in 0..8 {
                hv[i * 8 + j] += bytes[j] as i8 as i32;
            }
        }
    }
    hv
}

// ---------------------------------------------------------------------------
// Portable wide fallback baselines (2026-08-08): explicit copies of the wide
// branches in `src/libs/hv.rs` (used on non-x86_64 / non-AVX2 CPUs, aarch64
// NEON etc.), measured directly so hv.md §2.1 can compare AVX2 / wide /
// scalar three-way.
// ---------------------------------------------------------------------------

pub fn hash_hv_bit_wide(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    let num_seed = seed_vec.len();
    let num_chunk = hv_d / 64;
    let mut hv = vec![-(num_seed as i32); hv_d];
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);
        for i in 0..num_chunk {
            let rnd_bits = rng.next_u64();
            let halves = [(rnd_bits as u32), (rnd_bits >> 32) as u32];
            for (k, half) in halves.iter().enumerate() {
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

pub fn hash_hv_i8_wide(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    let mut hv = vec![0i32; hv_d];
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);
        for i in 0..(hv_d / 8) {
            let rnd_bits = rng.next_u64();
            let bytes = rnd_bits.to_ne_bytes();
            let mut arr = [0u8; 16];
            arr[..8].copy_from_slice(&bytes);
            let vec_u8 = u8x16::from(arr);
            let vec_u16 = u16x8::from_u8x16_low(vec_u8);
            let vec_i32 = i32x8::from_u16x8(vec_u16);
            let vec_vals = (vec_i32 << i32x8::splat(24)) >> i32x8::splat(24);
            let mut hv_simd = i32x8::from(&hv[i * 8..(i + 1) * 8]);
            hv_simd += vec_vals;
            hv[i * 8..(i + 1) * 8].copy_from_slice(&hv_simd.to_array());
        }
    }
    hv
}

/// AVX2 bit encoding with **i16 accumulators** (2026-08-08): same 64-dim
/// block-major structure as `hash_hv_bit_avx2`, but the eight i32 8-lane
/// accumulators become four i16 16-lane ones (bits 0-15 / 16-31 / 32-47 /
/// 48-63 of each `rnd_at` output). Deferred −N and final i16→i32 widening
/// keep the i32 HV layout. Safe for n ≤ 32767 without segmentation (n=10k
/// here); larger n requires the segmented scheme (08-07 experiment).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn hash_hv_bit_i16(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    use std::arch::x86_64::*;
    let n = seed_vec.len() as i32;
    let mut hv = vec![0i32; hv_d];
    let num_chunk = hv_d / 64;
    let one = _mm256_set1_epi16(1);
    let shifts = _mm256_set_epi16(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    let nv = _mm256_set1_epi16(n as i16);
    for b in 0..num_chunk {
        let j = b as u64 + 1;
        let mut a = [_mm256_setzero_si256(); 4];
        for &seed in seed_vec {
            let r = rnd_at(seed, j);
            for (k, acc) in a.iter_mut().enumerate() {
                let half = (r >> (16 * k)) as u16;
                let v = _mm256_set1_epi16(half as i16);
                let bits = _mm256_and_si256(_mm256_srlv_epi16(v, shifts), one);
                *acc = _mm256_add_epi16(*acc, _mm256_slli_epi16(bits, 1));
            }
        }
        let base = b * 64;
        for (k, acc) in a.iter().enumerate() {
            let sub = _mm256_sub_epi16(*acc, nv);
            let lo = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(sub));
            let hi = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(sub, 1));
            _mm256_storeu_si256(
                hv[base + k * 16..base + k * 16 + 8].as_mut_ptr() as *mut _,
                lo,
            );
            _mm256_storeu_si256(
                hv[base + k * 16 + 8..base + k * 16 + 16].as_mut_ptr() as *mut _,
                hi,
            );
        }
    }
    for v in &mut hv[num_chunk * 64..] {
        *v = -n;
    }
    hv
}

/// AVX2 bit encoding via **pshufb 4-bit LUT expansion** (2026-08-08): each
/// 64-bit `rnd_at` output is split into 16 nibbles; per nibble-bit position
/// one `pshufb` expands 16 nibbles into 16 0/1 bytes, widened with
/// `vpmovzxbd` and accumulated as ±1 (deferred −N). Tests whether table
/// lookup beats the srlv+and+slli expansion (expectation: more ops per bit,
/// since byte-level LUTs need extra widening).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn hash_hv_bit_pshufb(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    use std::arch::x86_64::*;
    let n = seed_vec.len() as i32;
    let mut hv = vec![0i32; hv_d];
    let num_chunk = hv_d / 64;
    // LUT for each nibble bit position: T[k] = bit value (0/1 byte).
    let t0 = _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1);
    let t1 = _mm_setr_epi8(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1);
    let t2 = _mm_setr_epi8(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1);
    let t3 = _mm_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1);
    let nv = _mm256_set1_epi32(n);
    for b in 0..num_chunk {
        let j = b as u64 + 1;
        let mut acc = [_mm256_setzero_si256(); 8];
        for &seed in seed_vec {
            let r = rnd_at(seed, j);
            let mut nidx = [0u8; 16];
            for (k, v) in nidx.iter_mut().enumerate() {
                *v = ((r >> (4 * k)) & 0xF) as u8;
            }
            let idx = _mm_loadu_si128(nidx.as_ptr() as *const __m128i);
            for (bit, tab) in [t0, t1, t2, t3].iter().enumerate() {
                let out = _mm_shuffle_epi8(*tab, idx);
                let lo = _mm256_cvtepu8_epi32(out);
                let hi = _mm256_cvtepu8_epi32(_mm_srli_si128(out, 8));
                acc[bit * 2] = _mm256_add_epi32(acc[bit * 2], _mm256_slli_epi32(lo, 1));
                acc[bit * 2 + 1] = _mm256_add_epi32(acc[bit * 2 + 1], _mm256_slli_epi32(hi, 1));
            }
        }
        let base = b * 64;
        for (k, a) in acc.iter().enumerate() {
            _mm256_storeu_si256(
                hv[base + k * 8..base + k * 8 + 8].as_mut_ptr() as *mut _,
                _mm256_sub_epi32(*a, nv),
            );
        }
    }
    for v in &mut hv[num_chunk * 64..] {
        *v = -n;
    }
    hv
}

// Generate a random k-mer hash set
fn generate_kmer_hash_set(size: usize) -> RapidHashSet<u64> {
    let mut rng = StdRng::seed_from_u64(42); // Fixed seed for reproducibility
    let mut kmer_hash_set = RapidHashSet::default();

    for _ in 0..size {
        kmer_hash_set.insert(rng.random::<u64>());
    }

    kmer_hash_set
}

// Benchmark function
fn bench_encode_hash_hd(c: &mut Criterion) {
    // Statistical sanity of the high 32 bits before benchmarking the
    // 64-bit-consuming variant.
    verify_64bit_halves();

    // Create test datasets of different sizes
    let kmer_hash_set_small = generate_kmer_hash_set(1000); // Small dataset
    let kmer_hash_set_medium = generate_kmer_hash_set(10_000); // Medium dataset
    let kmer_hash_set_large = generate_kmer_hash_set(100_000); // Large dataset

    let seed_vec_small: Vec<u64> = kmer_hash_set_small.iter().cloned().collect();
    let seed_vec_medium: Vec<u64> = kmer_hash_set_medium.iter().cloned().collect();
    let seed_vec_large: Vec<u64> = kmer_hash_set_large.iter().cloned().collect();

    let hv_d = 4096; // Hypervector dimension

    // Benchmark small dataset
    c.bench_function("encode_hash_hd_lib_small", |b| {
        b.iter(|| hash_hv_bit(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_simd_i8_small", |b| {
        b.iter(|| hash_hv_i8(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_rapid_small", |b| {
        b.iter(|| encode_hash_hd_rapid(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_std_small", |b| {
        b.iter(|| encode_hash_hd_std(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_small_small", |b| {
        b.iter(|| encode_hash_hd_small(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("avx512_ref_bit_small", |b| {
        b.iter(|| bit_avx512_ref(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("avx512_ref_i8_small", |b| {
        b.iter(|| i8_avx512_ref(black_box(&seed_vec_small), hv_d))
    });

    // Benchmark medium dataset
    c.bench_function("encode_hash_hd_lib_medium", |b| {
        b.iter(|| hash_hv_bit(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_simd_i8_medium", |b| {
        b.iter(|| hash_hv_i8(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_rapid_medium", |b| {
        b.iter(|| encode_hash_hd_rapid(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_std_medium", |b| {
        b.iter(|| encode_hash_hd_std(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_small_medium", |b| {
        b.iter(|| encode_hash_hd_small(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("avx512_ref_bit_medium", |b| {
        b.iter(|| bit_avx512_ref(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("avx512_ref_i8_medium", |b| {
        b.iter(|| i8_avx512_ref(black_box(&seed_vec_medium), hv_d))
    });

    // Scalar RNG baselines (classic / alternative generators, medium only:
    // the small vs medium ratio is stable across generators, see §2 tables).
    c.bench_function("encode_hash_hd_mt19937_medium", |b| {
        b.iter(|| encode_hash_hd_mt19937(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_pcg32_medium", |b| {
        b.iter(|| encode_hash_hd_pcg32(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_xorshift_medium", |b| {
        b.iter(|| encode_hash_hd_xorshift(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_wyrand_medium", |b| {
        b.iter(|| encode_hash_hd_wyrand(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_splitmix_medium", |b| {
        b.iter(|| encode_hash_hd_splitmix(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_i8_serial_medium", |b| {
        b.iter(|| encode_hash_hd_i8_serial(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("hash_hv_bit_wide_medium", |b| {
        b.iter(|| hash_hv_bit_wide(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("hash_hv_i8_wide_medium", |b| {
        b.iter(|| hash_hv_i8_wide(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("hash_hv_bit_i16_medium", |b| {
        b.iter(|| unsafe { hash_hv_bit_i16(black_box(&seed_vec_medium), hv_d) })
    });
    c.bench_function("hash_hv_bit_pshufb_medium", |b| {
        b.iter(|| unsafe { hash_hv_bit_pshufb(black_box(&seed_vec_medium), hv_d) })
    });

    // Sparse projection (`.hv` v2 path): s random dims per seed, splitmix64.
    for s in [1usize, 3, 8, 16, 64] {
        c.bench_function(&format!("hash_hv_sparse_s{}_medium", s), |b| {
            b.iter(|| hash_hv_sparse(black_box(&seed_vec_medium), hv_d, s))
        });
    }
    c.bench_function("hash_hv_sparse_s3_large", |b| {
        b.iter(|| hash_hv_sparse(black_box(&seed_vec_large), hv_d, 3))
    });
    // Encoding time vs D at fixed s: sparse cost is O(n·s), independent of D
    // (only the HV array grows). Key advantage over dense (O(n·D)) when
    // pushing precision via larger D.
    for d in [16384usize, 65536] {
        c.bench_function(&format!("hash_hv_sparse_s1_d{}_medium", d), |b| {
            b.iter(|| hash_hv_sparse(black_box(&seed_vec_medium), d, 1))
        });
    }

    // Large dataset (n=100k, D=4096)
    c.bench_function("encode_hash_hd_lib_large", |b| {
        b.iter(|| hash_hv_bit(black_box(&seed_vec_large), hv_d))
    });
    c.bench_function("encode_hash_hd_simd_i8_large", |b| {
        b.iter(|| hash_hv_i8(black_box(&seed_vec_large), hv_d))
    });
    c.bench_function("avx512_ref_bit_large", |b| {
        b.iter(|| bit_avx512_ref(black_box(&seed_vec_large), hv_d))
    });
    c.bench_function("avx512_ref_i8_large", |b| {
        b.iter(|| i8_avx512_ref(black_box(&seed_vec_large), hv_d))
    });

    // D = 16384 variants on the medium (10k) seed set
    let hv_d_16k = 16384;
    c.bench_function("encode_hash_hd_lib_d16k", |b| {
        b.iter(|| hash_hv_bit(black_box(&seed_vec_medium), hv_d_16k))
    });
    c.bench_function("encode_hash_hd_simd_i8_d16k", |b| {
        b.iter(|| hash_hv_i8(black_box(&seed_vec_medium), hv_d_16k))
    });
    c.bench_function("avx512_ref_bit_d16k", |b| {
        b.iter(|| bit_avx512_ref(black_box(&seed_vec_medium), hv_d_16k))
    });
    c.bench_function("avx512_ref_i8_d16k", |b| {
        b.iter(|| i8_avx512_ref(black_box(&seed_vec_medium), hv_d_16k))
    });

    // RNG candidate comparison vs the RapidRng baseline (`hash_hv_bit`):
    // constant (zero-cost upper bound), splitmix64 and wyrand jump-ahead.
    c.bench_function("bit_avx2_rng_const_medium", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_const(black_box(&seed_vec_medium), hv_d) })
    });
    c.bench_function("bit_avx2_rng_splitmix_medium", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_splitmix(black_box(&seed_vec_medium), hv_d) })
    });
    c.bench_function("bit_avx2_rng_wyrand_medium", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_wyrand(black_box(&seed_vec_medium), hv_d) })
    });
    c.bench_function("bit_avx2_rng_const_large", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_const(black_box(&seed_vec_large), hv_d) })
    });
    c.bench_function("bit_avx2_rng_splitmix_large", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_splitmix(black_box(&seed_vec_large), hv_d) })
    });
    c.bench_function("bit_avx2_rng_wyrand_large", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_wyrand(black_box(&seed_vec_large), hv_d) })
    });
    c.bench_function("bit_avx2_rng_raw_seed_medium", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_raw_seed(black_box(&seed_vec_medium), hv_d) })
    });
    c.bench_function("bit_avx2_rng_raw_seed_large", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_raw_seed(black_box(&seed_vec_large), hv_d) })
    });
    c.bench_function("bit_avx2_rng_rapid_medium", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_rapid(black_box(&seed_vec_medium), hv_d) })
    });
    c.bench_function("bit_avx2_rng_rapid_large", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_rapid(black_box(&seed_vec_large), hv_d) })
    });
    c.bench_function("bit_avx2_rng_lcg_medium", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_lcg(black_box(&seed_vec_medium), hv_d) })
    });
    c.bench_function("bit_avx2_rng_lcg_large", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_lcg(black_box(&seed_vec_large), hv_d) })
    });
    c.bench_function("bit_avx2_rng_pcg_medium", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_pcg(black_box(&seed_vec_medium), hv_d) })
    });
    c.bench_function("bit_avx2_rng_pcg_large", |b| {
        b.iter(|| unsafe { bit_avx2_rnd_pcg(black_box(&seed_vec_large), hv_d) })
    });
}

/// Sampling-hash throughput (2026-08-08): t1ha2 (HyperGen's FracMinHash
/// sampler) and wyhash vs pgr's minimizer hashers, on 21-mer byte slices.
fn bench_hash_throughput(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(42);
    let kmers: Vec<Vec<u8>> = (0..10_000)
        .map(|_| (0..21).map(|_| rng.random_range(0..4) as u8).collect())
        .collect();

    c.bench_function("hash_t1ha2_21mer", |b| {
        b.iter(|| {
            let mut acc = 0u64;
            for k in &kmers {
                acc = acc.wrapping_add(t1ha::t1ha2_atonce(k, 42));
            }
            black_box(acc);
        })
    });
    c.bench_function("hash_wyhash_21mer", |b| {
        b.iter(|| {
            let mut acc = 0u64;
            for k in &kmers {
                acc = acc.wrapping_add(wyhash::wyhash(k, 42));
            }
            black_box(acc);
        })
    });
    c.bench_function("hash_rapidhash_21mer", |b| {
        b.iter(|| {
            let mut acc = 0u64;
            for k in &kmers {
                acc = acc.wrapping_add(rapidhash::rapidhash(k));
            }
            black_box(acc);
        })
    });
    c.bench_function("hash_fxhash_21mer", |b| {
        b.iter(|| {
            let mut acc = 0u64;
            for k in &kmers {
                acc = acc.wrapping_add(fxhash::hash64(k));
            }
            black_box(acc);
        })
    });
    c.bench_function("hash_murmur3_21mer", |b| {
        b.iter(|| {
            let mut acc = 0u64;
            for k in &kmers {
                acc = acc.wrapping_add(murmurhash3::murmurhash3_x64_128(k, 42).0);
            }
            black_box(acc);
        })
    });
}

// Define benchmark group
criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10); // Set sample size
    targets = bench_encode_hash_hd, bench_hash_throughput
);

// Run benchmarks
criterion_main!(benches);

// Historical reference (hnsm/benches/hd.rs, 2026-01-30):
// Small (1000):  encode_hash_hd_simd_i8 ~421 µs | lib(bit) ~670 µs |
//                rapid ~1.66 ms | std ~1.75 ms | small ~1.64 ms
// Medium (10000): encode_hash_hd_simd_i8 ~4.20 ms | lib(bit) ~6.73 ms |
//                 rapid ~16.5 ms | std ~17.5 ms | small ~16.3 ms
// Note: the exploratory SIMD variants (u64x4/i16x4, u8x8/i32) from hnsm were
// merged into the current lib implementations (hash_hv_bit / hash_hv_i8);
// wide does not provide i16x4/u8x8, so only the lib + scalar-RNG baselines
// are reproduced here.
