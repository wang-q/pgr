use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::hv::{hash_hv_bit, hash_hv_i8};
use rand::rngs::{SmallRng, StdRng};
use rand::{Rng, RngCore, SeedableRng};
use rapidhash::{RapidHashSet, RapidRng};

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
    let num_chunk = hv_d / 32;
    let one = _mm512_set1_epi32(1);
    let shifts_lo = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    let shifts_hi = _mm512_set_epi32(
        31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16,
    );
    for b in 0..num_chunk {
        let j = b as u64 + 1;
        let mut acc_lo = _mm512_setzero_si512();
        let mut acc_hi = _mm512_setzero_si512();
        for &seed in seed_vec {
            let r = rnd_at(seed, j) as u32;
            let v = _mm512_set1_epi32(r as i32);
            let bits_lo = _mm512_and_si512(_mm512_srlv_epi32(v, shifts_lo), one);
            let signed_lo = _mm512_sub_epi32(_mm512_slli_epi32(bits_lo, 1), one);
            acc_lo = _mm512_add_epi32(acc_lo, signed_lo);
            let bits_hi = _mm512_and_si512(_mm512_srlv_epi32(v, shifts_hi), one);
            let signed_hi = _mm512_sub_epi32(_mm512_slli_epi32(bits_hi, 1), one);
            acc_hi = _mm512_add_epi32(acc_hi, signed_hi);
        }
        _mm512_storeu_si512(hv[b * 32..b * 32 + 16].as_mut_ptr() as *mut _, acc_lo);
        _mm512_storeu_si512(hv[b * 32 + 16..b * 32 + 32].as_mut_ptr() as *mut _, acc_hi);
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
    // Create test datasets of different sizes
    let kmer_hash_set_small = generate_kmer_hash_set(1000); // Small dataset
    let kmer_hash_set_medium = generate_kmer_hash_set(10_000); // Medium dataset

    let seed_vec_small: Vec<u64> = kmer_hash_set_small.iter().cloned().collect();
    let seed_vec_medium: Vec<u64> = kmer_hash_set_medium.iter().cloned().collect();

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
}

// Define benchmark group
criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10); // Set sample size
    targets = bench_encode_hash_hd
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
