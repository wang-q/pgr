use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::hv::{hash_hv_bit, hash_hv_i8};
use rand::rngs::{SmallRng, StdRng};
use rand::{Rng, RngCore, SeedableRng};
use rapidhash::{RapidHashSet, RapidRng};

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
