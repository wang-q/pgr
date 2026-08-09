use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::hv::{hash_hv_bit, hash_hv_i8, hash_hv_sparse};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rapidhash::RapidHashSet;

// Production-path HV benchmarks: `hash_hv_bit` / `hash_hv_i8` (dense) and
// `hash_hv_sparse` across seed-set sizes and dimensions. Historical
// reference implementations (AVX-512, RNG candidates, i16/pshufb, sampling
// hashes) live in `hv_benchmark_ref.rs`.

fn generate_kmer_hash_set(size: usize) -> RapidHashSet<u64> {
    let mut rng = StdRng::seed_from_u64(42); // Fixed seed for reproducibility
    let mut kmer_hash_set = RapidHashSet::default();

    for _ in 0..size {
        kmer_hash_set.insert(rng.random::<u64>());
    }

    kmer_hash_set
}

fn bench_encode_hash_hd(c: &mut Criterion) {
    let kmer_hash_set_small = generate_kmer_hash_set(1000); // Small dataset
    let kmer_hash_set_medium = generate_kmer_hash_set(10_000); // Medium dataset
    let kmer_hash_set_large = generate_kmer_hash_set(100_000); // Large dataset

    let seed_vec_small: Vec<u64> = kmer_hash_set_small.iter().cloned().collect();
    let seed_vec_medium: Vec<u64> = kmer_hash_set_medium.iter().cloned().collect();
    let seed_vec_large: Vec<u64> = kmer_hash_set_large.iter().cloned().collect();

    let hv_d = 4096; // Hypervector dimension

    c.bench_function("encode_hash_hd_lib_small", |b| {
        b.iter(|| hash_hv_bit(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_simd_i8_small", |b| {
        b.iter(|| hash_hv_i8(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_lib_medium", |b| {
        b.iter(|| hash_hv_bit(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_simd_i8_medium", |b| {
        b.iter(|| hash_hv_i8(black_box(&seed_vec_medium), hv_d))
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
    // Encoding time vs D at fixed s: sparse cost is O(n·s), independent of D.
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

    // D = 16384 variants on the medium (10k) seed set
    let hv_d_16k = 16384;
    c.bench_function("encode_hash_hd_lib_d16k", |b| {
        b.iter(|| hash_hv_bit(black_box(&seed_vec_medium), hv_d_16k))
    });
    c.bench_function("encode_hash_hd_simd_i8_d16k", |b| {
        b.iter(|| hash_hv_i8(black_box(&seed_vec_medium), hv_d_16k))
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10); // Set sample size
    targets = bench_encode_hash_hd
);
criterion_main!(benches);
