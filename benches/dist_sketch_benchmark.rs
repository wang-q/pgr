//! Throughput comparison of the three sketch samplers used by `pgr dist`
//! (`mini` minimizer, `mash` bottom-k MinHash, `frac` FracMinHash) on a
//! 1 Mb random DNA sequence. Mash and FracMinHash hash every k-mer; the
//! minimizer is windowed. See notes/benchmarks/bench-simd-hv-jaccard.md
//! for the broader HV/RNG benchmark context.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use pgr::libs::hash::{for_each_mash_hash, seq_fracminhash, seq_mins, BottomK};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

const SEQ_LEN: usize = 1_000_000;
const SEED: u64 = 20260808;

fn random_dna(len: usize) -> Vec<u8> {
    let mut rng = StdRng::seed_from_u64(SEED);
    (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
}

fn bench_mini(c: &mut Criterion) {
    let seq = random_dna(SEQ_LEN);
    let mut group = c.benchmark_group("dist_sketch_mini");
    group.throughput(Throughput::Bytes(SEQ_LEN as u64));
    group.bench_function("k21_w5", |b| {
        b.iter(|| black_box(seq_mins(&seq, "rapid", 21, 5).unwrap().len()))
    });
    group.finish();
}

fn bench_mash(c: &mut Criterion) {
    let seq = random_dna(SEQ_LEN);
    let mut group = c.benchmark_group("dist_sketch_mash");
    group.throughput(Throughput::Bytes(SEQ_LEN as u64));
    group.bench_function("k21_size1000", |b| {
        b.iter(|| {
            // The production path: streaming rolling-window hashes into a
            // bounded bottom-k accumulator (no full-length hash materialization).
            let mut acc = BottomK::new(1000);
            for_each_mash_hash(&seq, 21, 42, |h| acc.insert(h));
            black_box(acc.into_set().len())
        })
    });
    group.finish();
}

fn bench_frac(c: &mut Criterion) {
    let seq = random_dna(SEQ_LEN);
    let mut group = c.benchmark_group("dist_sketch_frac");
    group.throughput(Throughput::Bytes(SEQ_LEN as u64));
    group.bench_function("k21_scale1000", |b| {
        b.iter(|| black_box(seq_fracminhash(&seq, 21, 1000, false).unwrap().len()))
    });
    group.finish();
}

criterion_group!(benches, bench_mini, bench_mash, bench_frac);
criterion_main!(benches);
