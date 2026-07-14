//! Benchmark: BitMap vs IntSpan for dense range coverage on a single chromosome.
//!
//! Tests the claim in `notes/chain-algorithm-reuse.md` that BitMap's bit-vector
//! layout is more compact and potentially faster than IntSpan's interval list
//! for large chromosome coverage aggregation.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use intspan::IntSpan;
use pgr::libs::chain::BitMap;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

const SEED: u64 = 20260715;

/// A prepared workload for one chromosome size.
struct Workload {
    intervals: Vec<(u64, u64)>,
    query_intervals: Vec<(u64, u64)>,
}

/// Generate non-overlapping-ish random intervals covering roughly `coverage_frac`
/// of the chromosome. Each interval has length uniformly in [min_len, max_len].
fn generate_workload(
    chr_size: u64,
    coverage_frac: f64,
    min_len: u64,
    max_len: u64,
    rng: &mut StdRng,
) -> Workload {
    let mut intervals = Vec::new();
    let mut covered = 0u64;
    let target = (chr_size as f64 * coverage_frac) as u64;

    while covered < target {
        let len = rng.random_range(min_len..=max_len);
        let start = rng.random_range(0..chr_size.saturating_sub(len));
        intervals.push((start, start + len));
        covered += len;
    }

    // Pick a deterministic subset of intervals as query targets.
    let n_queries = intervals.len().min(10_000);
    let query_intervals: Vec<_> = intervals.iter().take(n_queries).copied().collect();

    Workload {
        intervals,
        query_intervals,
    }
}

fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("build");

    for chr_size in [1_000_000u64, 10_000_000u64, 100_000_000u64] {
        let workload = generate_workload(
            chr_size,
            0.10,
            100,
            10_000,
            &mut StdRng::seed_from_u64(SEED),
        );
        let intervals = workload.intervals.clone();

        group.bench_with_input(BenchmarkId::new("BitMap", chr_size), &chr_size, |b, _| {
            b.iter(|| {
                let mut bm = BitMap::new(black_box(chr_size));
                for (start, end) in &intervals {
                    bm.set_range(black_box(*start), black_box(*end - *start));
                }
                black_box(bm);
            })
        });

        group.bench_with_input(BenchmarkId::new("IntSpan", chr_size), &chr_size, |b, _| {
            b.iter(|| {
                let mut ints = IntSpan::new();
                for (start, end) in &intervals {
                    // IntSpan uses 1-based closed intervals.
                    ints.add_pair(black_box((*start + 1) as i32), black_box(*end as i32));
                }
                black_box(ints);
            })
        });
    }

    group.finish();
}

fn bench_query(c: &mut Criterion) {
    let mut group = c.benchmark_group("query");

    for chr_size in [1_000_000u64, 10_000_000u64, 100_000_000u64] {
        let workload = generate_workload(
            chr_size,
            0.10,
            100,
            10_000,
            &mut StdRng::seed_from_u64(SEED),
        );

        // Pre-build containers.
        let mut bm = BitMap::new(chr_size);
        let mut ints = IntSpan::new();
        for (start, end) in &workload.intervals {
            bm.set_range(*start, *end - *start);
            ints.add_pair((*start + 1) as i32, *end as i32);
        }

        // Pre-build IntSpan query objects for superset checks.
        let int_queries: Vec<IntSpan> = workload
            .query_intervals
            .iter()
            .map(|(start, end)| IntSpan::from_pair((*start + 1) as i32, *end as i32))
            .collect();

        group.bench_with_input(
            BenchmarkId::new("BitMap_is_fully_set", chr_size),
            &chr_size,
            |b, _| {
                b.iter(|| {
                    let mut count = 0usize;
                    for (start, end) in &workload.query_intervals {
                        if bm.is_fully_set(black_box(*start), black_box(*end - *start)) {
                            count += 1;
                        }
                    }
                    black_box(count);
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("IntSpan_superset", chr_size),
            &chr_size,
            |b, _| {
                b.iter(|| {
                    let mut count = 0usize;
                    for query in &int_queries {
                        if ints.superset(black_box(query)) {
                            count += 1;
                        }
                    }
                    black_box(count);
                })
            },
        );
    }

    group.finish();
}

fn bench_mixed(c: &mut Criterion) {
    let mut group = c.benchmark_group("mixed_build_then_query");

    for chr_size in [1_000_000u64, 10_000_000u64, 100_000_000u64] {
        let workload = generate_workload(
            chr_size,
            0.10,
            100,
            10_000,
            &mut StdRng::seed_from_u64(SEED),
        );
        let intervals = workload.intervals.clone();
        let queries = workload.query_intervals.clone();
        let int_queries: Vec<IntSpan> = queries
            .iter()
            .map(|(start, end)| IntSpan::from_pair((*start + 1) as i32, *end as i32))
            .collect();

        group.bench_with_input(BenchmarkId::new("BitMap", chr_size), &chr_size, |b, _| {
            b.iter(|| {
                let mut bm = BitMap::new(black_box(chr_size));
                for (start, end) in &intervals {
                    bm.set_range(*start, *end - *start);
                }
                let mut count = 0usize;
                for (start, end) in &queries {
                    if bm.is_fully_set(*start, *end - *start) {
                        count += 1;
                    }
                }
                black_box(count);
            })
        });

        group.bench_with_input(BenchmarkId::new("IntSpan", chr_size), &chr_size, |b, _| {
            b.iter(|| {
                let mut ints = IntSpan::new();
                for (start, end) in &intervals {
                    ints.add_pair((*start + 1) as i32, *end as i32);
                }
                let mut count = 0usize;
                for query in &int_queries {
                    if ints.superset(query) {
                        count += 1;
                    }
                }
                black_box(count);
            })
        });
    }

    group.finish();

    // One-time memory comparison for the largest workload.
    print_memory_estimates();
}

/// Print memory estimates for the largest workload so the human reader can
/// compare compactness without relying on Criterion's time measurements.
fn print_memory_estimates() {
    let workload = generate_workload(
        100_000_000,
        0.10,
        100,
        10_000,
        &mut StdRng::seed_from_u64(SEED),
    );

    let mut bm = BitMap::new(100_000_000);
    let mut ints = IntSpan::new();
    for (start, end) in &workload.intervals {
        bm.set_range(*start, *end - *start);
        ints.add_pair((*start + 1) as i32, *end as i32);
    }

    let bitmap_bytes = bm.memory_size();
    let intspan_internal_len = ints.to_vec().len();
    let intspan_bytes_estimate = intspan_internal_len * std::mem::size_of::<i32>();

    eprintln!(
        "Memory estimate (chr_size=100M, ~10% coverage, {} intervals): BitMap={} bytes, IntSpan~{} bytes ({} internal ints)",
        workload.intervals.len(),
        bitmap_bytes,
        intspan_bytes_estimate,
        intspan_internal_len,
    );
}

criterion_group!(benches, bench_build, bench_query, bench_mixed);
criterion_main!(benches);
