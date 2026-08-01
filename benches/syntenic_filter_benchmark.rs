//! Benchmark: linear-scan vs DupeTree for syntenic-filter overlap checks.
//!
//! `pgr paf query --syntenic-filter` originally checked each query interval
//! against every chain span linearly; the DupeTree rewrite turns that into a
//! tree query (see notes/chain-algorithm-reuse.md).

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pgr::libs::ds::DupeTree;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

const SEED: u64 = 20260801;

fn bench_linear(spans: &[(u64, u64)], queries: &[(u64, u64)]) -> usize {
    queries
        .iter()
        .filter(|&&(qs, qe)| spans.iter().any(|&(cs, ce)| qs < ce && qe > cs))
        .count()
}

fn bench_dupetree(spans: &[(u64, u64)], queries: &[(u64, u64)]) -> usize {
    let mut tree = DupeTree::new();
    for &(s, e) in spans {
        tree.add(s, e);
    }
    tree.build();
    queries
        .iter()
        .filter(|&&(qs, qe)| tree.count_over(qs, qe, 1) > 0)
        .count()
}

fn bench_filter(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(SEED);
    let mut group = c.benchmark_group("syntenic_filter");
    for span_count in [1_000usize, 10_000, 50_000] {
        // Non-overlapping chain spans covering a chromosome-like coordinate line.
        let mut spans = Vec::with_capacity(span_count);
        let mut pos = 0u64;
        for _ in 0..span_count {
            let len = rng.random_range(100..50_000);
            pos += rng.random_range(0..2_000);
            spans.push((pos, pos + len));
        }
        // Query intervals anchored inside random spans, so most overlap.
        let queries: Vec<(u64, u64)> = (0..span_count)
            .map(|_| {
                let (s, _) = spans[rng.random_range(0..spans.len())];
                let qs = s + rng.random_range(0..200);
                (qs, qs + rng.random_range(20..100))
            })
            .collect();

        group.bench_with_input(
            BenchmarkId::new("linear", span_count),
            &span_count,
            |b, _| b.iter(|| bench_linear(black_box(&spans), black_box(&queries))),
        );
        group.bench_with_input(
            BenchmarkId::new("dupetree", span_count),
            &span_count,
            |b, _| b.iter(|| bench_dupetree(black_box(&spans), black_box(&queries))),
        );
    }
    group.finish();
}

criterion_group!(benches, bench_filter);
criterion_main!(benches);
