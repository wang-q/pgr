//! Benchmark: IntSpan set operations and bulk construction.
//!
//! Covers the linearized set operations (`intersect`/`union`/`diff`/`xor`,
//! rewritten as two-pointer merges) and the sorted bulk builder
//! (`from_pairs`) against their previous implementations, on large synthetic
//! runlists. The old implementations are reconstructed here from the public
//! API (`complement`/`merge`/`subtract`/`invert`/`add_pair`).

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pgr::libs::ds::IntSpan;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

const SEED: u64 = 20260804;
const CHR_LEN: i32 = 100_000_000;

/// Generate a runlist of `n` random spans (length 100..2000) over a 100 Mb
/// chromosome, built with the sorted bulk builder.
fn random_runlist(n: usize, rng: &mut StdRng) -> IntSpan {
    let pairs: Vec<(i32, i32)> = (0..n)
        .map(|_| {
            let start = rng.random_range(0..CHR_LEN - 2000);
            (start, start + rng.random_range(100..=2000) - 1)
        })
        .collect();
    IntSpan::from_pairs(pairs)
}

// ── Previous implementations, kept as baselines ──────────────────────────

fn old_intersect(a: &IntSpan, b: &IntSpan) -> IntSpan {
    if a.is_empty() || b.is_empty() {
        IntSpan::new()
    } else {
        let mut new = a.complement();
        new.merge(&b.complement());
        new.invert();
        new
    }
}

fn old_union(a: &IntSpan, b: &IntSpan) -> IntSpan {
    let mut new = a.copy();
    new.merge(b);
    new
}

fn old_diff(a: &IntSpan, b: &IntSpan) -> IntSpan {
    if a.is_empty() {
        IntSpan::new()
    } else {
        let mut new = a.copy();
        new.subtract(b);
        new
    }
}

fn old_xor(a: &IntSpan, b: &IntSpan) -> IntSpan {
    let mut new = old_union(a, b);
    new.subtract(&old_intersect(a, b));
    new
}

fn bench_setops(c: &mut Criterion) {
    let mut group = c.benchmark_group("setops");
    for &n in &[5_000usize, 20_000] {
        let mut rng = StdRng::seed_from_u64(SEED ^ (n as u64));
        let a = random_runlist(n, &mut rng);
        let b = random_runlist(n, &mut rng);
        for (name, f) in [
            (
                "intersect",
                old_intersect as fn(&IntSpan, &IntSpan) -> IntSpan,
            ),
            ("union", old_union),
            ("diff", old_diff),
            ("xor", old_xor),
        ] {
            group.bench_with_input(
                BenchmarkId::new(format!("{name} old"), n),
                &(&a, &b),
                |bb, (x, y)| bb.iter(|| black_box(f(x, y))),
            );
        }
        group.bench_with_input(
            BenchmarkId::new("intersect new", n),
            &(&a, &b),
            |bb, (x, y)| bb.iter(|| black_box(x.intersect(y))),
        );
        group.bench_with_input(BenchmarkId::new("union new", n), &(&a, &b), |bb, (x, y)| {
            bb.iter(|| black_box(x.union(y)))
        });
        group.bench_with_input(BenchmarkId::new("diff new", n), &(&a, &b), |bb, (x, y)| {
            bb.iter(|| black_box(x.diff(y)))
        });
        group.bench_with_input(BenchmarkId::new("xor new", n), &(&a, &b), |bb, (x, y)| {
            bb.iter(|| black_box(x.xor(y)))
        });
    }
    group.finish();
}

fn bench_construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("construction");
    for &n in &[10_000usize, 100_000] {
        // Adversarial case for incremental insertion: 1 bp spans at random
        // positions, mostly disjoint (each insert lands mid-VecDeque).
        let mut rng = StdRng::seed_from_u64(SEED ^ (n as u64) << 8);
        let pairs: Vec<(i32, i32)> = (0..n)
            .map(|_| {
                let p = rng.random_range(0..CHR_LEN);
                (p, p)
            })
            .collect();
        group.bench_with_input(BenchmarkId::new("from_pairs", n), &pairs, |bb, p| {
            bb.iter(|| black_box(IntSpan::from_pairs(p.clone())))
        });
        group.bench_with_input(BenchmarkId::new("add_pair loop", n), &pairs, |bb, p| {
            bb.iter(|| {
                let mut s = IntSpan::new();
                for &(l, u) in p {
                    s.add_pair(l, u);
                }
                black_box(s)
            })
        });
    }
    group.finish();
}

fn bench_covered(c: &mut Criterion) {
    let mut group = c.benchmark_group("covered");
    // The `intersect+cardinality` baseline is O(n) per query, so keep both
    // n and the query count small enough for criterion to finish quickly.
    for &n in &[2_000usize, 5_000] {
        let mut rng = StdRng::seed_from_u64(SEED ^ (n as u64) << 4);
        let set = random_runlist(n, &mut rng);
        let queries: Vec<(i32, i32)> = (0..2_000)
            .map(|_| {
                let s = rng.random_range(0..CHR_LEN - 2000);
                (s, s + rng.random_range(100..=2000) - 1)
            })
            .collect();
        // SpanIndex-style: the same spans as a flat Vec, queried with
        // `partition_point` (the structure we merged back into IntSpan).
        let spans: Vec<(i32, i32)> = set.spans();

        group.bench_with_input(
            BenchmarkId::new("covered", n),
            &(&set, &queries),
            |bb, (s, q)| {
                bb.iter(|| {
                    let mut total = 0i64;
                    for &(a, b) in q.iter() {
                        total += i64::from(s.covered(a, b));
                    }
                    black_box(total)
                })
            },
        );
        group.bench_with_input(
            BenchmarkId::new("intersect+cardinality", n),
            &(&set, &queries),
            |bb, (s, q)| {
                bb.iter(|| {
                    let mut total = 0i64;
                    for &(a, b) in q.iter() {
                        let mut r = IntSpan::new();
                        r.add_pair(a, b);
                        total += i64::from(s.intersect(&r).cardinality());
                    }
                    black_box(total)
                })
            },
        );
        group.bench_with_input(
            BenchmarkId::new("partition_point vec", n),
            &(&spans, &queries),
            |bb, (sp, q)| {
                bb.iter(|| {
                    let mut total = 0i64;
                    for &(a, b) in q.iter() {
                        let first = sp.partition_point(|&(_, e)| e < a);
                        let last = sp.partition_point(|&(st, _)| st <= b);
                        for &(st, e) in &sp[first..last] {
                            total += i64::from(b.min(e) - a.max(st) + 1);
                        }
                    }
                    black_box(total)
                })
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_setops, bench_construction, bench_covered);
criterion_main!(benches);
