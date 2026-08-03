//! Benchmark: point-membership and overlap queries across the three interval
//! structures available in pgr.
//!
//! * `libs::ds::IntSpan` — merged runlist set: answers "is the point covered"
//!   via binary search, but discards per-interval identity.
//! * `coitrees::BasicCOITree` — COITree interval index (already used by the
//!   PAF/pbit indexes): enumerates all overlapping intervals.
//! * `rust-lapper` — sorted-start interval index (used by the external
//!   intspan project): enumerates overlapping intervals and counts via BITS.
//!
//! Data generation mirrors `rust-lapper-master/benches/lapper_benchmark.rs`
//! (random intervals on a 100 Mb chromosome, lengths 500..80 kb), plus a
//! chromosome-spanning interval that degrades sorted-start indexes.

use coitrees::{BasicCOITree, Interval as CoiInterval, IntervalTree};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use pgr::libs::ds::IntSpan;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rust_lapper::{Interval as LapInterval, Lapper};

const CHROM_SIZE: u32 = 100_000_000;
const MIN_LEN: u32 = 500;
const MAX_LEN: u32 = 80_000;
const SEED: u64 = 20260804;

/// A prepared workload: intervals (half-open `[start, end)`) plus query points
/// and query windows.
struct Workload {
    intervals: Vec<(u32, u32)>,
    points: Vec<u32>,
    windows: Vec<(u32, u32)>,
}

fn random_intervals(n: usize, rng: &mut StdRng) -> Vec<(u32, u32)> {
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let start = rng.random_range(0..CHROM_SIZE);
        let len = rng.random_range(MIN_LEN..=MAX_LEN);
        out.push((start, start.saturating_add(len).min(CHROM_SIZE)));
    }
    out
}

fn make_workload(n: usize, spanning: bool) -> Workload {
    let mut rng = StdRng::seed_from_u64(SEED ^ (n as u64) << 32 ^ u64::from(spanning));
    let mut intervals = random_intervals(n, &mut rng);
    if spanning {
        // One interval engulfing most of the chromosome: the pathological
        // case for sorted-start indexes (rust-lapper), which coitrees handles
        // with guaranteed bounds.
        intervals.push((0, CHROM_SIZE * 9 / 10));
    }
    let points = (0..n).map(|_| rng.random_range(0..CHROM_SIZE)).collect();
    let windows = (0..n)
        .map(|_| {
            let start = rng.random_range(0..CHROM_SIZE);
            (start, start.saturating_add(2_000).min(CHROM_SIZE))
        })
        .collect();
    Workload {
        intervals,
        points,
        windows,
    }
}

/// Build the three structures from the same half-open interval list.
fn build_structures(
    intervals: &[(u32, u32)],
) -> (IntSpan, BasicCOITree<bool, u32>, Lapper<u32, bool>) {
    // IntSpan: merged inclusive runlist (set semantics).
    let mut intspan = IntSpan::new();
    for &(s, e) in intervals {
        intspan.add_pair(s as i32, e.saturating_sub(1) as i32);
    }

    // coitrees: half-open `[start, end)` -> inclusive `[first, last]`.
    let coi: Vec<CoiInterval<bool>> = intervals
        .iter()
        .map(|&(s, e)| CoiInterval::new(s as i32, e.saturating_sub(1) as i32, true))
        .collect();
    let coitree = BasicCOITree::new(&coi);

    // rust-lapper: half-open `[start, stop)`.
    let lap: Vec<LapInterval<u32, bool>> = intervals
        .iter()
        .map(|&(s, e)| LapInterval {
            start: s,
            stop: e,
            val: true,
        })
        .collect();
    let lapper = Lapper::new(lap);

    (intspan, coitree, lapper)
}

fn construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("construction");
    for &n in &[1_000usize, 10_000, 100_000] {
        group.throughput(Throughput::Elements(n as u64));
        let intervals = random_intervals(n, &mut StdRng::seed_from_u64(n as u64));
        group.bench_with_input(
            BenchmarkId::new("intspan add_pair", n),
            &intervals,
            |b, iv| {
                b.iter(|| {
                    let mut s = IntSpan::new();
                    for &(lo, hi) in black_box(iv) {
                        s.add_pair(lo as i32, hi.saturating_sub(1) as i32);
                    }
                    black_box(s);
                })
            },
        );
        let coi: Vec<CoiInterval<bool>> = intervals
            .iter()
            .map(|&(s, e)| CoiInterval::new(s as i32, e.saturating_sub(1) as i32, true))
            .collect();
        group.bench_with_input(BenchmarkId::new("coitrees new", n), &coi, |b, iv| {
            b.iter(|| {
                let tree: BasicCOITree<bool, u32> = BasicCOITree::new(black_box(iv));
                black_box(tree)
            })
        });
        let lap: Vec<LapInterval<u32, bool>> = intervals
            .iter()
            .map(|&(s, e)| LapInterval {
                start: s,
                stop: e,
                val: true,
            })
            .collect();
        group.bench_with_input(BenchmarkId::new("lapper new", n), &lap, |b, iv| {
            b.iter(|| black_box(Lapper::new(black_box(iv.clone()))))
        });
    }
    group.finish();
}

fn point_membership(c: &mut Criterion) {
    let mut group = c.benchmark_group("point membership");
    for &n in &[1_000usize, 10_000, 100_000] {
        for spanning in [false, true] {
            let wl = make_workload(n, spanning);
            let (intspan, coitree, lapper) = build_structures(&wl.intervals);
            let suffix = if spanning { "spanning" } else { "normal" };
            group.bench_with_input(
                BenchmarkId::new(format!("intspan contains {suffix}"), n),
                &wl.points,
                |b, pts| {
                    b.iter(|| {
                        let mut hits = 0u64;
                        for &p in black_box(pts) {
                            hits += u64::from(intspan.contains(p as i32));
                        }
                        black_box(hits);
                    })
                },
            );
            group.bench_with_input(
                BenchmarkId::new(format!("coitree point query {suffix}"), n),
                &wl.points,
                |b, pts| {
                    b.iter(|| {
                        let mut hits = 0u64;
                        for &p in black_box(pts) {
                            let mut k = 0usize;
                            coitree.query(p as i32, p as i32, |_| k += 1);
                            hits += u64::from(k > 0);
                        }
                        black_box(hits);
                    })
                },
            );
            group.bench_with_input(
                BenchmarkId::new(format!("lapper count {suffix}"), n),
                &wl.points,
                |b, pts| {
                    b.iter(|| {
                        let mut hits = 0u64;
                        for &p in black_box(pts) {
                            hits += u64::from(lapper.count(p, p + 1) > 0);
                        }
                        black_box(hits);
                    })
                },
            );
        }
    }
    group.finish();
}

fn overlap_enumeration(c: &mut Criterion) {
    let mut group = c.benchmark_group("overlap enumeration");
    for &n in &[1_000usize, 10_000, 100_000] {
        for spanning in [false, true] {
            let wl = make_workload(n, spanning);
            let (_, coitree, lapper) = build_structures(&wl.intervals);
            let suffix = if spanning { "spanning" } else { "normal" };
            group.bench_with_input(
                BenchmarkId::new(format!("coitree query {suffix}"), n),
                &wl.windows,
                |b, wins| {
                    b.iter(|| {
                        let mut hits = 0usize;
                        for &(s, e) in black_box(wins) {
                            let mut k = 0usize;
                            coitree.query(s as i32, e.saturating_sub(1) as i32, |_| k += 1);
                            hits += k;
                        }
                        black_box(hits);
                    })
                },
            );
            group.bench_with_input(
                BenchmarkId::new(format!("lapper find {suffix}"), n),
                &wl.windows,
                |b, wins| {
                    b.iter(|| {
                        let mut hits = 0usize;
                        for &(s, e) in black_box(wins) {
                            hits += lapper.find(s, e).count();
                        }
                        black_box(hits);
                    })
                },
            );
        }
    }
    group.finish();
}

criterion_group!(benches, construction, point_membership, overlap_enumeration);
criterion_main!(benches);
