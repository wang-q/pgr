//! HNSW vs exact linear scan for 4096-d HV: recall@10 and per-query latency.
//!
//! Synthetic cohort: each genome = shared k-mer core + private k-mers,
//! encoded with `hash_hv_bit` (the real HV pipeline), then L2-normalized so
//! Euclidean distance ranks identically to cosine similarity (the ordering
//! pgr uses for HV near-neighbor search; see notes/design/genome-nn-query.md
//! §6.4). Ground truth is exact cosine top-k via `linalg::dot_product`.

use std::time::{Duration, Instant};

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use hnsw::{Hnsw, Params, Searcher};
use pgr::libs::hv::hash_hv_bit;
use pgr::libs::linalg::dot_product;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use space::{Metric, Neighbor};

const D: usize = 4096;
const K: usize = 10;
const QUERIES: usize = 50;
const EF_GRID: [usize; 6] = [10, 20, 50, 100, 200, 400];

// Each genome shares a fixed k-mer core and draws its private k-mers from a
// per-genome stream, so pairwise HV similarity is controlled by the
// core/private ratio (mimics related genomes within one species).
const CORE_HASHES: usize = 2048;
const PRIVATE_MIN: usize = 512;
const PRIVATE_MAX: usize = 4096;

/// Euclidean distance over L2-normalized vectors; ordering-equivalent to cosine.
struct CosineL2;

impl Metric<&'static [f32]> for CosineL2 {
    type Unit = u64;

    fn distance(&self, a: &&'static [f32], b: &&'static [f32]) -> u64 {
        // f64 accumulation keeps the metric stable at D=4096 (the crate docs
        // warn that naive f32 sums can break the triangle inequality).
        let mut acc = 0.0f64;
        for (&x, &y) in a.iter().zip(b.iter()) {
            let d = (x - y) as f64;
            acc += d * d;
        }
        acc.sqrt().to_bits()
    }
}

type HvHnsw = Hnsw<CosineL2, &'static [f32], rand_pcg::Pcg64, 12, 24>;

fn l2_normalize(v: &[i32]) -> &'static [f32] {
    let mut acc = 0.0f64;
    for &x in v {
        acc += (x as f64) * (x as f64);
    }
    let norm = acc.sqrt() as f32;
    let out: Vec<f32> = v.iter().map(|&x| x as f32 / norm).collect();
    Box::leak(out.into_boxed_slice())
}

fn genome_hv(core: &[u64], rng: &mut SmallRng) -> &'static [f32] {
    let private_len = PRIVATE_MIN + rng.random_range(0..=PRIVATE_MAX - PRIVATE_MIN);
    let mut hashes = Vec::with_capacity(CORE_HASHES + private_len);
    hashes.extend_from_slice(core);
    hashes.extend((0..private_len).map(|_| rng.random::<u64>()));
    l2_normalize(&hash_hv_bit(&hashes, D))
}

fn generate_cohort(seed: u64, n: usize) -> (Vec<&'static [f32]>, Vec<&'static [f32]>) {
    let mut core_rng = SmallRng::seed_from_u64(seed ^ 0xC0FFEE);
    let core: Vec<u64> = (0..CORE_HASHES).map(|_| core_rng.random()).collect();

    let mut db = Vec::with_capacity(n);
    for i in 0..n {
        let mut rng = SmallRng::seed_from_u64(seed ^ (i as u64 + 1));
        db.push(genome_hv(&core, &mut rng));
    }

    let mut queries = Vec::with_capacity(QUERIES);
    for j in 0..QUERIES {
        let mut rng = SmallRng::seed_from_u64(seed ^ (1_000_003 + j as u64));
        queries.push(genome_hv(&core, &mut rng));
    }
    (db, queries)
}

fn exact_topk(db: &[&'static [f32]], query: &[f32], k: usize) -> Vec<usize> {
    let mut scored: Vec<(f32, usize)> = db
        .iter()
        .enumerate()
        .map(|(i, v)| (dot_product(query, v), i))
        .collect();
    scored.sort_unstable_by(|a, b| b.0.total_cmp(&a.0));
    scored[..k].iter().map(|&(_, i)| i).collect()
}

fn hnsw_topk(
    hnsw: &HvHnsw,
    searcher: &mut Searcher<u64>,
    query: &'static [f32],
    ef: usize,
) -> Vec<usize> {
    let mut dest = vec![
        Neighbor {
            index: usize::MAX,
            distance: u64::MAX,
        };
        K
    ];
    let found = hnsw.nearest(&query, ef, searcher, &mut dest);
    found.iter().map(|n| n.index).collect()
}

fn build_hnsw(db: &[&'static [f32]]) -> (HvHnsw, Searcher<u64>) {
    let ef_c = ef_construction();
    let mut hnsw = HvHnsw::new_params(CosineL2, Params::new().ef_construction(ef_c));
    let mut searcher = Searcher::default();
    let t0 = Instant::now();
    for &point in db {
        hnsw.insert(point, &mut searcher);
    }
    eprintln!(
        "    build hnsw N={} ef_c={}: {:?}",
        db.len(),
        ef_c,
        t0.elapsed()
    );
    (hnsw, searcher)
}

fn parse_sizes() -> Vec<usize> {
    let parsed: Vec<usize> = std::env::var("PGR_HV_ANN_SIZES")
        .map(|s| s.split(',').filter_map(|x| x.trim().parse().ok()).collect())
        .unwrap_or_default();
    if parsed.is_empty() {
        vec![1000, 10_000, 30_000]
    } else {
        parsed
    }
}

fn ef_construction() -> usize {
    std::env::var("PGR_HV_ANN_EFC")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(64)
}

fn hv_ann_recall(c: &mut Criterion) {
    for &n in &parse_sizes() {
        let (db, queries) = generate_cohort(0x5EED_2026, n);
        let exact: Vec<Vec<usize>> = queries.iter().map(|q| exact_topk(&db, q, K)).collect();
        let (hnsw, mut searcher) = build_hnsw(&db);

        eprintln!("== N={n} ==");
        for &ef in &EF_GRID {
            let mut hits = 0usize;
            let mut total = Duration::ZERO;
            for (qi, q) in queries.iter().enumerate() {
                let t0 = Instant::now();
                let ann = hnsw_topk(&hnsw, &mut searcher, q, ef);
                total += t0.elapsed();
                hits += ann.iter().filter(|i| exact[qi].contains(i)).count();
            }
            let recall = hits as f64 / (QUERIES as f64 * K as f64);
            eprintln!(
                "    ef={ef:>3}  recall@10={recall:.3}  avg_query={:.2} us",
                total.as_secs_f64() * 1e6 / QUERIES as f64
            );

            let mut group = c.benchmark_group(format!("hnsw/{n}"));
            group
                .sample_size(10)
                .warm_up_time(Duration::from_millis(500))
                .measurement_time(Duration::from_secs(2));
            group.throughput(Throughput::Elements(QUERIES as u64));
            group.bench_function(format!("ef{ef}"), |b| {
                let mut dest = vec![
                    Neighbor {
                        index: usize::MAX,
                        distance: u64::MAX,
                    };
                    K
                ];
                b.iter(|| {
                    for q in &queries {
                        black_box(hnsw.nearest(q, ef, &mut searcher, &mut dest));
                    }
                });
            });
            group.finish();
        }

        let mut group = c.benchmark_group(format!("exact/{n}"));
        group
            .sample_size(10)
            .warm_up_time(Duration::from_millis(500))
            .measurement_time(Duration::from_secs(2));
        group.throughput(Throughput::Elements(QUERIES as u64));
        group.bench_function("scan", |b| {
            b.iter(|| {
                for q in &queries {
                    black_box(exact_topk(&db, q, K));
                }
            });
        });
        group.finish();
    }
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = hv_ann_recall
);

criterion_main!(benches);
