//! hnsw_rs HNSW vs single-layer HubNSW for 4096-d HV: recall@10 and latency.
//!
//! Same synthetic cohort and exact baseline as `hv_ann_recall.rs` (shared
//! k-mer core + private k-mers -> `hash_hv_bit` -> L2-normalized f32).
//! hnsw_rs supports GSearch's HubNSW trick: `modify_level_scale(0.2)` makes
//! level sampling collapse to layer 0 (flat NSW), which the GSearch README
//! recommends for high-dimensional data (arXiv 2412.01940).

use std::time::{Duration, Instant};

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use hnsw_rs::prelude::*;
use pgr::libs::hv::hash_hv_bit;
use pgr::libs::linalg::dot_product;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

const D: usize = 4096;
const K: usize = 10;
const QUERIES: usize = 50;
const EF_GRID: [usize; 6] = [10, 20, 50, 100, 200, 400];
const SCALES: [f64; 2] = [1.0, 0.2];
const M_CONN: usize = 24;

const CORE_HASHES: usize = 2048;
const PRIVATE_MIN: usize = 512;
const PRIVATE_MAX: usize = 4096;

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

fn build_hnsw_rs(db: &[&'static [f32]], scale: f64) -> Hnsw<'static, f32, DistL2> {
    let ef_c = ef_construction();
    let mut hnsw = Hnsw::new(M_CONN, db.len(), 16, ef_c, DistL2);
    hnsw.modify_level_scale(scale);
    let t0 = Instant::now();
    for (i, &point) in db.iter().enumerate() {
        hnsw.insert((point, i));
    }
    eprintln!(
        "    build hnsw_rs N={} ef_c={} scale={scale}: {:?}",
        db.len(),
        ef_c,
        t0.elapsed()
    );
    hnsw
}

fn hnsw_rs_topk(hnsw: &Hnsw<'static, f32, DistL2>, query: &[f32], ef: usize) -> Vec<usize> {
    let mut nbs = hnsw.search(query, K, ef);
    nbs.sort_unstable_by(|a, b| a.distance.total_cmp(&b.distance));
    nbs.iter().map(|n| n.d_id).collect()
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

fn parse_efs() -> Vec<usize> {
    let parsed: Vec<usize> = std::env::var("PGR_HV_ANN_EFS")
        .map(|s| s.split(',').filter_map(|x| x.trim().parse().ok()).collect())
        .unwrap_or_default();
    if parsed.is_empty() {
        EF_GRID.to_vec()
    } else {
        parsed
    }
}

fn hv_ann_hubnsw(c: &mut Criterion) {
    for &n in &parse_sizes() {
        let (db, queries) = generate_cohort(0x5EED_2026, n);
        let exact: Vec<Vec<usize>> = queries.iter().map(|q| exact_topk(&db, q, K)).collect();
        let efs = parse_efs();
        let mut summary = Vec::<(f64, usize, f64, f64)>::new();

        for &scale in &SCALES {
            let hnsw = build_hnsw_rs(&db, scale);
            eprintln!("== N={n} scale={scale} ==");
            for &ef in &efs {
                let mut hits = 0usize;
                let mut total = Duration::ZERO;
                for (qi, q) in queries.iter().enumerate() {
                    let t0 = Instant::now();
                    let ann = hnsw_rs_topk(&hnsw, q, ef);
                    total += t0.elapsed();
                    hits += ann.iter().filter(|i| exact[qi].contains(i)).count();
                }
                let recall = hits as f64 / (QUERIES as f64 * K as f64);
                let latency_us = total.as_secs_f64() * 1e6 / QUERIES as f64;
                summary.push((scale, ef, recall, latency_us));
                eprintln!(
                    "    ef={ef:>3}  recall@10={recall:.3}  avg_query={:.2} us",
                    latency_us
                );

                let mut group = c.benchmark_group(format!("hnsw_rs/{n}/scale{scale:.1}"));
                group
                    .sample_size(10)
                    .warm_up_time(Duration::from_millis(500))
                    .measurement_time(Duration::from_secs(2));
                group.throughput(Throughput::Elements(QUERIES as u64));
                group.bench_function(format!("ef{ef}"), |b| {
                    b.iter(|| {
                        for q in &queries {
                            black_box(hnsw_rs_topk(&hnsw, q, ef));
                        }
                    });
                });
                group.finish();
            }
        }

        eprintln!("== N={n} summary (scale, ef, recall@10, avg_query_us) ==");
        for (scale, ef, recall, latency_us) in &summary {
            eprintln!("    scale={scale:.1}  ef={ef:>3}  {recall:.3}  {latency_us:.2}");
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
    targets = hv_ann_hubnsw
);

criterion_main!(benches);
