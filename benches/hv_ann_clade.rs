//! Knowledge-routed per-clade HNSW vs global HNSW for 4096-d HV.
//!
//! Synthetic cohort with explicit clade structure: each genome shares a
//! small global k-mer core + a clade-specific core + private k-mers, so
//! clade mates are far more similar than cross-clade genomes. Compares a
//! single global HNSW against per-clade HNSW graphs with query routing by
//! exact distance to clade representatives (a stand-in for external
//! knowledge such as phylogeny/traits; see notes/design/genome-nn-query.md
//! §6.5). Ground truth is the global exact cosine top-10.

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
const R_GRID: [usize; 3] = [1, 2, 4];
const M_CONN: usize = 24;

const GLOBAL_CORE: usize = 256;
const CLADE_CORE: usize = 2048;
const PRIVATE_MIN: usize = 512;
const PRIVATE_MAX: usize = 4096;

type HvGraph = Hnsw<'static, f32, DistL2>;

fn l2_normalize(v: &[i32]) -> &'static [f32] {
    let mut acc = 0.0f64;
    for &x in v {
        acc += (x as f64) * (x as f64);
    }
    let norm = acc.sqrt() as f32;
    let out: Vec<f32> = v.iter().map(|&x| x as f32 / norm).collect();
    Box::leak(out.into_boxed_slice())
}

fn genome_hv(global_core: &[u64], clade_core: &[u64], rng: &mut SmallRng) -> &'static [f32] {
    let private_len = PRIVATE_MIN + rng.random_range(0..=PRIVATE_MAX - PRIVATE_MIN);
    let mut hashes = Vec::with_capacity(GLOBAL_CORE + CLADE_CORE + private_len);
    hashes.extend_from_slice(global_core);
    hashes.extend_from_slice(clade_core);
    hashes.extend((0..private_len).map(|_| rng.random::<u64>()));
    l2_normalize(&hash_hv_bit(&hashes, D))
}

struct Cohort {
    db: Vec<&'static [f32]>,
    queries: Vec<&'static [f32]>,
    query_clades: Vec<usize>,
}

fn generate_cohort(seed: u64, n: usize, clades: usize) -> Cohort {
    assert!(
        n.is_multiple_of(clades),
        "N must be divisible by clade count"
    );
    let n_per = n / clades;
    let mut core_rng = SmallRng::seed_from_u64(seed ^ 0xC0FFEE);
    let global_core: Vec<u64> = (0..GLOBAL_CORE).map(|_| core_rng.random()).collect();
    let mut clade_cores = Vec::with_capacity(clades);
    for c in 0..clades {
        let mut r = SmallRng::seed_from_u64(seed ^ (0x1000 + c as u64));
        clade_cores.push((0..CLADE_CORE).map(|_| r.random()).collect::<Vec<u64>>());
    }

    let mut db = Vec::with_capacity(n);
    for i in 0..n {
        let c = i / n_per;
        let mut rng = SmallRng::seed_from_u64(seed ^ (i as u64 + 1));
        db.push(genome_hv(&global_core, &clade_cores[c], &mut rng));
    }

    let mut queries = Vec::with_capacity(QUERIES);
    let mut query_clades = Vec::with_capacity(QUERIES);
    for j in 0..QUERIES {
        let c = j % clades;
        let mut rng = SmallRng::seed_from_u64(seed ^ (1_000_003 + j as u64));
        queries.push(genome_hv(&global_core, &clade_cores[c], &mut rng));
        query_clades.push(c);
    }
    Cohort {
        db,
        queries,
        query_clades,
    }
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

fn build_graph(data: &[&'static [f32]]) -> HvGraph {
    let ef_c = ef_construction();
    let hnsw = Hnsw::new(M_CONN, data.len(), 16, ef_c, DistL2);
    let t0 = Instant::now();
    for (i, &v) in data.iter().enumerate() {
        hnsw.insert((v, i));
    }
    eprintln!(
        "    build N={} ef_c={}: {:?}",
        data.len(),
        ef_c,
        t0.elapsed()
    );
    hnsw
}

fn search_topk(hnsw: &HvGraph, query: &[f32], ef: usize) -> Vec<(f32, usize)> {
    let mut nbs = hnsw.search(query, K, ef);
    nbs.sort_unstable_by(|a, b| a.distance.total_cmp(&b.distance));
    nbs.iter().map(|n| (n.distance, n.d_id)).collect()
}

fn routed_topk(
    clade_graphs: &[HvGraph],
    reps: &[&'static [f32]],
    query: &[f32],
    r: usize,
    n_per: usize,
    ef: usize,
) -> Vec<usize> {
    let mut scored: Vec<(f32, usize)> = reps
        .iter()
        .enumerate()
        .map(|(c, v)| (dot_product(query, v), c))
        .collect();
    scored.sort_unstable_by(|a, b| b.0.total_cmp(&a.0));

    let mut pool: Vec<(f32, usize)> = Vec::new();
    for &(_, c) in scored.iter().take(r) {
        for (dist, id) in search_topk(&clade_graphs[c], query, ef) {
            pool.push((dist, c * n_per + id));
        }
    }
    pool.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    pool.iter().take(K).map(|&(_, gid)| gid).collect()
}

fn parse_sizes() -> Vec<usize> {
    let parsed: Vec<usize> = std::env::var("PGR_HV_ANN_SIZES")
        .map(|s| s.split(',').filter_map(|x| x.trim().parse().ok()).collect())
        .unwrap_or_default();
    if parsed.is_empty() {
        vec![10_000, 30_000]
    } else {
        parsed
    }
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

fn ef_construction() -> usize {
    std::env::var("PGR_HV_ANN_EFC")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(64)
}

fn clades() -> usize {
    std::env::var("PGR_HV_ANN_CLADES")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(16)
}

fn hv_ann_clade(c: &mut Criterion) {
    for &n in &parse_sizes() {
        let clades = clades();
        let cohort = generate_cohort(0x5EED_2026, n, clades);
        let n_per = n / clades;
        let efs = parse_efs();
        let exact: Vec<Vec<usize>> = cohort
            .queries
            .iter()
            .map(|q| exact_topk(&cohort.db, q, K))
            .collect();

        let global = build_graph(&cohort.db);
        let mut clade_graphs = Vec::with_capacity(clades);
        for c in 0..clades {
            clade_graphs.push(build_graph(&cohort.db[c * n_per..(c + 1) * n_per]));
        }
        let reps: Vec<&'static [f32]> = (0..clades).map(|c| cohort.db[c * n_per]).collect();

        // Routing accuracy: how often the query's true clade ranks in top-R
        // by exact distance to clade representatives.
        let mut own_rank = Vec::new();
        for (qi, q) in cohort.queries.iter().enumerate() {
            let mut scored: Vec<(f32, usize)> = reps
                .iter()
                .enumerate()
                .map(|(c, v)| (dot_product(q, v), c))
                .collect();
            scored.sort_unstable_by(|a, b| b.0.total_cmp(&a.0));
            own_rank.push(
                scored
                    .iter()
                    .position(|&(_, c)| c == cohort.query_clades[qi])
                    .unwrap_or(usize::MAX)
                    + 1,
            );
        }
        for r in &R_GRID {
            let acc = own_rank.iter().filter(|&&rank| rank <= *r).count();
            eprintln!(
                "== N={n} routing accuracy: own clade in top-{r} reps = {}/{}",
                acc, QUERIES
            );
        }

        let mut summary = Vec::<(String, usize, f64, f64)>::new();

        // Global HNSW.
        eprintln!("== N={n} global HNSW ==");
        for &ef in &efs {
            let mut hits = 0usize;
            let mut total = Duration::ZERO;
            for (qi, q) in cohort.queries.iter().enumerate() {
                let t0 = Instant::now();
                let ann: Vec<usize> = search_topk(&global, q, ef)
                    .iter()
                    .map(|&(_, id)| id)
                    .collect();
                total += t0.elapsed();
                hits += ann.iter().filter(|i| exact[qi].contains(i)).count();
            }
            let recall = hits as f64 / (QUERIES as f64 * K as f64);
            let latency_us = total.as_secs_f64() * 1e6 / QUERIES as f64;
            summary.push(("global".into(), ef, recall, latency_us));
            eprintln!("    ef={ef:>3}  recall@10={recall:.3}  avg_query={latency_us:.2} us");
        }

        // Routed per-clade HNSW.
        for &r in &R_GRID {
            eprintln!("== N={n} routed top-{r} clades ==");
            for &ef in &efs {
                let mut hits = 0usize;
                let mut total = Duration::ZERO;
                for (qi, q) in cohort.queries.iter().enumerate() {
                    let t0 = Instant::now();
                    let ann = routed_topk(&clade_graphs, &reps, q, r, n_per, ef);
                    total += t0.elapsed();
                    hits += ann.iter().filter(|i| exact[qi].contains(i)).count();
                }
                let recall = hits as f64 / (QUERIES as f64 * K as f64);
                let latency_us = total.as_secs_f64() * 1e6 / QUERIES as f64;
                summary.push((format!("routed{r}"), ef, recall, latency_us));
                eprintln!("    ef={ef:>3}  recall@10={recall:.3}  avg_query={latency_us:.2} us");
            }
        }

        eprintln!("== N={n} summary (variant, ef, recall@10, avg_query_us) ==");
        for (variant, ef, recall, latency_us) in &summary {
            eprintln!("    {variant:<7}  ef={ef:>3}  {recall:.3}  {latency_us:.2}");
        }

        // Criterion timing: global + routed (R=1/2/4) at each ef.
        for &ef in &efs {
            let mut group = c.benchmark_group(format!("global/{n}"));
            group
                .sample_size(10)
                .warm_up_time(Duration::from_millis(500))
                .measurement_time(Duration::from_secs(2));
            group.throughput(Throughput::Elements(QUERIES as u64));
            group.bench_function(format!("ef{ef}"), |b| {
                b.iter(|| {
                    for q in &cohort.queries {
                        black_box(search_topk(&global, q, ef));
                    }
                });
            });
            group.finish();

            for &r in &R_GRID {
                let mut group = c.benchmark_group(format!("routed{r}/{n}"));
                group
                    .sample_size(10)
                    .warm_up_time(Duration::from_millis(500))
                    .measurement_time(Duration::from_secs(2));
                group.throughput(Throughput::Elements(QUERIES as u64));
                group.bench_function(format!("ef{ef}"), |b| {
                    b.iter(|| {
                        for q in &cohort.queries {
                            black_box(routed_topk(&clade_graphs, &reps, q, r, n_per, ef));
                        }
                    });
                });
                group.finish();
            }
        }

        let mut group = c.benchmark_group(format!("exact/{n}"));
        group
            .sample_size(10)
            .warm_up_time(Duration::from_millis(500))
            .measurement_time(Duration::from_secs(2));
        group.throughput(Throughput::Elements(QUERIES as u64));
        group.bench_function("scan", |b| {
            b.iter(|| {
                for q in &cohort.queries {
                    black_box(exact_topk(&cohort.db, q, K));
                }
            });
        });
        group.finish();
    }
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = hv_ann_clade
);

criterion_main!(benches);
