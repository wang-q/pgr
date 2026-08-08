//! Real-cohort HV nearest-neighbor retrieval with ANI ground truth.
//!
//! Loads `.hv` files (produced by `pgr pgi to-hv`), compares exact cosine
//! scan, global HNSW (`hnsw_rs`), and species-routed per-clade HNSW.
//! Ground truths: (a) skani ANI top-10 (biological), (b) exact HV cosine
//! top-10 (graph-only error). See notes/design/genome-nn-query.md §7.4 #6/#7.
//!
//! Env: PGR_HV_REAL_DIR (default /tmp/hv_calib/hv135), PGR_HV_REAL_ANI
//! (default /tmp/hv_calib/ani.full.tsv), PGR_HV_REAL_META
//! (default /tmp/hv_calib/cohort.meta.tsv), PGR_HV_REAL_EFS.

use std::collections::HashMap;
use std::fs;
use std::io::Read;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use hnsw_rs::prelude::*;

const K: usize = 10;
const M_CONN: usize = 24;
const EF_GRID: [usize; 4] = [10, 20, 50, 100];

type HvGraph = Hnsw<'static, f32, DistL2>;

fn env_or(key: &str, default: &str) -> String {
    std::env::var(key).unwrap_or_else(|_| default.to_string())
}

fn parse_efs() -> Vec<usize> {
    let parsed: Vec<usize> = std::env::var("PGR_HV_REAL_EFS")
        .map(|s| s.split(',').filter_map(|x| x.trim().parse().ok()).collect())
        .unwrap_or_default();
    if parsed.is_empty() {
        EF_GRID.to_vec()
    } else {
        parsed
    }
}

fn short_name(path: &str) -> String {
    PathBuf::from(path)
        .parent()
        .and_then(|p| p.file_name())
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default()
}

fn read_hv(path: &PathBuf) -> (String, Vec<i32>) {
    let mut f = fs::File::open(path).expect("open hv");
    let mut magic = [0u8; 4];
    f.read_exact(&mut magic).expect("magic");
    assert_eq!(&magic, b"PGV1", "bad hv magic");
    let mut u32buf = [0u8; 4];
    let mut rd_u32 = |f: &mut fs::File| {
        f.read_exact(&mut u32buf).unwrap();
        u32::from_le_bytes(u32buf)
    };
    let _version = rd_u32(&mut f);
    let _k = rd_u32(&mut f) as usize;
    let dim = rd_u32(&mut f) as usize;
    let _sparse = rd_u32(&mut f) as usize;
    let mut u64buf = [0u8; 8];
    f.read_exact(&mut u64buf).unwrap();
    let _n_kmer = u64::from_le_bytes(u64buf) as usize;
    let nb = rd_u32(&mut f) as usize;
    let mut name = vec![0u8; nb];
    f.read_exact(&mut name).unwrap();
    let name = String::from_utf8(name).unwrap();
    let mut hv = Vec::with_capacity(dim);
    for _ in 0..dim {
        f.read_exact(&mut u32buf).unwrap();
        hv.push(i32::from_le_bytes(u32buf));
    }
    (name, hv)
}

fn l2_normalize(v: &[i32]) -> Vec<f32> {
    let mut acc = 0.0f64;
    for &x in v {
        acc += (x as f64) * (x as f64);
    }
    let norm = acc.sqrt() as f32;
    v.iter().map(|&x| x as f32 / norm).collect()
}

fn l2_dist(a: &[f32], b: &[f32]) -> f32 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y) as f64 * (x - y) as f64)
        .sum::<f64>()
        .sqrt() as f32
}

fn exact_topk(db: &[Vec<f32>], query: &[f32], k: usize) -> Vec<usize> {
    let mut scored: Vec<(f32, usize)> = db
        .iter()
        .enumerate()
        .map(|(i, v)| (l2_dist(v, query), i))
        .collect();
    scored.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    scored[..k].iter().map(|&(_, i)| i).collect()
}

fn build_graph(data: &[Vec<f32>]) -> HvGraph {
    let hnsw = Hnsw::new(M_CONN, data.len(), 16, 64, DistL2);
    for (i, v) in data.iter().enumerate() {
        hnsw.insert((v, i));
    }
    hnsw
}

fn hnsw_topk(hnsw: &HvGraph, query: &[f32], ef: usize) -> Vec<usize> {
    let mut nbs = hnsw.search(query, K, ef);
    nbs.sort_unstable_by(|a, b| a.distance.total_cmp(&b.distance));
    nbs.iter().map(|n| n.d_id).collect()
}

fn routed_topk(
    clade_graphs: &[HvGraph],
    clade_ids: &[Vec<usize>],
    reps: &[Vec<f32>],
    query: &[f32],
    r: usize,
    ef: usize,
) -> Vec<usize> {
    let mut scored: Vec<(f32, usize)> = reps
        .iter()
        .enumerate()
        .map(|(c, v)| (l2_dist(v, query), c))
        .collect();
    scored.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let mut pool: Vec<(f32, usize)> = Vec::new();
    for &(_, c) in scored.iter().take(r) {
        let mut nbs = clade_graphs[c].search(query, K, ef);
        nbs.sort_unstable_by(|a, b| a.distance.total_cmp(&b.distance));
        for n in nbs {
            pool.push((n.distance, clade_ids[c][n.d_id]));
        }
    }
    pool.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    pool.iter().take(K).map(|&(_, gid)| gid).collect()
}

fn hv_ann_real(c: &mut Criterion) {
    let dir = env_or("PGR_HV_REAL_DIR", "/tmp/hv_calib/hv135");
    let ani_path = env_or("PGR_HV_REAL_ANI", "/tmp/hv_calib/ani.full.tsv");
    let meta_path = env_or("PGR_HV_REAL_META", "/tmp/hv_calib/cohort.meta.tsv");
    let efs = parse_efs();

    let mut paths: Vec<PathBuf> = fs::read_dir(&dir)
        .expect("hv dir")
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().map(|e| e == "hv").unwrap_or(false))
        .collect();
    paths.sort();
    let mut names = Vec::with_capacity(paths.len());
    let mut db: Vec<Vec<f32>> = Vec::with_capacity(paths.len());
    for p in &paths {
        let (name, hv) = read_hv(p);
        names.push(name);
        db.push(l2_normalize(&hv));
    }
    let n = names.len();
    let idx: HashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    eprintln!("loaded {n} real HV files from {dir}");

    let mut ani = vec![vec![0.0f32; n]; n];
    let ani_raw = fs::read_to_string(&ani_path).expect("ani file");
    for line in ani_raw.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let a = short_name(f[0]);
        let b = short_name(f[1]);
        if let (Some(&ia), Some(&ib)) = (idx.get(a.as_str()), idx.get(b.as_str())) {
            if let Ok(v) = f[2].parse::<f32>() {
                ani[ia][ib] = v;
                ani[ib][ia] = v;
            }
        }
    }

    let meta_raw = fs::read_to_string(&meta_path).expect("meta file");
    let mut sp_of = vec![0usize; n];
    let mut species: Vec<String> = Vec::new();
    for line in meta_raw.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 2 {
            continue;
        }
        if let Some(&i) = idx.get(f[0]) {
            let s = f[1].to_string();
            let ci = species.iter().position(|x| *x == s).unwrap_or_else(|| {
                species.push(s.clone());
                species.len() - 1
            });
            sp_of[i] = ci;
        }
    }
    eprintln!("species (clades): {}", species.len());

    let ani_truth: Vec<Option<Vec<usize>>> = (0..n)
        .map(|i| {
            let mut known: Vec<usize> = (0..n).filter(|&j| j != i && ani[i][j] > 0.0).collect();
            if known.len() < K {
                return None;
            }
            known.sort_unstable_by(|&a, &b| ani[i][b].total_cmp(&ani[i][a]));
            Some(known[..K].to_vec())
        })
        .collect();
    let exact_truth: Vec<Vec<usize>> = (0..n).map(|i| exact_topk(&db, &db[i], K)).collect();

    let global = build_graph(&db);
    let mut clade_graphs: Vec<HvGraph> = Vec::new();
    let mut clade_ids: Vec<Vec<usize>> = Vec::new();
    let mut reps: Vec<Vec<f32>> = Vec::new();
    {
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); species.len()];
        for i in 0..n {
            members[sp_of[i]].push(i);
        }
        for c in 0..species.len() {
            let data: Vec<Vec<f32>> = members[c].iter().map(|&i| db[i].clone()).collect();
            clade_graphs.push(build_graph(&data));
            clade_ids.push(members[c].clone());
            reps.push(db[members[c][0]].clone());
        }
    }

    let mut summary = Vec::<(String, f64, f64, f64)>::new();
    let eval = |ann: &Vec<usize>, truth: &Option<Vec<usize>>| -> (f64, bool) {
        match truth {
            Some(tr) => (
                ann.iter().filter(|x| tr.contains(x)).count() as f64 / K as f64,
                true,
            ),
            None => (0.0, false),
        }
    };

    let measure = |label: String,
                   arg: &dyn Fn(&[f32], usize) -> Vec<usize>,
                   ef: usize,
                   summary: &mut Vec<(String, f64, f64, f64)>| {
        let mut hits_ani = 0.0f64;
        let mut hits_hv = 0.0f64;
        let mut counted = 0usize;
        let mut total = Duration::ZERO;
        for i in 0..n {
            let t0 = Instant::now();
            let ann = arg(&db[i], ef);
            total += t0.elapsed();
            let (ha, ok) = eval(&ann, &ani_truth[i]);
            hits_ani += ha;
            counted += ok as usize;
            hits_hv += ann.iter().filter(|x| exact_truth[i].contains(x)).count() as f64 / K as f64;
        }
        let ra = hits_ani / counted.max(1) as f64;
        let rh = hits_hv / n as f64;
        let us = total.as_secs_f64() * 1e6 / n as f64;
        summary.push((label.clone(), ra, rh, us));
        eprintln!("{label}: recall_ANI@10={ra:.3} recall_HV@10={rh:.3} avg_query={us:.1} us");
    };

    measure(
        "exact".into(),
        &|q, _| exact_topk(&db, q, K),
        0,
        &mut summary,
    );
    for &ef in &efs {
        measure(
            format!("hnsw ef{ef}"),
            &|q, _| hnsw_topk(&global, q, ef),
            ef,
            &mut summary,
        );
        for r in [1usize, 2] {
            measure(
                format!("routed{r} ef{ef}"),
                &|q, _| routed_topk(&clade_graphs, &clade_ids, &reps, q, r, ef),
                ef,
                &mut summary,
            );
        }
    }

    eprintln!("== summary (variant, recall_ANI@10, recall_HV@10, avg_query_us) ==");
    for (v, ra, rh, us) in &summary {
        eprintln!("    {v:<14} {ra:.3}  {rh:.3}  {us:.1}");
    }

    let mut g = c.benchmark_group("real/hnsw");
    g.sample_size(10)
        .warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_secs(1));
    g.throughput(Throughput::Elements(n as u64));
    g.bench_function("ef50", |b| {
        b.iter(|| {
            for v in &db {
                black_box(hnsw_topk(&global, v, 50));
            }
        });
    });
    g.finish();

    let mut g = c.benchmark_group("real/exact");
    g.sample_size(10)
        .warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_secs(1));
    g.throughput(Throughput::Elements(n as u64));
    g.bench_function("scan", |b| {
        b.iter(|| {
            for v in &db {
                black_box(exact_topk(&db, v, K));
            }
        });
    });
    g.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = hv_ann_real
);

criterion_main!(benches);
