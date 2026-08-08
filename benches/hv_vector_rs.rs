//! Real-cohort sqlite-vector-rs (usearch HNSW via PGVector-like vtab) vs exact scan.
//!
//! Loads the real Enterobacterales HV cohort (i32 4096, produced by
//! `pgr pgi to-hv`) and measures:
//!   (a) exact f32 L2 top-10 scan — latency and recall baseline (#10);
//!   (b) sqlite-vector-rs `VectorTable` vtab (usearch HNSW) — build time,
//!       per-query latency, top-10 recall vs exact, DB file size, and
//!       reload-from-shadow-table cost.
//!
//! The crate's library-mode `register()` is `todo!()` (0.1.0), so the vtab is
//! registered manually with sqlite3-ext-vtab exactly as its loadable-extension
//! entry point does. HNSW parameters are tunable per table via
//! `CREATE VIRTUAL TABLE ... USING vector(..., m=, ef_construction=, ef_search=)`.
//!
//! Env: PGR_HV_REAL_DIR (default /tmp/hv_calib/hv2115),
//!      PGR_HV_REAL_EFS (comma list, default 10,20,50,100,200),
//!      PGR_HV_REAL_NQ  (queries measured, default min(n, 200)).

use std::error::Error;
use std::fs;
use std::io::Read;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use sqlite3_ext::vtab::Module;
use sqlite3_ext::Blob;
use sqlite3_ext::FallibleIteratorMut;
use sqlite3_ext::FromValue;
use sqlite3_ext::Value;
use sqlite_vector_rs::distance::DistanceMetric;
use sqlite_vector_rs::index::{HnswIndex, HnswParams};
use sqlite_vector_rs::scalar;
use sqlite_vector_rs::types::VectorType;
use sqlite_vector_rs::vtab::VectorTable;

const K: usize = 10;
const DIM: usize = 4096;
const EF_GRID: [usize; 5] = [10, 20, 50, 100, 200];

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

fn read_hv(path: &PathBuf) -> Vec<i32> {
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
    let mut hv = Vec::with_capacity(dim);
    for _ in 0..dim {
        f.read_exact(&mut u32buf).unwrap();
        hv.push(i32::from_le_bytes(u32buf));
    }
    hv
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
        .sum::<f64>() as f32
}

fn exact_topk(db: &[Vec<f32>], query: &[f32], k: usize) -> Vec<u64> {
    let mut scored: Vec<(f32, u64)> = db
        .iter()
        .enumerate()
        .map(|(i, v)| (l2_dist(v, query), i as u64))
        .collect();
    scored.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    scored[..k].iter().map(|&(_, id)| id).collect()
}

fn f32_blob(v: &[f32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

/// Register the `vector` vtab module + scalar functions, exactly like the
/// crate's loadable-extension entry point (library `register()` is a stub).
fn register_extension(db: &sqlite3_ext::Connection) -> Result<(), Box<dyn Error>> {
    let module = sqlite3_ext::vtab::StandardModule::<VectorTable<'_>>::new()
        .with_update()
        .with_transactions()
        .with_find_function();
    db.create_module("vector", module, ())?;
    scalar::register_scalar_functions(db)?;
    Ok(())
}

fn create_vtab(
    db: &sqlite3_ext::Connection,
    table: &str,
    ef_search: usize,
) -> Result<(), Box<dyn Error>> {
    let sql = format!(
        "CREATE VIRTUAL TABLE {table} USING vector(\
         dim={DIM}, type=float4, metric=l2, m=16, ef_construction=200, ef_search={ef_search})"
    );
    db.execute(&sql, ())?;
    Ok(())
}

/// Insert all cohort vectors, wrapped in an explicit transaction so the HNSW
/// graph is serialized to the shadow table once (at COMMIT), not per row.
fn insert_all(
    db: &sqlite3_ext::Connection,
    table: &str,
    data: &[Vec<f32>],
) -> Result<Duration, Box<dyn Error>> {
    let t0 = Instant::now();
    db.execute("BEGIN", ())?;
    {
        // The vtab's shadow table is `id INTEGER PRIMARY KEY AUTOINCREMENT`;
        // user-supplied ids are ignored, so insert only the vector column and
        // map returned ids back via `id - 1` (ids start at 1).
        let sql = format!("INSERT INTO {table}(vector) VALUES(?)");
        let mut stmt = db.prepare(&sql)?;
        for v in data.iter() {
            let blob = f32_blob(v);
            stmt.query([Value::Blob(blob.as_slice().into())])?;
            while stmt.next()?.is_some() {}
        }
    }
    db.execute("COMMIT", ())?;
    Ok(t0.elapsed())
}

fn knn_query(
    db: &sqlite3_ext::Connection,
    table: &str,
    query_blob: &[u8],
    k: usize,
) -> Result<Vec<(i64, f64)>, Box<dyn Error>> {
    let sql = format!("SELECT id, distance FROM {table} WHERE knn_match(distance, ?) LIMIT {k}");
    let mut stmt = db.prepare(&sql)?;
    stmt.query([Value::Blob(Blob::from(query_blob))])?;
    let mut out = Vec::with_capacity(k);
    while let Some(row) = stmt.next()? {
        out.push((row[0].get_i64(), row[1].get_f64()));
    }
    Ok(out)
}

fn recall_hv(ann: &[(i64, f64)], truth: &[u64]) -> f64 {
    ann.iter()
        .filter(|(id, _)| truth.contains(&(*id as u64)))
        .count() as f64
        / K as f64
}

/// vtab ids are 1-based (shadow AUTOINCREMENT); map back to 0-based cohort idx.
fn recall_hv_vtab(ann: &[(i64, f64)], truth: &[u64]) -> f64 {
    ann.iter()
        .filter(|(id, _)| truth.contains(&(*id as u64 - 1)))
        .count() as f64
        / K as f64
}

fn hv_vector_rs(c: &mut Criterion) {
    let dir = env_or("PGR_HV_REAL_DIR", "/tmp/hv_calib/hv2115");
    let efs = parse_efs();

    let mut paths: Vec<PathBuf> = fs::read_dir(&dir)
        .expect("hv dir")
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().map(|e| e == "hv").unwrap_or(false))
        .collect();
    paths.sort();
    let db: Vec<Vec<f32>> = paths.iter().map(|p| l2_normalize(&read_hv(p))).collect();
    let n = db.len();
    assert_eq!(db[0].len(), DIM, "HV dimension != 4096");
    eprintln!("loaded {n} real HV files from {dir}");

    let nq = std::env::var("PGR_HV_REAL_NQ")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(n.min(200));
    let truth: Vec<Vec<u64>> = db.iter().take(nq).map(|q| exact_topk(&db, q, K)).collect();
    eprintln!("exact top-{K} computed for {nq} queries");

    // --- exact scan latency (baseline; #10's numpy path was ~2.46 ms/query) ---
    let mut total = Duration::ZERO;
    for i in 0..nq {
        let t0 = Instant::now();
        black_box(exact_topk(&db, &db[i], K));
        total += t0.elapsed();
    }
    let exact_us = total.as_secs_f64() * 1e6 / nq as f64;
    eprintln!("exact scan: avg_query = {exact_us:.1} us");

    // --- HnswIndex ef_search grid (usearch core, no SQLite layer) ---
    for &ef in &efs {
        let params = HnswParams {
            m: 16,
            ef_construction: 200,
            ef_search: ef,
        };
        let idx = HnswIndex::new(DIM, VectorType::Float4, DistanceMetric::L2, Some(params))
            .expect("usearch index");
        let t0 = Instant::now();
        for (i, v) in db.iter().enumerate() {
            let blob = f32_blob(v);
            idx.add(i as u64, &blob).expect("usearch add");
        }
        let build_s = t0.elapsed().as_secs_f64();
        let mut qus = 0.0f64;
        let mut hits = 0.0f64;
        for i in 0..nq {
            let t0 = Instant::now();
            let blob = f32_blob(&db[i]);
            let res = idx.search(&blob, K).expect("usearch search");
            qus += t0.elapsed().as_secs_f64() * 1e6;
            hits += recall_hv(
                &res.iter()
                    .map(|&(id, d)| (id as i64, d as f64))
                    .collect::<Vec<_>>(),
                &truth[i],
            );
        }
        eprintln!(
            "usearch ef{ef:<3}: build = {build_s:.2} s, avg_query = {:.1} us, recall_HV@10 = {:.3}",
            qus / nq as f64,
            hits / nq as f64
        );
    }

    // --- sqlite-vector-rs vtab: file DB, transaction-wrapped bulk insert ---
    let tmp = tempfile::NamedTempFile::new().expect("temp db");
    let db_path = tmp.path().to_owned();
    {
        eprintln!("[vtab] opening db");
        let conn = sqlite3_ext::Database::open(&db_path).expect("open sqlite db");
        eprintln!("[vtab] registering extension");
        register_extension(&conn).expect("register vector extension");
        eprintln!("[vtab] creating vtab");
        create_vtab(&conn, "emb", 64).expect("create vtab");
        eprintln!("[vtab] inserting");
        let build = insert_all(&conn, "emb", &db).expect("bulk insert");
        eprintln!(
            "vtab(ef64) build: {:.2} s (BEGIN/COMMIT wrapped; HNSW serialized once at COMMIT)",
            build.as_secs_f64()
        );

        let mut qus = 0.0f64;
        let mut hits = 0.0f64;
        for i in 0..nq {
            let t0 = Instant::now();
            let blob = f32_blob(&db[i]);
            let res = knn_query(&conn, "emb", &blob, K).expect("knn query");
            qus += t0.elapsed().as_secs_f64() * 1e6;
            hits += recall_hv_vtab(&res, &truth[i]);
        }
        let vtab_us = qus / nq as f64;
        eprintln!(
            "vtab(ef64) query: avg = {vtab_us:.1} us, recall_HV@10 = {:.3}",
            hits / nq as f64
        );
    }
    let db_size = fs::metadata(&db_path).map(|m| m.len()).unwrap_or(0);
    eprintln!("vtab(ef64) db file size: {:.1} MB", db_size as f64 / 1e6);

    // --- reopen: reload HNSW from the shadow table, then query again ---
    {
        eprintln!("[vtab] reopening db");
        let t0 = Instant::now();
        let conn = sqlite3_ext::Database::open(&db_path).expect("reopen sqlite db");
        eprintln!("[vtab] re-registering extension");
        register_extension(&conn).expect("register vector extension");
        let load_s = t0.elapsed().as_secs_f64();
        let mut qus = 0.0f64;
        let mut hits = 0.0f64;
        for i in 0..nq {
            let t0 = Instant::now();
            let blob = f32_blob(&db[i]);
            let res = knn_query(&conn, "emb", &blob, K).expect("knn query");
            qus += t0.elapsed().as_secs_f64() * 1e6;
            hits += recall_hv_vtab(&res, &truth[i]);
        }
        eprintln!(
            "vtab(ef64) reload: {load_s:.2} s (connect+register), avg_query = {:.1} us, recall = {:.3}",
            qus / nq as f64,
            hits / nq as f64
        );
        // Second pass on the same reopened connection: excludes one-time reload cost.
        let mut qus2 = 0.0f64;
        for v in db.iter().take(nq.min(50)) {
            let t0 = Instant::now();
            let blob = f32_blob(v);
            black_box(knn_query(&conn, "emb", &blob, K).expect("knn query"));
            qus2 += t0.elapsed().as_secs_f64() * 1e6;
        }
        eprintln!(
            "vtab(ef64) warm: avg_query = {:.1} us",
            qus2 / nq.min(50) as f64
        );
    }

    let mut g = c.benchmark_group("vector_rs/query");
    g.sample_size(10)
        .warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_secs(1));
    g.throughput(Throughput::Elements(nq as u64));
    g.bench_function("exact", |b| {
        b.iter(|| {
            for i in 0..nq {
                black_box(exact_topk(&db, &db[i], K));
            }
        });
    });
    g.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10).warm_up_time(Duration::from_millis(300)).measurement_time(Duration::from_secs(1));
    targets = hv_vector_rs
);
criterion_main!(benches);
