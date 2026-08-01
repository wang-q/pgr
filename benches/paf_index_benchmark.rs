//! Benchmarks for the PAF interval index: in-memory build, single-hop query,
//! and transitive BFS, on the real 10-genome E. coli dataset
//! (`benches/data/ecoli10.paf.gz`, produced by `scripts/verify-pangenome.sh`
//! with the fixture variant; see notes/paf-pangenome.md §5.0).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use flate2::read::GzDecoder;
use pgr::libs::paf::index::PafIndex;
use std::io::Cursor;
use std::io::Read;

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/benches/data/ecoli10.paf.gz");
const SEED: (&str, i32, i32) = ("mg1655.NC_000913", 100_000, 110_000);

fn fixture_bytes() -> Vec<u8> {
    let f = std::fs::File::open(FIXTURE).unwrap_or_else(|_| {
        panic!(
            "benchmark fixture missing: regenerate with the fixture variant of \
             scripts/verify-pangenome.sh into {FIXTURE}"
        )
    });
    let mut bytes = Vec::new();
    GzDecoder::new(f)
        .read_to_end(&mut bytes)
        .expect("gzip fixture");
    bytes
}

fn bench_index_build(c: &mut Criterion) {
    let bytes = fixture_bytes();
    let mut group = c.benchmark_group("paf_index_build");
    group.bench_function("ecoli10", |b| {
        b.iter(|| PafIndex::build(Cursor::new(black_box(&bytes))).unwrap())
    });
    group.finish();
}

fn bench_query(c: &mut Criterion) {
    let bytes = fixture_bytes();
    let idx = PafIndex::build(Cursor::new(&bytes)).unwrap();
    let tid = idx
        .name_to_id(SEED.0)
        .expect("seed contig present in index");
    let mut group = c.benchmark_group("paf_index_query");
    group.bench_function("ecoli10_mg1655_100k", |b| {
        b.iter(|| {
            let hits = idx.query(tid, SEED.1, SEED.2, 0.0, 0);
            black_box(hits.len())
        })
    });
    group.finish();
}

fn bench_bfs(c: &mut Criterion) {
    let bytes = fixture_bytes();
    let idx = PafIndex::build(Cursor::new(&bytes)).unwrap();
    let tid = idx
        .name_to_id(SEED.0)
        .expect("seed contig present in index");
    let mut group = c.benchmark_group("paf_index_bfs");
    group.bench_function("ecoli10_mg1655_100k_depth2", |b| {
        b.iter(|| {
            let hits = idx.query_transitive_bfs(tid, SEED.1, SEED.2, 2, 10, 10, 0.0, 0, 0, None);
            black_box(hits.len())
        })
    });
    group.finish();
}

criterion_group!(benches, bench_index_build, bench_query, bench_bfs);
criterion_main!(benches);
