//! Benchmark for the coarse pangenome graph build (seqwish-style DSU closure)
//! on the real 10-genome E. coli dataset (`benches/data/ecoli10.paf.gz`).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use flate2::read::GzDecoder;
use pgr::libs::paf::graph::PafGraph;
use std::io::Cursor;
use std::io::Read;

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/benches/data/ecoli10.paf.gz");

fn bench_graph_build(c: &mut Criterion) {
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
    let mut group = c.benchmark_group("paf_graph_build");
    group.bench_function("ecoli10_100", |b| {
        b.iter(|| {
            let g = PafGraph::build(Cursor::new(black_box(&bytes)), None, 100).unwrap();
            black_box(g.node_seqs.len() + g.edges.len())
        })
    });
    group.finish();
}

criterion_group!(benches, bench_graph_build);
criterion_main!(benches);
