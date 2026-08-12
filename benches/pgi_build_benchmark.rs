//! Benchmarks for `pgi build` (syncmer sampling + packed radix sort +
//! grouping) on local E. coli genomes.
//!
//! Baseline (2026-08-12, after the packed-key refactor): see
//! `notes/benchmarks/bench-pgi-vs-gixmake.md`. These criterion numbers are
//! the regression gate for any future pgi-build change (the project rule is
//! that performance-sensitive changes must establish a benchmark baseline
//! first).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::pgi::build::{build_from_seqs, read_fasta};

const MG1655: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz");
const SAKAI: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/sakai.fa.gz");

fn bench_build(c: &mut Criterion) {
    let mg = read_fasta(MG1655).unwrap();
    let sk = read_fasta(SAKAI).unwrap();
    let mut two = mg.clone();
    two.extend(sk.clone());

    let mut group = c.benchmark_group("pgi_build");
    group.sample_size(10);
    group.bench_function("mg1655_k40", |b| {
        b.iter(|| {
            let idx = build_from_seqs(mg.clone(), 40, 8, 5, false, false).unwrap();
            black_box(idx.entries.len())
        })
    });
    group.bench_function("mg1655_sakai_k40", |b| {
        b.iter(|| {
            let idx = build_from_seqs(two.clone(), 40, 8, 5, false, false).unwrap();
            black_box(idx.entries.len())
        })
    });
    group.bench_function("mg1655_k21", |b| {
        b.iter(|| {
            let idx = build_from_seqs(mg.clone(), 21, 8, 5, false, false).unwrap();
            black_box(idx.entries.len())
        })
    });
}

criterion_group!(benches, bench_build);
criterion_main!(benches);
