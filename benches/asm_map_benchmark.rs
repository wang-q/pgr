//! Benchmarks for `pgr asm map` (perfect-match mapping): reference index
//! build and the full mapping of the BBTools Lambda 2k paired reads against
//! the tadpole contig golden (the same reads assembled).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::map::{build_index, map_files, read_fasta, MapOptions, RefRecord};
use std::sync::OnceLock;

const REF: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/bbtools/Lambda/golden/tadpole_contigs31.fasta.gz"
);
const R1: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/bbtools/Lambda/R1.2k.fq.gz"
);
const R2: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/bbtools/Lambda/R2.2k.fq.gz"
);

fn refs() -> &'static Vec<RefRecord> {
    static REFS: OnceLock<Vec<RefRecord>> = OnceLock::new();
    REFS.get_or_init(|| read_fasta(&[REF.to_string()]).unwrap())
}

fn read_paths() -> &'static [String] {
    static PATHS: OnceLock<Vec<String>> = OnceLock::new();
    PATHS.get_or_init(|| vec![R1.to_string(), R2.to_string()])
}

fn bench_map(c: &mut Criterion) {
    let mut group = c.benchmark_group("asm_map");
    group.sample_size(10);
    group.bench_function("build_index_lambda_contigs", |b| {
        b.iter(|| {
            let index = build_index(refs(), 31).unwrap();
            black_box(index)
        })
    });
    group.bench_function("map_lambda_2k_no_outputs", |b| {
        let opts = MapOptions {
            k: 31,
            outm: None,
            outu: None,
            paired: false,
            max_reads: None,
        };
        b.iter(|| {
            let stats = map_files(refs(), read_paths(), &opts).unwrap();
            black_box(stats.mapped)
        })
    });
    group.finish();
}

criterion_group!(benches, bench_map);
criterion_main!(benches);
