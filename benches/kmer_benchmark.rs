//! Benchmarks for the native k-mer repeat pipeline: canonical count-table
//! build and per-position profile generation on MG1655
//! (`tests/genome/mg1655.fa.gz`, 4.6 Mb).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::kmer;
use pgr::libs::pgi::build::read_fasta;

const GENOME: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz");
const LIB: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/pgr/tncentral.fa.gz");

fn bench_build_and_profiles(c: &mut Criterion) {
    let k = 17usize;
    let seqs: Vec<Vec<u8>> = read_fasta(GENOME)
        .unwrap()
        .into_iter()
        .map(|(_, s)| s)
        .collect();
    let lib_seqs: Vec<Vec<u8>> = read_fasta(LIB)
        .unwrap()
        .into_iter()
        .map(|(_, s)| s)
        .collect();

    let mut group = c.benchmark_group("kmer");
    group.sample_size(10);
    group.bench_function("count_mg1655", |b| {
        b.iter(|| {
            let table = kmer::count::build_table(&seqs, k).unwrap();
            black_box(table.keys.len())
        })
    });
    group.bench_function("self_profiles_mg1655", |b| {
        let table = kmer::count::build_table(&seqs, k).unwrap();
        b.iter(|| {
            let profiles = kmer::profile::self_profiles(&seqs, k, &table);
            black_box(profiles[0].len())
        })
    });
    group.bench_function("relative_profiles_mg1655", |b| {
        let table = kmer::count::build_table(&lib_seqs, k).unwrap();
        b.iter(|| {
            let profiles = kmer::profile::relative_profiles(&seqs, k, &table);
            black_box(profiles[0].len())
        })
    });
    group.finish();
}

criterion_group!(benches, bench_build_and_profiles);
criterion_main!(benches);
