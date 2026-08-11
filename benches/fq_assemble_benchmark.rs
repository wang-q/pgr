//! Benchmarks for `pgr fq assemble` (tadpole-compatible contig assembly)
//! on the BBTools Lambda 20k-read paired dataset: the quality-gated k-mer
//! count-table build and the full assembly (bubbles on and off).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use pgr::libs::fq::assemble::{assemble, AssembleOptions};
use pgr::libs::fq::tadpole::TadpoleTable;
use std::sync::OnceLock;

const R1: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/bbtools/Lambda/R1.fq.gz");
const R2: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/bbtools/Lambda/R2.fq.gz");

/// The paired 20k-read Lambda reads, loaded once into memory.
fn reads() -> &'static Vec<(Vec<u8>, Vec<u8>)> {
    static READS: OnceLock<Vec<(Vec<u8>, Vec<u8>)>> = OnceLock::new();
    READS.get_or_init(|| {
        let mut out = Vec::new();
        for path in [R1, R2] {
            let mut reader = SeqReader::new(path).unwrap();
            let mut rec = SeqRecord::new();
            while reader.read_record(&mut rec).unwrap() {
                out.push((rec.sequence().to_vec(), rec.quality_scores().to_vec()));
            }
        }
        out
    })
}

fn infiles() -> &'static [String] {
    static INFILES: OnceLock<Vec<String>> = OnceLock::new();
    INFILES.get_or_init(|| vec![R1.to_string(), R2.to_string()])
}

fn bench_build_and_assemble(c: &mut Criterion) {
    let reads = reads();
    let mut group = c.benchmark_group("fq_assemble");
    group.sample_size(10);
    group.bench_function("tadpole_table_build_k31_20k", |b| {
        b.iter(|| {
            let table = TadpoleTable::build(reads, 31, 0.5);
            black_box(table)
        })
    });
    group.bench_function("assemble_full_k31_20k_no_bubbles", |b| {
        let opts = AssembleOptions {
            k: 31,
            pop_bubbles: false,
            ..AssembleOptions::default()
        };
        b.iter(|| {
            let mut out = Vec::new();
            let stats = assemble(infiles(), &mut out, &opts).unwrap();
            black_box(stats.contigs_built)
        })
    });
    group.bench_function("assemble_full_k31_20k_bubbles", |b| {
        let opts = AssembleOptions {
            k: 31,
            ..AssembleOptions::default()
        };
        b.iter(|| {
            let mut out = Vec::new();
            let stats = assemble(infiles(), &mut out, &opts).unwrap();
            black_box(stats.contigs_built)
        })
    });
    group.finish();
}

criterion_group!(benches, bench_build_and_assemble);
criterion_main!(benches);
