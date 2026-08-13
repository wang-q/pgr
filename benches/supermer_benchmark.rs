//! Super-mer / minimizer two-stage counting vs the direct radix path.
//!
//! Two datasets: the MG1655 genome itself (low redundancy) and a simulated
//! ~20x coverage read set sampled from a 1 Mb region (high redundancy, the
//! FastK design target). The stage-1 collapse should shrink the stage-2
//! sorting volume roughly by the coverage factor on the read set, while the
//! genome (mostly unique spans) shows the worst-case overhead.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pgr::libs::kmer::{count, supermer};
use pgr::libs::pgi::build::read_fasta;

const GENOME: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz");

fn genome_seqs() -> Vec<Vec<u8>> {
    read_fasta(GENOME)
        .unwrap()
        .into_iter()
        .map(|(_, s)| s)
        .collect()
}

/// ~20x coverage: 150 bp windows at a 7 bp stride over the first 1 Mb.
fn read_seqs() -> Vec<Vec<u8>> {
    let genome = genome_seqs();
    let contig = genome
        .iter()
        .map(Vec::as_slice)
        .find(|s| s.len() >= 1_000_000)
        .expect("mg1655 has a >= 1 Mb contig");
    let end = 1_000_000usize;
    (0..end - 150)
        .step_by(7)
        .map(|p| contig[p..p + 150].to_vec())
        .collect()
}

fn bench_supermer(c: &mut Criterion) {
    let genome = genome_seqs();
    let reads = read_seqs();
    let mut group = c.benchmark_group("kmer_count");
    group.sample_size(10);

    for (name, seqs) in [("genome", &genome), ("reads20x", &reads)] {
        for k in [17usize, 31, 100] {
            group.bench_with_input(
                BenchmarkId::new("direct", format!("{name}_k{k}")),
                &k,
                |b, &k| b.iter(|| black_box(count::build_table(seqs, k).unwrap().counts.len())),
            );
            group.bench_with_input(
                BenchmarkId::new("supermer", format!("{name}_k{k}")),
                &k,
                |b, &k| b.iter(|| black_box(supermer::build_table(seqs, k).unwrap().counts.len())),
            );
        }
    }
    group.finish();
}

criterion_group!(benches, bench_supermer);
criterion_main!(benches);
