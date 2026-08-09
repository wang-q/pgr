use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::hash::seq_mins;
use pgr::libs::pgi::build::read_fasta;

const GENOME: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz");

fn bench_minimizer(c: &mut Criterion) {
    let seqs: Vec<Vec<u8>> = read_fasta(GENOME)
        .unwrap()
        .into_iter()
        .map(|(_, s)| s)
        .collect();
    let mut group = c.benchmark_group("minimizer");
    group.sample_size(10);
    group.bench_function("mg1655_k21_w5_rapid", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for seq in &seqs {
                let set = seq_mins(seq, "rapid", 21, 5).unwrap();
                total += set.len();
            }
            black_box(total)
        })
    });
    group.bench_function("mg1655_k21_w5_fx", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for seq in &seqs {
                let set = seq_mins(seq, "fx", 21, 5).unwrap();
                total += set.len();
            }
            black_box(total)
        })
    });
    group.finish();
}

criterion_group!(benches, bench_minimizer);
criterion_main!(benches);
