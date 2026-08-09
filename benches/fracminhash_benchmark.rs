use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::hash::seq_fracminhash;
use pgr::libs::pgi::build::read_fasta;

const GENOME: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz");

/// Canonical 2-bit k-mer keys only (rolling pack + rc), no rapidhash: the
/// fraction of `seq_fracminhash` spent before the hasher. Returns an
/// accumulator over canonical keys so the loop cannot be folded.
fn canonical_keys_only(seq: &[u8], k: usize) -> u128 {
    let mut acc = 0u128;
    for key in pgr::libs::nt::rolling_kmer_keys(seq, k) {
        let Some(key) = key else { continue };
        let canonical = key.min(pgr::libs::nt::rc_key(key, k));
        acc = acc.wrapping_add(canonical);
    }
    acc
}

fn bench_fracminhash(c: &mut Criterion) {
    let seqs: Vec<Vec<u8>> = read_fasta(GENOME)
        .unwrap()
        .into_iter()
        .map(|(_, s)| s)
        .collect();
    let mut group = c.benchmark_group("fracminhash");
    group.sample_size(10);
    group.bench_function("mg1655_k21_scale1000", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for seq in &seqs {
                let set = seq_fracminhash(seq, 21, 1000, false).unwrap();
                total += set.len();
            }
            black_box(total)
        })
    });
    group.bench_function("mg1655_k21_scale100", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for seq in &seqs {
                let set = seq_fracminhash(seq, 21, 100, false).unwrap();
                total += set.len();
            }
            black_box(total)
        })
    });
    group.bench_function("canonical_keys_only_k21", |b| {
        b.iter(|| {
            let mut total = 0u128;
            for seq in &seqs {
                total = total.wrapping_add(canonical_keys_only(seq, 21));
            }
            black_box(total)
        })
    });
    group.finish();
}

criterion_group!(benches, bench_fracminhash);
criterion_main!(benches);
