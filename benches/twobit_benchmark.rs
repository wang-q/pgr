use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use pgr::libs::fmt::twobit::{classify_dna_with, Blocks, SimdPath};
use rand::Rng;

fn random_dna(rng: &mut impl Rng, len: usize) -> String {
    let alphabet = *b"ACGTacgtNMRWSYKVHDXn-*_ 0123456789Uu";
    (0..len)
        .map(|_| alphabet[rng.random_range(0..alphabet.len())] as char)
        .collect()
}

fn bench_twobit(c: &mut Criterion) {
    let mut rng = rand::rng();
    let dna = random_dna(&mut rng, 10_000_000);
    let mut group = c.benchmark_group("twobit_from_dna");
    for (name, path) in [
        ("scalar", SimdPath::Scalar),
        ("wide", SimdPath::Wide),
        ("avx2", SimdPath::Avx2),
    ] {
        group.bench_function(name, |b| {
            b.iter_batched(
                || black_box(dna.clone()),
                |d| {
                    let classes = classify_dna_with(path, d.as_bytes());
                    black_box((classes.n_mask.len(), classes.codes.len()))
                },
                BatchSize::LargeInput,
            )
        });
    }
    group.bench_function("from_dna_auto", |b| {
        b.iter_batched(
            || black_box(dna.clone()),
            |d| {
                let r = Blocks::from_dna(&d, true).unwrap();
                black_box((r.0.len(), r.3))
            },
            BatchSize::LargeInput,
        )
    });
    group.finish();
}

criterion_group!(benches, bench_twobit);
criterion_main!(benches);
