use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use pgr::libs::poa::align::{
    AlignmentEngine, AlignmentParams, AlignmentType, ScalarAlignmentEngine,
};
use pgr::libs::poa::graph::PoaGraph;
use pgr::libs::poa::simd::{SimdAlignmentEngine, SimdPath};
use pgr::libs::poa::Poa;
use rand::Rng;

// Consensus-like scenario: a set of similar sequences (10% mutations) forms
// the POA graph, then one more sequence is aligned to it.
fn make_seqs(rng: &mut impl Rng, n: usize, len: usize) -> Vec<Vec<u8>> {
    let bases = *b"ACGT";
    let ancestor: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();
    (0..n)
        .map(|_| {
            ancestor
                .iter()
                .map(|&b| {
                    if rng.random_bool(0.1) {
                        bases[rng.random_range(0..4)]
                    } else {
                        b
                    }
                })
                .collect()
        })
        .collect()
}

fn bench_align(c: &mut Criterion) {
    let mut rng = rand::rng();
    for (name, seq_len) in [("short_120bp", 120usize), ("long_600bp", 600usize)] {
        let seqs = make_seqs(&mut rng, 21, seq_len);
        let mut poa = Poa::new(AlignmentParams::default(), AlignmentType::Global);
        for s in &seqs[..20] {
            poa.add_sequence(s);
        }
        let graph: &PoaGraph = poa.graph();
        let seq = &seqs[20];
        let scalar = ScalarAlignmentEngine::new(AlignmentParams::default(), AlignmentType::Global);
        let simd = SimdAlignmentEngine::new(AlignmentParams::default(), AlignmentType::Global);

        // Same score/path sanity check across all engines/paths.
        let a = scalar.align(seq, graph);
        let b = simd.align(seq, graph);
        let w = simd.align_with(SimdPath::Wide, seq, graph);
        assert_eq!(a.score, b.score);
        assert_eq!(a.path, b.path);
        assert_eq!(a.score, w.score);
        assert_eq!(a.path, w.path);

        let mut group = c.benchmark_group(format!("poa_align_{name}"));
        group.bench_function("scalar", |b| {
            b.iter_batched(
                || black_box(seq.clone()),
                |s| scalar.align(&s, black_box(graph)),
                BatchSize::SmallInput,
            )
        });
        group.bench_function("wide", |b| {
            b.iter_batched(
                || black_box(seq.clone()),
                |s| simd.align_with(SimdPath::Wide, &s, black_box(graph)),
                BatchSize::SmallInput,
            )
        });
        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") {
            group.bench_function("avx2", |b| {
                b.iter_batched(
                    || black_box(seq.clone()),
                    |s| simd.align_with(SimdPath::Avx2, &s, black_box(graph)),
                    BatchSize::SmallInput,
                )
            });
        }
        group.bench_function("simd", |b| {
            b.iter_batched(
                || black_box(seq.clone()),
                |s| simd.align(&s, black_box(graph)),
                BatchSize::SmallInput,
            )
        });
        group.finish();
    }
}

criterion_group!(benches, bench_align);
criterion_main!(benches);
