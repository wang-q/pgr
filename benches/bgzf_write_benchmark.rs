use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pgr::libs::bgzf::{BgzfWriter, ParallelBgzfWriter};
use std::io::{sink, Write};
use std::sync::OnceLock;

/// Deterministic ~50 MB FASTA (same generator as `bgzf_benchmark.rs`).
fn fasta_50m() -> &'static Vec<u8> {
    static FA: OnceLock<Vec<u8>> = OnceLock::new();
    FA.get_or_init(|| {
        let bases = b"ACGT";
        let mut x = 0x1234_5678_9abc_def0u64;
        let mut fa = Vec::with_capacity(50 * (1_000_000 + 1_000_000 / 80 + 20));
        for s in 0..50 {
            fa.extend_from_slice(format!(">seq{s}\n").as_bytes());
            let mut col = 0usize;
            for _ in 0..1_000_000 {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                fa.push(bases[(x >> 33) as usize & 3]);
                col += 1;
                if col == 80 {
                    fa.push(b'\n');
                    col = 0;
                }
            }
            if col != 0 {
                fa.push(b'\n');
            }
        }
        fa
    })
}

fn bench_write_group(c: &mut Criterion) {
    let fa = fasta_50m();
    let mut group = c.benchmark_group("bgzf_write_50m");
    group.sample_size(10);

    group.bench_function("single_threaded", |b| {
        b.iter(|| {
            let mut w = BgzfWriter::new(sink()).expect("writer");
            w.write_all(black_box(fa)).expect("write");
            w.finish().expect("finish");
        })
    });

    for workers in [1usize, 2, 3, 4, 6, 8] {
        group.bench_with_input(
            BenchmarkId::new("parallel", workers),
            &workers,
            |b, &workers| {
                b.iter(|| {
                    let mut w = ParallelBgzfWriter::new(sink(), workers).expect("writer");
                    w.write_all(black_box(fa)).expect("write");
                    w.finish().expect("finish");
                })
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_write_group);
criterion_main!(benches);
