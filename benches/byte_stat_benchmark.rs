use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use pgr::libs::nt::{Nt, NT_VAL};
use pgr::libs::nt_simd::{self, SimdPath};
use rand::Rng;

// Mixed realistic-ish input: uppercase ACGT + lowercase soft-mask runs + N/IUPAC
// + occasional invalid bytes.
fn random_seq(rng: &mut impl Rng, len: usize) -> Vec<u8> {
    let bases = *b"ACGT";
    let mut seq: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();
    for _ in 0..len / 20 {
        let i = rng.random_range(0..len);
        seq[i] = match rng.random_range(0..4) {
            0 => b'n',
            1 => b'a',
            2 => b'M',
            _ => b'-',
        };
    }
    seq
}

fn scalar_count_valid(seq: &[u8]) -> usize {
    seq.iter()
        .filter(|&&b| !matches!(NT_VAL[b as usize] as u8, 4 | 255))
        .count()
}

fn scalar_count_n(seq: &[u8]) -> usize {
    seq.iter()
        .filter(|&&b| NT_VAL[b as usize] == Nt::N as usize)
        .count()
}

fn scalar_masked_bitmap(seq: &[u8], gap_only: bool) -> Vec<u32> {
    seq.chunks(32)
        .map(|chunk| {
            chunk.iter().enumerate().fold(0u32, |w, (k, &b)| {
                let is_n = NT_VAL[b as usize] == Nt::N as usize;
                let masked = if gap_only {
                    is_n
                } else {
                    is_n || b.is_ascii_lowercase()
                };
                w | ((masked as u32) << k)
            })
        })
        .collect()
}

fn bench_stats(c: &mut Criterion) {
    let mut rng = rand::rng();
    for (name, len) in [("1mb", 1_000_000usize), ("10mb", 10_000_000usize)] {
        let seq = random_seq(&mut rng, len);
        assert_eq!(
            nt_simd::count_valid(&seq),
            scalar_count_valid(&seq),
            "count_valid sanity"
        );
        assert_eq!(
            nt_simd::count_n(&seq),
            scalar_count_n(&seq),
            "count_n sanity"
        );
        assert_eq!(
            nt_simd::masked_bitmap(&seq, false),
            scalar_masked_bitmap(&seq, false),
            "bitmap sanity"
        );
        assert_eq!(
            nt_simd::masked_bitmap(&seq, true),
            scalar_masked_bitmap(&seq, true),
            "gap bitmap sanity"
        );

        let mut group = c.benchmark_group(format!("byte_stat_{name}"));

        group.bench_function("count_valid_scalar", |b| {
            b.iter_batched(
                || black_box(seq.clone()),
                |s| scalar_count_valid(&s),
                BatchSize::LargeInput,
            )
        });
        group.bench_function("count_valid_wide", |b| {
            b.iter_batched(
                || black_box(seq.clone()),
                |s| nt_simd::count_valid_with(SimdPath::Wide, &s),
                BatchSize::LargeInput,
            )
        });
        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") {
            group.bench_function("count_valid_avx2", |b| {
                b.iter_batched(
                    || black_box(seq.clone()),
                    |s| nt_simd::count_valid_with(SimdPath::Avx2, &s),
                    BatchSize::LargeInput,
                )
            });
        }

        group.bench_function("count_n_scalar", |b| {
            b.iter_batched(
                || black_box(seq.clone()),
                |s| scalar_count_n(&s),
                BatchSize::LargeInput,
            )
        });
        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") {
            group.bench_function("count_n_avx2", |b| {
                b.iter_batched(
                    || black_box(seq.clone()),
                    |s| nt_simd::count_n_with(SimdPath::Avx2, &s),
                    BatchSize::LargeInput,
                )
            });
        }

        group.bench_function("masked_bitmap_scalar", |b| {
            b.iter_batched(
                || black_box(seq.clone()),
                |s| scalar_masked_bitmap(&s, false),
                BatchSize::LargeInput,
            )
        });
        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") {
            group.bench_function("masked_bitmap_avx2", |b| {
                b.iter_batched(
                    || black_box(seq.clone()),
                    |s| nt_simd::masked_bitmap_with(SimdPath::Avx2, &s, false),
                    BatchSize::LargeInput,
                )
            });
        }
        group.finish();
    }
}

criterion_group!(benches, bench_stats);
criterion_main!(benches);
