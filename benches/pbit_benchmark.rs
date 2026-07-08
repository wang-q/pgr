//! Benchmarks for the pbit compressor/decompressor.
//!
//! Measures the core library APIs directly (not the CLI subprocess) to avoid
//! process-startup overhead dominating the measurements. Uses `tempfile` for
//! automatic cleanup and `rand::StdRng` for deterministic data.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pgr::libs::pbit::compressor::Compressor;
use pgr::libs::pbit::decompressor::Decompressor;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use std::fs;
use std::io::Write;
use std::path::Path;
use tempfile::TempDir;

/// Generate deterministic random DNA of the given length.
fn random_dna(len: usize, seed: u64) -> String {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..len)
        .map(|_| match rng.random_range(0u8..4) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        })
        .collect()
}

/// Introduce a SNP at every 100th position.
fn introduce_snps(seq: &str, seed: u64) -> String {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut out: Vec<char> = seq.chars().collect();
    for i in (0..out.len()).step_by(100) {
        let orig = out[i];
        let new = match orig {
            'A' => match rng.random_range(0u8..3) {
                0 => 'C',
                1 => 'G',
                _ => 'T',
            },
            'C' => match rng.random_range(0u8..3) {
                0 => 'A',
                1 => 'G',
                _ => 'T',
            },
            'G' => match rng.random_range(0u8..3) {
                0 => 'A',
                1 => 'C',
                _ => 'T',
            },
            _ => match rng.random_range(0u8..3) {
                0 => 'A',
                1 => 'C',
                _ => 'G',
            },
        };
        out[i] = new;
    }
    out.into_iter().collect()
}

/// Write a FASTA file with the given records.
fn write_fasta(path: &Path, records: &[(&str, &str)]) {
    let mut f = fs::File::create(path).unwrap();
    for (name, seq) in records {
        writeln!(f, ">{name}").unwrap();
        // Wrap at 60 chars per line (matches pgr's FASTA writer convention).
        for chunk in seq.as_bytes().chunks(60) {
            writeln!(f, "{}", std::str::from_utf8(chunk).unwrap()).unwrap();
        }
    }
}

/// Build a `.pbit` archive with `n_samples` samples derived from a 1 Mb
/// reference. Returns (temp_dir, pbit_path) so the caller keeps the temp dir
/// alive for the duration of the benchmark.
fn build_archive(n_samples: usize) -> (TempDir, std::path::PathBuf) {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let ref_seq = random_dna(1_000_000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);

    let pbit_path = temp.path().join("bench.pbit");
    let mut comp = Compressor::create(&pbit_path, ref_fa.to_str().unwrap(), 4096, 15, 18).unwrap();
    for i in 0..n_samples {
        let sample_fa = temp.path().join(format!("sample{i}.fa"));
        let sample_seq = introduce_snps(&ref_seq, 100 + i as u64);
        write_fasta(&sample_fa, &[("chr1", &sample_seq)]);
        comp.append_sample(&format!("sample{i}"), sample_fa.to_str().unwrap())
            .unwrap();
    }
    comp.finish().unwrap();
    (temp, pbit_path)
}

/// 1. Compress speed: 1 Mb reference + N samples (N = 1, 10, 100).
///
/// Measures the full `create` + `append_sample` * N + `finish` pipeline.
fn bench_compress(c: &mut Criterion) {
    let mut group = c.benchmark_group("Compress");
    // N=100 is slow (~100 Mb input), reduce sample size to keep runtime sane.
    group.sample_size(10);

    let ref_seq = random_dna(1_000_000, 42);

    for n_samples in [1usize, 10, 100].iter() {
        group.bench_with_input(
            BenchmarkId::new("1Mb_ref", n_samples),
            n_samples,
            |b, &n| {
                b.iter_with_setup(
                    || {
                        // Setup: create temp dir, write ref + N sample FASTAs.
                        let temp = TempDir::new().unwrap();
                        let ref_fa = temp.path().join("ref.fa");
                        write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
                        let sample_paths: Vec<(String, std::path::PathBuf)> = (0..n)
                            .map(|i| {
                                let p = temp.path().join(format!("sample{i}.fa"));
                                let s = introduce_snps(&ref_seq, 100 + i as u64);
                                write_fasta(&p, &[("chr1", &s)]);
                                (format!("sample{i}"), p)
                            })
                            .collect();
                        let pbit_path = temp.path().join("out.pbit");
                        (temp, ref_fa, sample_paths, pbit_path)
                    },
                    |(_temp, ref_fa, sample_paths, pbit_path)| {
                        let mut comp =
                            Compressor::create(&pbit_path, ref_fa.to_str().unwrap(), 4096, 15, 18)
                                .unwrap();
                        for (name, path) in &sample_paths {
                            comp.append_sample(name, path.to_str().unwrap()).unwrap();
                        }
                        comp.finish().unwrap();
                        black_box(pbit_path);
                    },
                )
            },
        );
    }
    group.finish();
}

/// 2. Decompress speed: full sample vs full contig vs interval slice.
///
/// Pre-builds an archive with 10 samples, then measures three extraction modes.
fn bench_decompress(c: &mut Criterion) {
    let (temp, pbit_path) = build_archive(10);

    let mut group = c.benchmark_group("Decompress");

    // a) get_sample: full extraction of one sample (1 Mb output).
    group.bench_function("get_sample_1Mb", |b| {
        b.iter_with_setup(
            || Decompressor::open(&pbit_path).unwrap(),
            |mut dec| {
                let mut out = Vec::new();
                dec.get_sample("sample0", &mut out).unwrap();
                black_box(out);
            },
        )
    });

    // b) get_contig: full contig from all samples (10 * 1 Mb output).
    group.bench_function("get_contig_full", |b| {
        b.iter_with_setup(
            || Decompressor::open(&pbit_path).unwrap(),
            |mut dec| {
                let mut out = Vec::new();
                dec.get_contig("chr1", None, None, "+", &mut out).unwrap();
                black_box(out);
            },
        )
    });

    // c) get_contig: interval slice [1000, 2000) from all samples.
    group.bench_function("get_contig_slice_1kb", |b| {
        b.iter_with_setup(
            || Decompressor::open(&pbit_path).unwrap(),
            |mut dec| {
                let mut out = Vec::new();
                dec.get_contig("chr1", Some(1000), Some(2000), "+", &mut out)
                    .unwrap();
                black_box(out);
            },
        )
    });

    group.finish();
    drop(temp);
}

/// 3. Decompressor::open: parse header + footer + indexes (no decompression).
fn bench_open(c: &mut Criterion) {
    let (temp, pbit_path) = build_archive(10);

    c.bench_function("Decompressor::open", |b| {
        b.iter(|| {
            let dec = Decompressor::open(&pbit_path).unwrap();
            black_box(dec);
        })
    });

    drop(temp);
}

/// 4. delta_cache hit rate: first get_contig is cold, second is warm.
///
/// Measures both in the same iteration to compare.
fn bench_cache_hit(c: &mut Criterion) {
    let (temp, pbit_path) = build_archive(10);

    let mut group = c.benchmark_group("CacheHit");

    // Cold: fresh Decompressor, no cache.
    group.bench_function("cold", |b| {
        b.iter_with_setup(
            || Decompressor::open(&pbit_path).unwrap(),
            |mut dec| {
                let mut out = Vec::new();
                dec.get_contig("chr1", Some(0), Some(1000), "+", &mut out)
                    .unwrap();
                black_box(out);
            },
        )
    });

    // Warm: same Decompressor, second call hits the delta cache.
    group.bench_function("warm", |b| {
        b.iter_with_setup(
            || {
                let mut dec = Decompressor::open(&pbit_path).unwrap();
                // Prime the cache with a first call.
                let mut out = Vec::new();
                dec.get_contig("chr1", Some(0), Some(1000), "+", &mut out)
                    .unwrap();
                dec
            },
            |mut dec| {
                let mut out = Vec::new();
                dec.get_contig("chr1", Some(0), Some(1000), "+", &mut out)
                    .unwrap();
                black_box(out);
            },
        )
    });

    group.finish();
    drop(temp);
}

criterion_group!(
    benches,
    bench_compress,
    bench_decompress,
    bench_open,
    bench_cache_hit
);
criterion_main!(benches);
