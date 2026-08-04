//! Benchmark: `rg merge` mapping construction over large single-chromosome
//! `.rg` files.
//!
//! Two workload shapes: disjoint ranges (the dedup-heavy case that used to
//! be O(n²) via `Vec::contains`) and clusters of overlapping ranges (the
//! COITree query + union-find path).

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pgr::libs::runlist::rg_merge_mapping;

/// Write an `.rg` file with `n` ranges on `chr1` and return its path.
fn write_rg(dir: &std::path::Path, n: usize, clustered: bool) -> String {
    let mut content = String::with_capacity(n * 18);
    if clustered {
        for i in 0..n {
            let start = (i / 4) * 300 + (i % 4) * 60 + 1;
            content.push_str(&format!("chr1:{}-{}\n", start, start + 100));
        }
    } else {
        let mut pos = 1i32;
        for _ in 0..n {
            content.push_str(&format!("chr1:{}-{}\n", pos, pos + 100));
            pos += 150;
        }
    }
    let path = dir.join("in.rg");
    std::fs::write(&path, content).unwrap();
    path.to_string_lossy().into_owned()
}

fn bench_merge(c: &mut Criterion) {
    let mut group = c.benchmark_group("rg_merge");
    group.sample_size(10);
    for &n in &[10_000usize, 50_000] {
        for &clustered in &[false, true] {
            let dir = tempfile::tempdir().unwrap();
            let path = write_rg(dir.path(), n, clustered);
            let label = format!("{}-{}", if clustered { "clustered" } else { "disjoint" }, n);
            group.bench_with_input(BenchmarkId::new("mapping", label), &path, |b, p| {
                b.iter(|| {
                    black_box(rg_merge_mapping(std::slice::from_ref(p), 0.95).unwrap());
                });
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_merge);
criterion_main!(benches);
