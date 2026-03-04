use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pgr::libs::clust::hier::{linkage_with_algo, Algorithm, Method};
use pgr::libs::pairmat::NamedMatrix;
use rand::Rng;

fn create_random_matrix(size: usize) -> NamedMatrix {
    let mut rng = rand::rng();
    let names: Vec<String> = (0..size).map(|i| i.to_string()).collect();
    let mut matrix = NamedMatrix::new(names);

    for i in 0..size {
        for j in (i + 1)..size {
            let val: f32 = rng.random();
            matrix.set(i, j, val);
        }
    }
    matrix
}

fn bench_hier(c: &mut Criterion) {
    // 1. Primitive vs NN-Chain (Small N)
    // Demonstrate the O(N^3) vs O(N^2) difference
    let mut group = c.benchmark_group("Algo Comparison");
    for size in [100, 200, 400].iter() {
        let matrix = create_random_matrix(*size);

        group.bench_with_input(
            BenchmarkId::new("Primitive (Average)", size),
            size,
            |b, &_| {
                b.iter(|| {
                    linkage_with_algo(
                        black_box(&matrix),
                        black_box(Method::Average),
                        black_box(Algorithm::Primitive),
                    )
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("NN-Chain (Average)", size),
            size,
            |b, &_| {
                b.iter(|| {
                    linkage_with_algo(
                        black_box(&matrix),
                        black_box(Method::Average),
                        black_box(Algorithm::NnChain),
                    )
                })
            },
        );
    }
    group.finish();

    // 2. Scalability (Large N, NN-Chain only)
    // Assess performance on larger datasets
    let mut group = c.benchmark_group("Scalability (NN-Chain)");
    group.sample_size(10); // Reduce sample size for slow benchmarks
    for size in [1000, 2000, 4000].iter() {
        let matrix = create_random_matrix(*size);

        group.bench_with_input(BenchmarkId::new("Average", size), size, |b, &_| {
            b.iter(|| {
                linkage_with_algo(
                    black_box(&matrix),
                    black_box(Method::Average),
                    black_box(Algorithm::NnChain),
                )
            })
        });

        // Also test Ward to verify optimization scalability
        group.bench_with_input(BenchmarkId::new("Ward", size), size, |b, &_| {
            b.iter(|| {
                linkage_with_algo(
                    black_box(&matrix),
                    black_box(Method::Ward),
                    black_box(Algorithm::NnChain),
                )
            })
        });
    }
    group.finish();

    // 3. Method Comparison (Fixed N=1000)
    // Compare relative cost of different linkage criteria
    let mut group = c.benchmark_group("Method Comparison (N=1000)");
    let size = 1000;
    let matrix = create_random_matrix(size);

    for method in [
        Method::Single,
        Method::Complete,
        Method::Average,
        Method::Weighted,
        Method::Ward,
        // Centroid/Median might fallback to Primitive if not Reducible?
        // Actually NN-Chain works for them if we ignore reducibility issues or if they are reducible in geometric sense?
        // Wait, Centroid/Median are NOT reducible, so pgr might fallback to Primitive automatically?
        // Let's check: Algorithm::Auto does that. But here we force NN-Chain or use Auto?
        // Let's use Auto to see real-world performance.
    ]
    .iter()
    {
        group.bench_with_input(
            BenchmarkId::new(format!("{:?}", method), size),
            method,
            |b, &m| {
                b.iter(|| {
                    linkage_with_algo(black_box(&matrix), black_box(m), black_box(Algorithm::Auto))
                })
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_hier);
criterion_main!(benches);
