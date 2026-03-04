use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use pgr::libs::pairmat::NamedMatrix;
use pgr::libs::clust::hier::{linkage_with_algo, Method, Algorithm};
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
    let mut group = c.benchmark_group("Hierarchical Clustering (Average Linkage)");
    
    // Test sizes
    // Keep sizes small enough for Primitive O(N^3) to finish reasonably fast
    // 100^3 = 1e6 ops (fast)
    // 500^3 = 1.25e8 ops (slow but ok)
    // 1000^3 = 1e9 ops (very slow)
    for size in [100, 200, 400].iter() {
        let matrix = create_random_matrix(*size);
        
        group.bench_with_input(BenchmarkId::new("Primitive", size), size, |b, &_| {
            b.iter(|| linkage_with_algo(black_box(&matrix), black_box(Method::Average), black_box(Algorithm::Primitive)))
        });
        
        group.bench_with_input(BenchmarkId::new("NN-Chain", size), size, |b, &_| {
            b.iter(|| linkage_with_algo(black_box(&matrix), black_box(Method::Average), black_box(Algorithm::NnChain)))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_hier);
criterion_main!(benches);
