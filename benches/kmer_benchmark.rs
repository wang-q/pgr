//! Benchmarks for the native k-mer repeat pipeline: canonical count-table
//! build and per-position profile generation on MG1655
//! (`tests/genome/mg1655.fa.gz`, 4.6 Mb).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::kmer;
use pgr::libs::pgi::build::read_fasta;

const GENOME: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz");
const LIB: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/pgr/tncentral.fa.gz");

fn bench_build_and_profiles(c: &mut Criterion) {
    let k = 17usize;
    let seqs: Vec<Vec<u8>> = read_fasta(GENOME)
        .unwrap()
        .into_iter()
        .map(|(_, s)| s)
        .collect();
    let lib_seqs: Vec<Vec<u8>> = read_fasta(LIB)
        .unwrap()
        .into_iter()
        .map(|(_, s)| s)
        .collect();

    let mut group = c.benchmark_group("kmer");
    group.sample_size(10);
    group.bench_function("count_mg1655", |b| {
        b.iter(|| {
            let table = kmer::count::build_table(&seqs, k).unwrap();
            black_box(table.counts.len())
        })
    });
    group.bench_function("self_profiles_mg1655", |b| {
        let table = kmer::count::build_table(&seqs, k).unwrap();
        b.iter(|| {
            let profiles = kmer::profile::self_profiles(&seqs, k, &table);
            black_box(profiles[0].len())
        })
    });
    group.bench_function("relative_profiles_mg1655", |b| {
        let table = kmer::count::build_table(&lib_seqs, k).unwrap();
        b.iter(|| {
            let profiles = kmer::profile::relative_profiles(&seqs, k, &table);
            black_box(profiles[0].len())
        })
    });
    group.bench_function("canonical_keys_only_mg1655", |b| {
        b.iter(|| {
            let mut n = 0usize;
            for seq in &seqs {
                kmer::canonical_keys(seq, k, |_, _| n += 1);
            }
            black_box(n)
        })
    });
    group.finish();
}

// Isolated lookup comparison kept as a negative result (2026-08-09): on the
// 73 MB MG1655 table, both global partition_point and a prefix-bucket index
// take ~1.1 s for 4.6M lookups -- random-access latency dominates, not the
// number of comparisons. The production path switched to sort+merge instead
// (see notes/benchmarks/bench-profile-hotspots.md).
struct BenchPrefixIndex {
    key_bytes: usize,
    offsets: Vec<u32>,
}

impl BenchPrefixIndex {
    fn new(keys: &[u8], key_bytes: usize) -> Self {
        let n_buckets = 1usize << 16;
        let mut offsets = vec![0u32; n_buckets + 1];
        let n = keys.len() / key_bytes;
        let mut prev = 0usize;
        for i in 0..n {
            let bucket = prefix_bytes(keys, key_bytes, i);
            while prev < bucket {
                offsets[prev + 1] = i as u32;
                prev += 1;
            }
        }
        for b in prev..n_buckets {
            offsets[b + 1] = n as u32;
        }
        Self { key_bytes, offsets }
    }

    fn lookup(&self, keys: &[u8], counts: &[u32], key: &pgr::libs::kmer::key::Kmer) -> Option<u32> {
        let kb = self.key_bytes;
        let bucket = prefix_bytes(key.to_bytes(), kb, 0);
        let start = self.offsets[bucket] as usize;
        let end = self.offsets[bucket + 1] as usize;
        let mut lo = start;
        let mut hi = end;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if &keys[mid * kb..(mid + 1) * kb] < key.to_bytes() {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if lo < end && &keys[lo * kb..(lo + 1) * kb] == key.to_bytes() {
            Some(counts[lo])
        } else {
            None
        }
    }
}

/// 16-bit prefix of packed record `i` (high bytes, zero padded on the left
/// when `key_bytes < 2`).
fn prefix_bytes(keys: &[u8], key_bytes: usize, i: usize) -> usize {
    let mut b = 0usize;
    for j in 0..key_bytes.min(2) {
        b = (b << 8) | keys[i * key_bytes + j] as usize;
    }
    b << ((2 - key_bytes.min(2)) * 8)
}

fn bench_lookups(c: &mut Criterion) {
    let k = 17usize;
    let seqs: Vec<Vec<u8>> = read_fasta(GENOME)
        .unwrap()
        .into_iter()
        .map(|(_, s)| s)
        .collect();
    let table = kmer::count::build_table(&seqs, k).unwrap();
    let key_bytes = table.key_bytes();
    let mut windows: Vec<pgr::libs::kmer::key::Kmer> = Vec::new();
    for seq in &seqs {
        kmer::canonical_keys(seq, k, |_, key| windows.push(key));
    }
    let index = BenchPrefixIndex::new(&table.keys, key_bytes);
    let mut group = c.benchmark_group("kmer_lookup");
    group.bench_function("global_partition_point", |b| {
        b.iter(|| {
            let mut hits = 0usize;
            for key in &windows {
                let mut lo = 0usize;
                let mut hi = table.counts.len();
                while lo < hi {
                    let mid = (lo + hi) / 2;
                    if &table.keys[mid * key_bytes..(mid + 1) * key_bytes] < key.to_bytes() {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                if lo < table.counts.len()
                    && &table.keys[lo * key_bytes..(lo + 1) * key_bytes] == key.to_bytes()
                {
                    hits += 1;
                }
            }
            black_box(hits)
        })
    });
    group.bench_function("prefix_index", |b| {
        b.iter(|| {
            let mut hits = 0usize;
            for key in &windows {
                if index.lookup(&table.keys, &table.counts, key).is_some() {
                    hits += 1;
                }
            }
            black_box(hits)
        })
    });
    group.finish();
}

criterion_group!(benches, bench_build_and_profiles);
criterion_group!(benches_lookup, bench_lookups);
criterion_main!(benches, benches_lookup);
