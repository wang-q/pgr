use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::pgi::build::read_fasta;
use pgr::libs::syncmer::{syncmer_dna, SyncmerParams};

const GENOME: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz");

/// Splitmix64-style odd factor, mirroring `libs::syncmer::hash_factor`.
fn hash_factor(seed: u64) -> u64 {
    let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    (z ^ (z >> 31)) | 1
}

/// Rolling canonical s-mer hashes only (no sliding-window minimum), mirroring
/// `dna_canonical_hashes`. Returns an accumulator over all canonical hashes
/// so the compiler cannot fold the loop away.
fn rolling_hashes_only(seq: &[u8], params: &SyncmerParams) -> u64 {
    let k = params.smer;
    let n = seq.len();
    if n < k {
        return 0;
    }
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let shift: u32 = (64 - 2 * k) as u32;
    let factor = hash_factor(params.seed);
    let pattern_rc: [u64; 4] = std::array::from_fn(|i| ((3 - i) as u64) << (2 * (k - 1)));
    let encode_base = |b: u8| {
        let v = pgr::libs::nt::NT_VAL[b as usize];
        if v <= 3 {
            v as u64
        } else {
            0
        }
    };
    let mut h: u64 = 0;
    let mut h_rc: u64 = 0;
    for &byte in seq.iter().take(k) {
        let b = encode_base(byte);
        h = (h << 2) | b;
        h_rc = (h_rc >> 2) | pattern_rc[b as usize];
    }
    let mut acc = h.wrapping_mul(factor) >> shift;
    for &byte in seq.iter().skip(k) {
        let b = encode_base(byte);
        h = ((h << 2) & mask) | b;
        h_rc = (h_rc >> 2) | pattern_rc[b as usize];
        let (hf, hr) = (
            h.wrapping_mul(factor) >> shift,
            h_rc.wrapping_mul(factor) >> shift,
        );
        acc = acc.wrapping_add(hf.min(hr));
    }
    acc
}

fn bench_syncmer(c: &mut Criterion) {
    // pgi build defaults: smer=8, window=5 (span = 12).
    let params = SyncmerParams {
        smer: 8,
        window: 5,
        seed: 7,
    };
    let seqs: Vec<Vec<u8>> = read_fasta(GENOME)
        .unwrap()
        .into_iter()
        .map(|(_, s)| s)
        .collect();
    let mut group = c.benchmark_group("syncmer_dna");
    group.sample_size(10);
    group.bench_function("mg1655_s8_w5", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for seq in &seqs {
                let syncs = syncmer_dna(seq, &params).unwrap();
                total += syncs.len();
            }
            black_box(total)
        })
    });
    group.bench_function("rolling_hashes_only_s8_w5", |b| {
        b.iter(|| {
            let mut total = 0u64;
            for seq in &seqs {
                total = total.wrapping_add(rolling_hashes_only(seq, &params));
            }
            black_box(total)
        })
    });
    group.finish();
}

criterion_group!(benches, bench_syncmer);
criterion_main!(benches);
