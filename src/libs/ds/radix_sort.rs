//! In-place MSD radix sort for `(u128 key, payload)` pairs.
//!
//! Port of Gene Myers' American-flag radix sort (FastGA `MSDsort.c`): at each
//! key byte the records are permuted into their buckets by following
//! displacement cycles, so no auxiliary arrays are needed and each record
//! moves once per byte level. Small segments fall back to insertion sort.
//!
//! The sort is deterministic: the same input always yields the same output
//! order (including the payload order within equal keys, which the algorithm
//! does not promise to stabilize).

/// Small-segment cutoff: below this record count a comparison sort is used.
const SMALL: usize = 16;
/// Record count above which the parallel variant distributes by top byte.
const PAR_SMALL: usize = 1 << 18;

/// Sort `keys` ascending (low `key_bits` bits significant) and permute
/// `payloads` to match. MSD radix, big-endian byte order, in place.
pub fn radix_sort_u128<T: Copy>(keys: &mut [u128], payloads: &mut [T], key_bits: u32) {
    assert_eq!(keys.len(), payloads.len());
    if keys.len() < 2 {
        return;
    }
    let key_bytes = key_bits.div_ceil(8) as usize;
    let mut stack = Vec::new();
    msd(keys, payloads, key_bytes - 1, &mut stack);
}

/// Sort a `u128` key array in place (no payload to permute).
pub fn radix_sort_keys(keys: &mut [u128], key_bits: u32) {
    let mut no_payload: Vec<()> = vec![(); keys.len()];
    radix_sort_u128(keys, &mut no_payload, key_bits);
}

/// Parallel variant of [`radix_sort_u128`]: distributes the records into 256
/// buckets by the most significant key byte in place, then sorts each bucket
/// in parallel (rayon) with the same MSD radix sort.
pub fn radix_sort_u128_par<T: Copy + Send + Sync>(
    keys: &mut [u128],
    payloads: &mut [T],
    key_bits: u32,
) {
    assert_eq!(keys.len(), payloads.len());
    let n = keys.len();
    if n < 2 {
        return;
    }
    if n < PAR_SMALL {
        radix_sort_u128(keys, payloads, key_bits);
        return;
    }
    let key_bytes = key_bits.div_ceil(8) as usize;
    let mut byte = key_bytes - 1;
    // Skip leading bytes identical across the whole array (common prefix).
    while byte > 0 {
        let first = key_byte(keys[0], byte);
        if keys[1..].iter().any(|&k| key_byte(k, byte) != first) {
            break;
        }
        byte -= 1;
    }
    if keys[1..]
        .iter()
        .all(|&k| key_byte(k, byte) == key_byte(keys[0], byte))
    {
        return; // every key is equal
    }
    let mut stack = Vec::new();
    let offsets = partition_at(keys, payloads, byte, &mut stack);
    if byte > 0 {
        sort_buckets_par(keys, payloads, &offsets, 0, 256, byte - 1, 0);
    }
}

/// Byte `byte` (0 = least significant) of `key`.
#[inline]
fn key_byte(key: u128, byte: usize) -> usize {
    ((key >> (8 * byte)) & 0xff) as usize
}

/// Sort a small segment with a comparison sort (payloads follow keys).
fn insertion_sort<T: Copy>(keys: &mut [u128], payloads: &mut [T]) {
    for i in 1..keys.len() {
        let k = keys[i];
        let p = payloads[i];
        let mut j = i;
        while j > 0 && keys[j - 1] > k {
            keys[j] = keys[j - 1];
            payloads[j] = payloads[j - 1];
            j -= 1;
        }
        keys[j] = k;
        payloads[j] = p;
    }
}

/// Sort `keys[..]` by the bytes `byte`, `byte-1`, ..., 0 of their `u128` keys.
fn msd<T: Copy>(keys: &mut [u128], payloads: &mut [T], byte: usize, stack: &mut Vec<usize>) {
    let n = keys.len();
    if n <= SMALL {
        insertion_sort(keys, payloads);
        return;
    }

    // Skip leading bytes that are identical across the whole segment (the
    // caller has already partitioned on every byte above `byte`).
    let mut b = byte;
    while b > 0 {
        let first = key_byte(keys[0], b);
        if keys[1..].iter().any(|&k| key_byte(k, b) != first) {
            break;
        }
        b -= 1;
    }
    let first = key_byte(keys[0], b);
    if keys[1..].iter().all(|&k| key_byte(k, b) == first) {
        return; // every key in this segment is equal
    }

    let offsets = partition_at(keys, payloads, b, stack);
    if b == 0 {
        return;
    }
    for v in 0..256 {
        let c = offsets[v + 1] - offsets[v];
        if c > 1 {
            let (s, e) = (offsets[v], offsets[v + 1]);
            msd(&mut keys[s..e], &mut payloads[s..e], b - 1, stack);
        }
    }
}

/// Partition `keys[..]` in place by `byte` (counting pass + American-flag
/// cycle permutation) and return the bucket start offsets (`offsets[256]` is
/// the segment length).
fn partition_at<T: Copy>(
    keys: &mut [u128],
    payloads: &mut [T],
    byte: usize,
    stack: &mut Vec<usize>,
) -> [usize; 257] {
    let n = keys.len();
    let mut counts = [0usize; 256];
    for &k in &keys[..n] {
        counts[key_byte(k, byte)] += 1;
    }
    let mut offsets = [0usize; 256];
    let mut cum = 0usize;
    for (o, &c) in offsets.iter_mut().zip(counts.iter()) {
        *o = cum;
        cum += c;
    }
    let mut boundaries = [0usize; 257];
    boundaries[..256].copy_from_slice(&offsets);
    boundaries[256] = n;

    // Permute each bucket into its region by following displacement cycles
    // (American-flag sort; every record moves once per byte level).
    let mut next = offsets;
    for v in 0..256 {
        let end = offsets[v] + counts[v];
        while next[v] < end {
            let t = key_byte(keys[next[v]], byte);
            if t == v {
                next[v] += 1;
                continue;
            }
            stack.clear();
            stack.push(next[v]);
            let mut t = t;
            loop {
                if t == v {
                    next[v] += 1;
                    break;
                }
                let mut u = next[t];
                while key_byte(keys[u], byte) == t {
                    u += 1;
                }
                next[t] = u + 1;
                stack.push(u);
                t = key_byte(keys[u], byte);
            }
            let last = stack[stack.len() - 1];
            let lk = keys[last];
            let lp = payloads[last];
            for k in (1..stack.len()).rev() {
                keys[stack[k]] = keys[stack[k - 1]];
                payloads[stack[k]] = payloads[stack[k - 1]];
            }
            keys[stack[0]] = lk;
            payloads[stack[0]] = lp;
        }
    }
    boundaries
}

/// Sort the bucket ranges `[b_lo, b_hi)` in parallel by recursively splitting
/// the bucket index range at bucket boundaries (`rayon::join`). `base` is the
/// absolute index of `keys[0]` so `offsets` stay absolute while the slices
/// shrink.
fn sort_buckets_par<T: Copy + Send + Sync>(
    keys: &mut [u128],
    payloads: &mut [T],
    offsets: &[usize; 257],
    b_lo: usize,
    b_hi: usize,
    byte: usize,
    base: usize,
) {
    if b_hi - b_lo == 1 {
        let (s, e) = (offsets[b_lo] - base, offsets[b_lo + 1] - base);
        if e - s > 1 {
            let mut stack = Vec::new();
            msd(&mut keys[s..e], &mut payloads[s..e], byte, &mut stack);
        }
        return;
    }
    let mid = (b_lo + b_hi) / 2;
    let cut = offsets[mid] - base;
    let (k1, k2) = keys.split_at_mut(cut);
    let (p1, p2) = payloads.split_at_mut(cut);
    rayon::join(
        || sort_buckets_par(k1, p1, offsets, b_lo, mid, byte, base),
        || sort_buckets_par(k2, p2, offsets, mid, b_hi, byte, base + cut),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn radix_sorts_keys_and_payloads() {
        let mut keys = vec![5u128, 3, 8, 1, 3, 0];
        let mut payloads = vec![0u32, 1, 2, 3, 4, 5];
        radix_sort_u128(&mut keys, &mut payloads, 40);
        assert_eq!(keys, vec![0, 1, 3, 3, 5, 8]);
        assert_eq!(payloads, vec![5, 3, 1, 4, 0, 2]);
    }

    #[test]
    fn radix_matches_sort_unstable() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        use rand::{Rng, SeedableRng};
        for &bits in &[8u32, 20, 40, 80, 112] {
            for &n in &[0usize, 1, 2, 17, 500, 5000] {
                for _ in 0..3 {
                    let mut keys: Vec<u128> = (0..n)
                        .map(|_| rng.random_range(0..(1u128 << bits.min(80))))
                        .collect();
                    let mut payloads: Vec<u32> = (0..n as u32).collect();
                    // The American-flag sort is not stable within equal keys,
                    // so compare the sorted (key, payload) multiset instead of
                    // a positional match.
                    let mut expect: Vec<(u128, u32)> =
                        keys.iter().copied().zip(payloads.iter().copied()).collect();
                    expect.sort_unstable();

                    radix_sort_u128(&mut keys, &mut payloads, bits);
                    assert!(
                        keys.windows(2).all(|w| w[0] <= w[1]),
                        "keys not ascending at {bits} bits, n={n}"
                    );
                    let mut got: Vec<(u128, u32)> =
                        keys.iter().copied().zip(payloads.iter().copied()).collect();
                    got.sort_unstable();
                    assert_eq!(got, expect, "payloads lost/mixed at {bits} bits, n={n}");
                }
            }
        }
    }

    #[test]
    fn radix_handles_duplicates_and_skew() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        use rand::{Rng, SeedableRng};
        // Mostly duplicate keys: exercise common-prefix skipping and cycles.
        for _ in 0..20 {
            let n = 2000usize;
            let mut keys: Vec<u128> = (0..n)
                .map(|_| {
                    if rng.random_range(0..4) == 0 {
                        rng.random_range(0..8)
                    } else {
                        rng.random_range(0..(1u128 << 40))
                    }
                })
                .collect();
            let mut payloads: Vec<u64> = (0..n as u64).collect();
            let mut expect: Vec<(u128, u64)> =
                keys.iter().copied().zip(payloads.iter().copied()).collect();
            expect.sort_unstable();
            radix_sort_u128(&mut keys, &mut payloads, 40);
            assert!(keys.windows(2).all(|w| w[0] <= w[1]));
            let mut got: Vec<(u128, u64)> =
                keys.iter().copied().zip(payloads.iter().copied()).collect();
            got.sort_unstable();
            assert_eq!(got, expect);
        }
    }

    #[test]
    fn radix_sort_keys_matches() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(99);
        use rand::{Rng, SeedableRng};
        let mut keys: Vec<u128> = (0..10_000u32)
            .map(|_| rng.random_range(0..(1u128 << 112)))
            .collect();
        let mut expect = keys.clone();
        expect.sort_unstable();
        radix_sort_keys(&mut keys, 112);
        assert_eq!(keys, expect);
    }

    #[test]
    fn radix_par_matches_sequential() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1234);
        use rand::{Rng, SeedableRng};
        // Above PAR_SMALL so the top-byte distribution path runs.
        let n = (1 << 19) + 17usize;
        let mut keys: Vec<u128> = (0..n).map(|_| rng.random_range(0..(1u128 << 80))).collect();
        let mut payloads: Vec<(u32, u32, u8)> = (0..n as u32)
            .map(|i| (i % 7, i / 7, (i & 1) as u8))
            .collect();
        let mut expect: Vec<(u128, (u32, u32, u8))> =
            keys.iter().copied().zip(payloads.iter().copied()).collect();
        expect.sort_unstable();

        radix_sort_u128_par(&mut keys, &mut payloads, 80);
        assert!(keys.windows(2).all(|w| w[0] <= w[1]));
        let mut got: Vec<(u128, (u32, u32, u8))> =
            keys.iter().copied().zip(payloads.iter().copied()).collect();
        got.sort_unstable();
        assert_eq!(got, expect);
    }
}
