//! FastK-compatible k-mer byte keys.
//!
//! Layout matches FastK (`count.c` / `libfastk.c`): 2 bits per base, bytes
//! 5'->3', and within each byte the 5'-most base sits in the high 2 bits
//! (`kclip` keeps the high bits, verified against `FastK -t` table bytes).

/// FastK-encoded k-mer key: `key_bytes = ceil(k/4)` significant bytes in a
/// fixed-size array (no heap in the window-rolling hot path).
const MAX_K: usize = 128;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Kmer {
    k: usize,
    bytes: [u8; MAX_K / 4],
}

impl Kmer {
    /// Maximum k-mer length representable by the fixed-size value.
    pub const MAX_K: usize = crate::libs::kmer::key::MAX_K;

    /// New empty k-mer of length `k` (`k <= MAX_K`).
    pub fn new(k: usize) -> anyhow::Result<Self> {
        anyhow::ensure!(k > 0 && k <= MAX_K, "k must be in 1..={MAX_K}, got {k}");
        Ok(Self {
            k,
            bytes: [0; MAX_K / 4],
        })
    }

    /// K-mer length in bases.
    pub fn k(&self) -> usize {
        self.k
    }

    /// Encoding size in bytes (`(k+3)>>2`, same as FastK `kbyte`).
    pub fn key_bytes(&self) -> usize {
        (self.k + 3) >> 2
    }

    /// Encode the first `k` bases of `seq` (None if any base is not A/C/G/T).
    pub fn from_bases(seq: &[u8], k: usize) -> Option<Self> {
        if seq.len() < k || k > MAX_K {
            return None;
        }
        let mut km = Self::new(k).ok()?;
        for (i, &b) in seq[..k].iter().enumerate() {
            let c = match b {
                b'A' | b'a' => 0u8,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => return None,
            };
            km.set_base(i, c);
        }
        Some(km)
    }

    /// View the packed bytes (ascending = FastK table order).
    pub fn to_bytes(&self) -> &[u8] {
        &self.bytes[..self.key_bytes()]
    }

    /// 2-bit key as a big-endian integer (`2*k` significant bits).
    pub fn to_u128(&self) -> u128 {
        let mut x = 0u128;
        for &b in self.to_bytes() {
            x = (x << 8) | b as u128;
        }
        x >> (8 * self.key_bytes() - 2 * self.k)
    }

    /// Build from `key_bytes` packed bytes (same layout as `to_bytes`).
    pub fn from_bytes(k: usize, bytes: &[u8]) -> Self {
        let mut km = Self::new(k).expect("valid k");
        let n = km.key_bytes();
        debug_assert_eq!(bytes.len(), n);
        km.bytes[..n].copy_from_slice(&bytes[..n]);
        km
    }

    /// Base `i` (0-based from the 5' end) as a 2-bit code.
    pub fn base_at(&self, i: usize) -> u8 {
        debug_assert!(i < self.k);
        (self.bytes[i / 4] >> (2 * (3 - i % 4))) & 3
    }

    /// Packed byte `i` (0 = 5'-most byte).
    pub fn byte_at(&self, i: usize) -> u8 {
        self.bytes[i]
    }

    /// Advance the window by one base: drop the 5' base, append `base` at
    /// the 3' end (FastK rolling semantics).
    pub fn push_right(&mut self, base: u8) {
        let n = self.key_bytes();
        let s = 8 * n - 2 * self.k;
        let mut carry = 0u8;
        for i in (0..n).rev() {
            let old = self.bytes[i];
            self.bytes[i] = (old << 2) | carry;
            carry = old >> 6;
        }
        // The new base occupies bits s..s+1 of the packed window (the low
        // `s` bits are the FastK zero pad, not the integer low bits).
        let sh = s as u32;
        self.bytes[n - 1] = (self.bytes[n - 1] & !(0b11 << sh)) | ((base & 3) << sh);
    }

    /// Prepend `base` at the 5' end, dropping the 3' base.
    pub fn push_left(&mut self, base: u8) {
        let n = self.key_bytes();
        let s = 8 * n - 2 * self.k;
        let mut carry = 0u8;
        for i in 0..n {
            let old = self.bytes[i];
            self.bytes[i] = (old >> 2) | carry;
            carry = old << 6;
        }
        // After the shift the dropped 3' base lingers in bits s-2..s-1.
        if s >= 2 {
            let sh = (s - 2) as u32;
            self.bytes[n - 1] &= !(0b11 << sh);
        }
        self.bytes[0] = (self.bytes[0] & 0x3f) | ((base & 3) << 6);
    }

    /// Reverse complement (each base complemented and the order reversed).
    pub fn rc(&self) -> Self {
        let mut r = Self {
            k: self.k,
            bytes: [0; MAX_K / 4],
        };
        for i in 0..self.k {
            r.set_base(self.k - 1 - i, 3 - self.base_at(i));
        }
        r
    }

    /// Canonical key: the smaller of the forward key and its reverse
    /// complement (byte order, same representative as FastK).
    pub fn canonical(&self) -> Self {
        let r = self.rc();
        if self.to_bytes() <= r.to_bytes() {
            *self
        } else {
            r
        }
    }

    fn set_base(&mut self, i: usize, x: u8) {
        debug_assert!(i < self.k);
        let b = i / 4;
        let sh = 2 * (3 - i % 4);
        self.bytes[b] = (self.bytes[b] & !(0b11 << sh)) | ((x & 3) << sh);
    }
}

impl PartialOrd for Kmer {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Kmer {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.to_bytes().cmp(other.to_bytes())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encoding_matches_fastk_bytes() {
        // FastK `-t1 -k5` table bytes for ACGTACGTACGTACGT: acgta = 1b 00.
        let km = Kmer::from_bases(b"ACGT", 4).unwrap();
        assert_eq!(km.to_bytes(), &[0x1b]);
        let km = Kmer::from_bases(b"ACGTA", 5).unwrap();
        assert_eq!(km.to_bytes(), &[0x1b, 0x00]);
        let km = Kmer::from_bases(b"TACGT", 5).unwrap();
        assert_eq!(km.to_bytes(), &[0xc6, 0xc0]);
    }

    #[test]
    fn rc_reverses_and_complements() {
        let km = Kmer::from_bases(b"ACGTA", 5).unwrap();
        let rc = km.rc();
        assert_eq!(rc.to_bytes(), &[0xc6, 0xc0]); // TACGT
        assert_eq!(rc.rc(), km);
        // Palindromes map to themselves.
        let pal = Kmer::from_bases(b"ACGT", 4).unwrap();
        assert_eq!(pal.rc(), pal);
    }

    #[test]
    fn canonical_picks_smaller_strand() {
        let km = Kmer::from_bases(b"ACGTA", 5).unwrap();
        assert_eq!(km.canonical(), km);
        let rev = Kmer::from_bases(b"TACGT", 5).unwrap();
        assert_eq!(rev.canonical(), km);
    }

    /// Deterministic pseudo-random DNA block (same LCG as pgi tests).
    fn random_block(len: usize, seed: u64) -> Vec<u8> {
        let bases = *b"ACGT";
        let mut x = seed;
        (0..len)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                bases[(x >> 33) as usize & 3]
            })
            .collect()
    }

    #[test]
    fn rolling_matches_direct_encoding() {
        for &k in &[4usize, 5, 7, 8, 9, 21, 51, 81] {
            let seq = random_block(300, 1000 + k as u64);
            let mut win = Kmer::from_bases(&seq, k).unwrap();
            for i in 1..=seq.len() - k {
                win.push_right(crate::libs::kmer::base_codes()[seq[i + k - 1] as usize] as u8);
                let expect = Kmer::from_bases(&seq[i..], k).unwrap();
                assert_eq!(win, expect, "push_right k={k} at {i}");
            }
        }
    }

    #[test]
    fn push_left_matches_extension() {
        for &k in &[4usize, 5, 8, 9, 81] {
            let seq = random_block(200, 2000 + k as u64);
            // Extend leftwards one base at a time from the k-mer window at
            // position 100: prepend seq[99], seq[98], ...
            let mut win = Kmer::from_bases(&seq[100..100 + k], k).unwrap();
            for i in (0..100).rev() {
                win.push_left(crate::libs::kmer::base_codes()[seq[i] as usize] as u8);
                let expect = Kmer::from_bases(&seq[i..i + k], k).unwrap();
                assert_eq!(win, expect, "push_left k={k} at {i}");
            }
        }
    }

    #[test]
    fn from_bytes_roundtrip() {
        let km = Kmer::from_bases(b"ACGTACGTACGTTTTT", 15).unwrap();
        let back = Kmer::from_bytes(15, km.to_bytes());
        assert_eq!(back, km);
    }

    fn kmer_to_u128(km: &Kmer) -> u128 {
        let mut x = 0u128;
        for &b in km.to_bytes() {
            x = (x << 8) | b as u128;
        }
        x >> (8 * km.key_bytes() - 2 * km.k())
    }

    #[test]
    fn canonical_matches_existing_u128_semantics() {
        // The FastK byte key must select the same canonical representative
        // as the existing u128 `min(pack_kmer, rc_key)` for k <= 64.
        use crate::libs::nt;
        for &k in &[5usize, 21, 31, 63] {
            let seq = random_block(400, 3000 + k as u64);
            let mut win = Kmer::from_bases(&seq, k).unwrap();
            for i in 0..=seq.len() - k {
                if i > 0 {
                    win.push_right(crate::libs::kmer::base_codes()[seq[i + k - 1] as usize] as u8);
                }
                let fwd = nt::pack_kmer(&seq[i..i + k], k).unwrap();
                let expect = fwd.min(nt::rc_key(fwd, k));
                assert_eq!(kmer_to_u128(&win.canonical()), expect, "k={k} window {i}");
            }
        }
    }

    /// Emit canonical keys for every N-free window of `seq` (FastK splits
    /// on gaps), mirroring the future `kmer::canonical_keys` byte path.
    fn canonical_windows(seq: &[u8], k: usize) -> Vec<Kmer> {
        let codes = crate::libs::kmer::base_codes();
        let mut keys = Vec::new();
        let n = seq.len();
        let mut start = 0;
        while start + k <= n {
            while start < n && codes[seq[start] as usize] == 4 {
                start += 1;
            }
            if start + k > n {
                break;
            }
            let mut end = start;
            while end < n && codes[seq[end] as usize] != 4 {
                end += 1;
            }
            if end - start < k {
                start = end;
                continue;
            }
            let mut win = Kmer::from_bases(&seq[start..start + k], k).unwrap();
            keys.push(win.canonical());
            for i in start + 1..=end - k {
                win.push_right(codes[seq[i + k - 1] as usize] as u8);
                keys.push(win.canonical());
            }
            start = end;
        }
        keys
    }

    /// Read a gzipped test asset as text.
    fn read_gz(path: &str) -> String {
        use std::io::Read;
        let mut s = String::new();
        flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap())
            .read_to_string(&mut s)
            .unwrap();
        s
    }

    #[test]
    fn golden_table_matches_fastk() {
        // End-to-end check against the real FastK (`-t1`) table for the
        // shared input `tests/kmer/m1.fa.gz`: same canonical byte keys, same
        // ascending byte order, same counts (k=21/51/81).
        use crate::libs::ds::radix_sort::radix_sort_bytes_par;
        let fasta = read_gz(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/kmer/m1.fa.gz"));
        let seq: Vec<u8> = fasta.lines().skip(1).flat_map(|l| l.bytes()).collect();
        for &k in &[21usize, 51, 81] {
            let key_bytes = k.div_ceil(4);
            let keys = canonical_windows(&seq, k);
            // Sanity: the rolling path must agree with direct per-window
            // encoding (N-free windows only).
            let direct: Vec<Kmer> = (0..=seq.len() - k)
                .filter_map(|i| Kmer::from_bases(&seq[i..i + k], k))
                .map(|km| km.canonical())
                .collect();
            let mut a = keys
                .iter()
                .map(|km| km.to_bytes().to_vec())
                .collect::<Vec<_>>();
            let mut b = direct
                .iter()
                .map(|km| km.to_bytes().to_vec())
                .collect::<Vec<_>>();
            a.sort_unstable();
            b.sort_unstable();
            assert_eq!(a, b, "rolling vs direct disagree at k={k}");
            let mut raw: Vec<u8> = Vec::with_capacity(keys.len() * key_bytes);
            for km in &keys {
                raw.extend_from_slice(km.to_bytes());
            }
            let mut no_payload: Vec<()> = vec![(); keys.len()];
            radix_sort_bytes_par(&mut raw, key_bytes, &mut no_payload);
            let mut got: Vec<(String, u32)> = Vec::new();
            for chunk in raw.chunks(key_bytes) {
                let km = Kmer::from_bytes(k, chunk);
                let s: String = (0..k)
                    .map(|i| b"acgt"[km.base_at(i) as usize] as char)
                    .collect();
                match got.last_mut() {
                    Some((last, c)) if *last == s => *c += 1,
                    _ => got.push((s, 1)),
                }
            }
            let golden = read_gz(&format!(
                "{}/tests/kmer/fastk_k{k}.golden.gz",
                env!("CARGO_MANIFEST_DIR")
            ));
            let expect: Vec<(String, u32)> = golden
                .lines()
                .map(|l| {
                    let (s, c) = l.split_once('\t').unwrap();
                    (s.to_string(), c.parse().unwrap())
                })
                .collect();
            assert_eq!(got.len(), expect.len(), "k={k} entry count");
            for (i, (g, e)) in got.iter().zip(&expect).enumerate() {
                assert_eq!(g, e, "k={k} entry {i}");
            }
        }
    }
}
