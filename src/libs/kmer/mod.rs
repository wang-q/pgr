//! Canonical k-mer counting, profile generation, and repeat-run extraction.
//!
//! Native replacement for the FastK (`-p` / `-t` / `-p:<table>`) and Profex
//! (`-z`) steps of `pgr rept s-kmer` / `e-kmer`, with the `.pkt` compact
//! table cache. Design: `notes/design/kmer.md`.

pub mod count;
pub mod extract;
pub mod gc;
pub mod genomescope;
pub mod hist;
pub mod key;
pub mod khist;
pub mod nbinom;
pub mod profile;
pub mod qcheck;
pub mod quality;
pub mod supermer;

/// Sorted canonical k-mer table with parallel counts.
#[derive(Debug, Clone, Default)]
pub struct KmerTable {
    /// K-mer length (bp).
    pub k: usize,
    /// Packed canonical 2-bit k-mers (FastK byte layout), `key_bytes` bytes
    /// per entry, ascending and duplicate-free.
    pub keys: Vec<u8>,
    /// Dataset-wide counts, parallel to `keys`.
    pub counts: Vec<u32>,
}

impl KmerTable {
    /// Packed bytes per key (`(k+3)>>2`, FastK `kbyte`).
    pub fn key_bytes(&self) -> usize {
        self.k.div_ceil(4)
    }

    /// Key at entry `i`.
    pub fn key_at(&self, i: usize) -> key::Kmer {
        let kb = self.key_bytes();
        key::Kmer::from_bytes(self.k, &self.keys[i * kb..(i + 1) * kb])
    }
}

/// Emit `(position, canonical FastK byte key)` for every N-free k-mer
/// window of `seq`.
///
/// Bases are 2-bit encoded (A=0, C=1, G=2, T=3, case-insensitive); a window
/// containing N or any other non-ACGT base is skipped (FastK splits on gaps
/// and its profile has a 0 at such positions). `k` must be `<= Kmer::MAX_K`
/// (the caller validates); larger values emit nothing.
pub fn canonical_keys(seq: &[u8], k: usize, mut emit: impl FnMut(usize, &key::Kmer)) {
    let n = seq.len();
    if n < k || k > key::Kmer::MAX_K {
        return;
    }
    let codes = base_codes();
    let mut start = 0usize;
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
        let mut win = key::Kmer::from_bases(&seq[start..start + k], k).expect("N-free window");
        // Rolling canonical pair: the forward key advances at the 3' end
        // while the reverse complement advances at its 5' end (each new
        // forward base `x` prepends `3-x` to the rc), so the canonical key
        // costs one byte compare per window instead of a per-window rc.
        let mut win_rc = win.rc();
        emit(
            start,
            if canonical_le(&win, &win_rc) {
                &win
            } else {
                &win_rc
            },
        );
        for i in start + 1..=end - k {
            let x = codes[seq[i + k - 1] as usize] as u8;
            win.push_right(x);
            win_rc.push_left(3 - x);
            emit(
                i,
                if canonical_le(&win, &win_rc) {
                    &win
                } else {
                    &win_rc
                },
            );
        }
        start = end;
    }
}

/// FastK canonical comparison: only the first half of the packed bytes
/// matter, because forward and reverse complement are mirror-symmetric
/// (`count.c` compares `KMd2 = (KMER_BYTES+1)>>1` bytes).
fn canonical_le(a: &key::Kmer, b: &key::Kmer) -> bool {
    let half = a.key_bytes().div_ceil(2);
    a.to_bytes()[..half] <= b.to_bytes()[..half]
}

/// Base -> 2-bit code (0..3) or 4 (N / ambiguity); indexed by byte value.
pub fn base_codes() -> [u64; 256] {
    let mut codes = [4u64; 256];
    for b in *b"Aa" {
        codes[b as usize] = 0;
    }
    for b in *b"Cc" {
        codes[b as usize] = 1;
    }
    for b in *b"Gg" {
        codes[b as usize] = 2;
    }
    for b in *b"Tt" {
        codes[b as usize] = 3;
    }
    codes
}
