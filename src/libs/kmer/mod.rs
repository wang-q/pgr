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
pub mod khist;
pub mod nbinom;
pub mod profile;
pub mod qcheck;
pub mod quality;

/// Sorted canonical k-mer table with parallel counts.
#[derive(Debug, Clone, Default)]
pub struct KmerTable {
    /// K-mer length (bp).
    pub k: usize,
    /// Canonical 2-bit k-mers (forward vs reverse-complement, smaller),
    /// ascending and duplicate-free.
    pub keys: Vec<u128>,
    /// Dataset-wide counts, parallel to `keys`.
    pub counts: Vec<u32>,
}

/// Emit `(position, canonical key)` for every N-free k-mer window of `seq`.
///
/// Bases are 2-bit encoded (A=0, C=1, G=2, T=3, case-insensitive); a window
/// containing N or any other non-ACGT base is skipped (FastK splits on gaps
/// and its profile has a 0 at such positions). The rolling forward key and
/// its reverse complement (same rolling scheme as `pgi/build.rs`) select the
/// lexicographically smaller of the two strands.
pub fn canonical_keys(seq: &[u8], k: usize, mut emit: impl FnMut(usize, u128)) {
    let n = seq.len();
    if n < k {
        return;
    }
    let kmask = if 2 * k >= 128 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    };
    let rc_top = (2 * k - 2) as u32;
    let codes = base_codes();
    let mut kx: u128 = 0;
    let mut kxr: u128 = 0;
    let mut valid = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        let code = codes[b as usize];
        if code == 4 {
            kx = 0;
            kxr = 0;
            valid = 0;
        } else {
            kx = ((kx << 2) | code as u128) & kmask;
            kxr = (kxr >> 2) | (((3 - code) as u128) << rc_top);
            valid += 1;
        }
        if i + 1 >= k && valid >= k {
            emit(i + 1 - k, kx.min(kxr));
        }
    }
}

/// Base -> 2-bit code (0..3) or 4 (N / ambiguity); indexed by byte value.
pub(crate) fn base_codes() -> [u64; 256] {
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
