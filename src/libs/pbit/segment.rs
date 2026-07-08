//! Segment-level delta compression/decompression.
//!
//! Thin wrapper around [`LzDiff`]. Handles the empty-delta special case
//! (text == reference) and 2-bit → ASCII conversion. Does NOT handle
//! flate2 compression, archive I/O, packing, or delta deduplication —
//! those are the responsibility of `Compressor`/`Decompressor`.

use anyhow::Result;

use super::lz_diff::{decode_base, LzDiff};

/// Segment-level delta compression/decompression.
pub struct Segment {
    lz_diff: LzDiff,
}

impl Segment {
    /// Create a new Segment with the given minimum match length.
    pub fn new(min_match_len: u32) -> Self {
        Self {
            lz_diff: LzDiff::new(min_match_len),
        }
    }

    /// Store the reference (2-bit encoded internally). Call before add/get.
    pub fn prepare(&mut self, ref_dna: &[u8]) {
        self.lz_diff.prepare(ref_dna);
    }

    /// Build the LZ-diff hash table. Call before add (not needed for get).
    pub fn prepare_index(&mut self) {
        self.lz_diff.prepare_index();
    }

    /// Encode `seq` against the prepared reference, return uncompressed delta.
    /// Empty delta means seq == reference (equal-sequences optimization).
    pub fn add(&mut self, seq: &[u8]) -> Vec<u8> {
        let mut delta = Vec::new();
        self.lz_diff.encode(seq, &mut delta);
        delta
    }

    /// Decode `delta` (uncompressed) back to the original ASCII sequence.
    /// Empty delta returns the reference as ASCII.
    pub fn get(&self, delta: &[u8]) -> Result<Vec<u8>> {
        let decoded_2bit: Vec<u8> = if delta.is_empty() {
            self.lz_diff.reference_2bit().to_vec()
        } else {
            let mut d = Vec::new();
            self.lz_diff.decode(delta, &mut d)?;
            d
        };
        Ok(decoded_2bit.iter().map(|&c| decode_base(c)).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_dna(len: usize, seed: u64) -> Vec<u8> {
        use rand::rngs::StdRng;
        use rand::Rng;
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(seed);
        (0..len)
            .map(|_| match rng.random_range(0u8..4) {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect()
    }

    fn random_dna_with_n(len: usize, seed: u64, n_freq: f64) -> Vec<u8> {
        use rand::rngs::StdRng;
        use rand::Rng;
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(seed);
        (0..len)
            .map(|_| {
                if rng.random_range(0.0..1.0) < n_freq {
                    b'N'
                } else {
                    match rng.random_range(0u8..4) {
                        0 => b'A',
                        1 => b'C',
                        2 => b'G',
                        _ => b'T',
                    }
                }
            })
            .collect()
    }

    #[test]
    fn test_roundtrip_single_seq() {
        let reference = random_dna(2000, 42);
        let text = random_dna(2000, 100);
        let mut seg = Segment::new(18);
        seg.prepare(&reference);
        seg.prepare_index();
        let delta = seg.add(&text);
        let decoded = seg.get(&delta).expect("decode failed");
        assert_eq!(decoded, text);
    }

    #[test]
    fn test_roundtrip_identical() {
        let reference = random_dna(2000, 42);
        let mut seg = Segment::new(18);
        seg.prepare(&reference);
        seg.prepare_index();
        let delta = seg.add(&reference);
        assert!(
            delta.is_empty(),
            "identical sequence should produce empty delta"
        );
        let decoded = seg.get(&delta).expect("decode failed");
        // 2-bit encoding is case-insensitive, so output is uppercase
        let expected: Vec<u8> = reference.iter().map(|&c| c.to_ascii_uppercase()).collect();
        assert_eq!(decoded, expected);
    }

    #[test]
    fn test_roundtrip_multiple_seqs() {
        let reference = random_dna(2000, 42);
        let seq1 = random_dna(2000, 100);
        let seq2 = random_dna(2000, 200);
        let mut seg = Segment::new(18);
        seg.prepare(&reference);
        seg.prepare_index();
        let delta1 = seg.add(&seq1);
        let delta2 = seg.add(&seq2);
        assert_eq!(seg.get(&delta1).unwrap(), seq1);
        assert_eq!(seg.get(&delta2).unwrap(), seq2);
    }

    #[test]
    fn test_roundtrip_with_n() {
        let reference = random_dna_with_n(2000, 42, 0.05);
        let text = random_dna_with_n(2000, 100, 0.05);
        let mut seg = Segment::new(18);
        seg.prepare(&reference);
        seg.prepare_index();
        let delta = seg.add(&text);
        let decoded = seg.get(&delta).expect("decode failed");
        // N and IUPAC codes both map to 4, decode_base maps 4 to 'N'
        let expected: Vec<u8> = text
            .iter()
            .map(|&c| {
                if c == b'N' || c == b'n' {
                    b'N'
                } else {
                    c.to_ascii_uppercase()
                }
            })
            .collect();
        assert_eq!(decoded, expected);
    }

    #[test]
    fn test_delta_size_similarity() {
        let reference = random_dna(2000, 42);
        // Similar: only 2 SNPs
        let mut similar = reference.clone();
        similar[100] = b'G';
        similar[500] = b'C';
        // Different: entirely random
        let different = random_dna(2000, 999);

        let mut seg = Segment::new(18);
        seg.prepare(&reference);
        seg.prepare_index();
        let delta_similar = seg.add(&similar);
        let delta_different = seg.add(&different);
        assert!(
            delta_similar.len() < delta_different.len(),
            "similar sequence should have smaller delta ({} < {})",
            delta_similar.len(),
            delta_different.len()
        );
    }

    #[test]
    fn test_get_without_index() {
        let reference = random_dna(2000, 42);
        let text = random_dna(2000, 100);
        let mut seg = Segment::new(18);
        seg.prepare(&reference);
        seg.prepare_index();
        let delta = seg.add(&text);

        // New Segment with same reference but no index
        let mut seg2 = Segment::new(18);
        seg2.prepare(&reference);
        // Do NOT call prepare_index
        let decoded = seg2.get(&delta).expect("decode failed");
        assert_eq!(decoded, text);
    }

    #[test]
    fn test_lowercase_input() {
        let mut reference = random_dna(2000, 42);
        for b in reference.iter_mut() {
            *b = b.to_ascii_lowercase();
        }
        let text = random_dna(2000, 100);
        let mut seg = Segment::new(18);
        seg.prepare(&reference);
        seg.prepare_index();
        let delta = seg.add(&text);
        let decoded = seg.get(&delta).expect("decode failed");
        // Output is always uppercase (2-bit encoding loses case info)
        let expected: Vec<u8> = text.iter().map(|&c| c.to_ascii_uppercase()).collect();
        assert_eq!(decoded, expected);
    }

    #[test]
    fn test_empty_reference_prepare() {
        let mut seg = Segment::new(18);
        seg.prepare(b"");
        // Should not panic; reference_2bit returns empty slice
        assert!(seg.lz_diff.reference_2bit().is_empty());
    }
}
