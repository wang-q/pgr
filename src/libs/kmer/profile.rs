//! Per-sequence k-mer profiles (FastK `-p` / `-p:<table>` equivalents).

use super::KmerTable;

/// FastK caps per-k-mer counts at 32767 (`0x7fff`); profile values are `u16`.
const PROFILE_CAP: u16 = 0x7fff;

/// Profile value per k-mer position from the dataset-wide count (self).
///
/// `profiles[i][p]` is the count of the canonical k-mer at position `p` of
/// sequence `i`, or 0 when the window contains N (FastK splits on gaps).
pub fn self_profiles(seqs: &[Vec<u8>], k: usize, table: &KmerTable) -> Vec<Vec<u16>> {
    table_profiles(seqs, k, table)
}

/// Profile value per k-mer position from the repeat-table count (relative).
///
/// Values are the table count of the k-mer, or 0 when the window contains N
/// or the k-mer is absent from the table (FastK `-p:<table>` semantics).
pub fn relative_profiles(seqs: &[Vec<u8>], k: usize, table: &KmerTable) -> Vec<Vec<u16>> {
    table_profiles(seqs, k, table)
}

/// Shared scan: each valid window's canonical key is looked up in `table` and
/// its count (capped to the FastK 32767 limit) becomes the profile value;
/// N-containing windows and missing keys profile as 0.
fn table_profiles(seqs: &[Vec<u8>], k: usize, table: &KmerTable) -> Vec<Vec<u16>> {
    seqs.iter()
        .map(|seq| {
            let n = seq.len();
            if n < k {
                return Vec::new();
            }
            let mut out = vec![0u16; n - k + 1];
            super::canonical_keys(seq, k, |p, key| {
                let idx = table.keys.partition_point(|&x| x < key);
                let value = if idx < table.keys.len() && table.keys[idx] == key {
                    table.counts[idx].min(PROFILE_CAP as u32) as u16
                } else {
                    0
                };
                out[p] = value;
            });
            out
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::kmer::count::build_table;

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

    /// A random block whose k-mer windows are all canonical-unique.
    fn unique_block(k: usize, seed0: u64) -> Vec<u8> {
        (0..100u64)
            .map(|i| random_block(80, seed0 + i))
            .find(|b| {
                build_table(std::slice::from_ref(b), k).unwrap().keys.len() == b.len() - k + 1
            })
            .expect("a collision-free block must exist")
    }

    #[test]
    fn self_profile_counts_across_sequences() {
        // Two identical copies as separate sequences: every block k-mer
        // counts 2 in the combined table, so both profiles are all 2.
        let block = unique_block(6, 42);
        let seqs = vec![block.clone(), block.clone()];
        let table = build_table(&seqs, 6).unwrap();
        let profiles = self_profiles(&seqs, 6, &table);
        assert_eq!(profiles.len(), 2);
        assert!(profiles[0].iter().all(|&v| v == 2));
        assert!(profiles[1].iter().all(|&v| v == 2));
    }

    #[test]
    fn profile_zero_at_n_windows() {
        let seq = b"ACGTACGTNNACGTACGT".to_vec();
        let table = build_table(std::slice::from_ref(&seq), 4).unwrap();
        let profiles = self_profiles(std::slice::from_ref(&seq), 4, &table);
        // 15 windows; the 5 touching N (starts 5..9) must be 0, others > 0.
        assert_eq!(profiles[0].len(), 15);
        for (p, &v) in profiles[0].iter().enumerate() {
            assert_eq!(v == 0, (5..10).contains(&p), "window {p} value {v}");
        }
    }

    #[test]
    fn relative_profile_uses_table_count() {
        // Lib holds one copy; genome has two copies. The relative profile
        // must report the table count (1), not the genomic count (2), and 0
        // for windows overlapping the N gap (absent from the table).
        let lib = unique_block(6, 7);
        let k = 6;
        let table = build_table(std::slice::from_ref(&lib), k).unwrap();
        let mut genome = lib.clone();
        genome.extend(std::iter::repeat_n(b'N', k));
        genome.extend_from_slice(&lib);
        let profiles = relative_profiles(&[genome], k, &table);
        let first = lib.len() - k + 1; // windows fully inside the first copy
        let zeros = 2 * k - 1; // windows overlapping the N run
        assert!(profiles[0][..first].iter().all(|&v| v == 1));
        assert!(profiles[0][first..first + zeros].iter().all(|&v| v == 0));
        assert!(profiles[0][first + zeros..].iter().all(|&v| v == 1));

        // A sequence sharing no canonical k-mer with the table profiles as 0.
        let other = (0..100u64)
            .map(|i| random_block(22, 900 + i))
            .find(|b| {
                relative_profiles(std::slice::from_ref(b), k, &table)[0]
                    .iter()
                    .all(|&v| v == 0)
            })
            .expect("an absent block must exist");
        assert_eq!(
            relative_profiles(&[other], k, &table)[0]
                .iter()
                .sum::<u16>(),
            0
        );
    }
}
