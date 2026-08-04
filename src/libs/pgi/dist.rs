//! Two-index merge distance (deterministic Jaccard/containment/Mash).

use super::{PgiIndex, PgiQuery};

/// Distance metrics between two indexes (by unique k-mer set).
#[derive(Debug, Clone, Copy)]
pub struct PgiDist {
    pub total1: u64,
    pub total2: u64,
    pub inter: u64,
    pub union: u64,
    pub mash: f32,
    pub jaccard: f32,
    pub containment: f32,
}

/// Validate that two indexes use identical sampling parameters.
pub fn validate_compatible(a: &impl PgiQuery, b: &impl PgiQuery) -> anyhow::Result<()> {
    anyhow::ensure!(
        a.k() == b.k(),
        "k-mer size mismatch: {} vs {}",
        a.k(),
        b.k()
    );
    anyhow::ensure!(
        a.smer() == b.smer(),
        "syncmer smer mismatch: {} vs {}",
        a.smer(),
        b.smer()
    );
    anyhow::ensure!(
        a.window() == b.window(),
        "syncmer window mismatch: {} vs {}",
        a.window(),
        b.window()
    );
    Ok(())
}

/// Merge two sorted k-mer key arrays; returns the intersection size.
fn merge_inter(a: &PgiIndex, b: &PgiIndex) -> u64 {
    let (mut i, mut j) = (0usize, 0usize);
    let mut inter = 0u64;
    while i < a.entries.len() && j < b.entries.len() {
        let ka = a.entries[i].kmer;
        let kb = b.entries[j].kmer;
        if ka == kb {
            inter += 1;
            i += 1;
            j += 1;
        } else if ka < kb {
            i += 1;
        } else {
            j += 1;
        }
    }
    inter
}

/// Compute distance metrics between two compatible indexes.
pub fn dist_between(a: &PgiIndex, b: &PgiIndex) -> anyhow::Result<PgiDist> {
    validate_compatible(a, b)?;
    let total1 = a.n_unique();
    let total2 = b.n_unique();
    let inter = merge_inter(a, b);
    let union = total1 + total2 - inter;
    let jaccard = if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    };
    let containment = if total1.min(total2) == 0 {
        0.0
    } else {
        inter as f64 / total1.min(total2) as f64
    };
    let mash = crate::libs::hash::mash_distance(jaccard, a.k) as f32;
    Ok(PgiDist {
        total1,
        total2,
        inter,
        union,
        mash,
        jaccard: jaccard as f32,
        containment: containment as f32,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::pgi::build::build_from_seqs;

    fn idx(seq: &[u8]) -> PgiIndex {
        build_from_seqs(
            vec![(String::from("c"), seq.to_vec())],
            10,
            4,
            2,
            true,
            false,
        )
        .unwrap()
    }

    #[test]
    fn identical_indexes_zero_distance() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let a = idx(&seq);
        let b = idx(&seq);
        let d = dist_between(&a, &b).unwrap();
        assert_eq!(d.jaccard, 1.0);
        assert!(d.mash < 1e-6);
    }

    #[test]
    fn disjoint_indexes_full_distance() {
        let a = idx(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let b = idx(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        let d = dist_between(&a, &b).unwrap();
        assert_eq!(d.inter, 0);
        assert_eq!(d.jaccard, 0.0);
    }

    #[test]
    fn parameter_mismatch_rejected() {
        let a = idx(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let b = build_from_seqs(
            vec![(
                String::from("c"),
                b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            )],
            12,
            4,
            2,
            true,
            false,
        )
        .unwrap();
        assert!(dist_between(&a, &b).is_err());
    }
}
