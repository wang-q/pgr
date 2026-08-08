//! Two-index merge distance (deterministic Jaccard/containment/Mash).

use super::PgiQuery;

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

/// Merge two sorted k-mer key streams in one pass, returning
/// `(unique_a, unique_b, intersection)`.
fn merge_stats(a: &impl PgiQuery, b: &impl PgiQuery) -> (u64, u64, u64) {
    let (a0, a1) = a.entry_range(0, u128::MAX);
    let (b0, b1) = b.entry_range(0, u128::MAX);
    let (mut i, mut j) = (a0, b0);
    let (mut total1, mut total2, mut inter) = (0u64, 0u64, 0u64);
    while i < a1 && j < b1 {
        let ka = a.entry_kmer(i);
        let kb = b.entry_kmer(j);
        if ka == kb {
            inter += 1;
            total1 += 1;
            total2 += 1;
            i = a.entry_next(i);
            j = b.entry_next(j);
        } else if ka < kb {
            total1 += 1;
            i = a.entry_next(i);
        } else {
            total2 += 1;
            j = b.entry_next(j);
        }
    }
    while i < a1 {
        total1 += 1;
        i = a.entry_next(i);
    }
    while j < b1 {
        total2 += 1;
        j = b.entry_next(j);
    }
    (total1, total2, inter)
}

/// Compute distance metrics between two compatible indexes (resident or
/// memory-mapped).
pub fn dist_between(a: &impl PgiQuery, b: &impl PgiQuery) -> anyhow::Result<PgiDist> {
    validate_compatible(a, b)?;
    let (total1, total2, inter) = merge_stats(a, b);
    let union = total1 + total2 - inter;
    // Empty indexes: two empty sets are identical (jaccard 1, distance 0);
    // containment is directional (first set as denominator), matching the
    // sketch-distance family.
    let jaccard = if union == 0 {
        1.0
    } else {
        inter as f64 / union as f64
    };
    let containment = if total1 == 0 {
        0.0
    } else {
        inter as f64 / total1 as f64
    };
    let mash = crate::libs::hash::mash_distance(jaccard, a.k()) as f32;
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
    use crate::libs::pgi::PgiIndex;

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

    #[test]
    fn dist_and_hv_mmap_equivalent_to_resident() {
        use crate::libs::pgi::to_hv::index_to_hv;
        use crate::libs::pgi::PgiMmap;

        let a = idx(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let b = idx(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let temp = tempfile::TempDir::new().unwrap();
        let pa = temp.path().join("a.pgi");
        let pb = temp.path().join("b.pgi");
        let mut fa = std::fs::File::create(&pa).unwrap();
        a.write(&mut fa).unwrap();
        drop(fa);
        let mut fb = std::fs::File::create(&pb).unwrap();
        b.write(&mut fb).unwrap();
        drop(fb);
        let ma = PgiMmap::open(&pa).unwrap();
        let mb = PgiMmap::open(&pb).unwrap();

        let dr = dist_between(&a, &b).unwrap();
        let dm = dist_between(&ma, &mb).unwrap();
        assert_eq!(dr.total1, dm.total1);
        assert_eq!(dr.total2, dm.total2);
        assert_eq!(dr.inter, dm.inter);
        assert_eq!(dr.union, dm.union);
        assert_eq!(dr.mash, dm.mash);

        assert_eq!(
            crate::libs::pgi::count_unique(&a),
            crate::libs::pgi::count_unique(&ma)
        );
        assert_eq!(index_to_hv(&a, 1024, 16), index_to_hv(&ma, 1024, 16));
    }
}
