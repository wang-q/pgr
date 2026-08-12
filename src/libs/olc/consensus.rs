//! Layout stitching into consensus contigs (OLC stage S3).
//!
//! Overlaps are exact, so consensus is an exact stitch: walk each layout in
//! order, orient every unitig by its strand, append only the bases beyond
//! the exact overlap with the previous step, and verify the overlapping
//! prefix matches the already-stitched suffix.

use super::layout::Layout;
use super::overlap::Unitig;
use crate::libs::nt::rev_comp;
use anyhow::Result;

/// One consensus contig.
#[derive(Debug, Clone, PartialEq)]
pub struct Contig {
    /// Consensus sequence (5' -> 3' as laid out).
    pub seq: Vec<u8>,
    /// Approximate unitig depth (`sum(unitig lengths) / contig length`).
    pub coverage: f64,
}

/// Stitches every layout into a consensus contig.
///
/// Layouts shorter than `min_contig_len` are dropped. A layout whose
/// overlapping bases disagree with the already-stitched contig is an error
/// (exact overlaps must agree); the contig index is reported for debugging.
pub fn consensus(
    unitigs: &[Unitig],
    layouts: &[Layout],
    min_contig_len: usize,
) -> Result<Vec<Contig>> {
    let mut contigs = Vec::new();
    for (ci, layout) in layouts.iter().enumerate() {
        let mut seq: Vec<u8> = Vec::new();
        let mut total = 0usize;
        for (si, step) in layout.steps.iter().enumerate() {
            let mut piece: Vec<u8> = if step.strand == '+' {
                unitigs[step.unitig].seq.clone()
            } else {
                rev_comp(&unitigs[step.unitig].seq).collect()
            };
            total += piece.len();
            if si == 0 {
                seq.append(&mut piece);
                continue;
            }
            let overlap = step.overlap_len;
            anyhow::ensure!(
                overlap <= seq.len() && overlap <= piece.len(),
                "layout contig_{} step {si}: overlap {overlap} exceeds step lengths",
                ci + 1
            );
            let start = seq.len() - overlap;
            anyhow::ensure!(
                seq[start..] == piece[..overlap],
                "layout contig_{} step {si}: overlapping bases disagree \
                 (exact overlaps must match)",
                ci + 1
            );
            seq.extend_from_slice(&piece[overlap..]);
        }
        if seq.len() >= min_contig_len {
            let coverage = total as f64 / seq.len() as f64;
            contigs.push(Contig { seq, coverage });
        }
    }
    Ok(contigs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::nt::rev_comp;
    use crate::libs::olc::layout::LayoutStep;

    fn unitigs(names: &[&str], seqs: &[&str]) -> Vec<Unitig> {
        names
            .iter()
            .zip(seqs)
            .map(|(n, s)| Unitig {
                name: (*n).to_string(),
                seq: s.as_bytes().to_vec(),
            })
            .collect()
    }

    fn layout(steps: Vec<LayoutStep>) -> Layout {
        Layout { steps }
    }

    /// Forward chain stitches into the full contig.
    #[test]
    fn stitches_forward_chain() {
        let us = unitigs(
            &["u0", "u1", "u2"],
            &["AAAAAAAAACGTACGT", "ACGTACGTCCCCCCCC", "CCCCCCCCGGGGGGGG"],
        );
        let layouts = vec![layout(vec![
            LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: 16,
                overlap_len: 0,
            },
            LayoutStep {
                unitig: 1,
                strand: '+',
                q_start: 8,
                q_end: 24,
                overlap_len: 8,
            },
            LayoutStep {
                unitig: 2,
                strand: '+',
                q_start: 16,
                q_end: 32,
                overlap_len: 8,
            },
        ])];
        let contigs = consensus(&us, &layouts, 1).unwrap();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].seq, b"AAAAAAAAACGTACGTCCCCCCCCGGGGGGGG");
        // 48 unitig bases over a 32 bp contig.
        assert!((contigs[0].coverage - 1.5).abs() < 1e-9);
    }

    /// Reverse-strand step stitches via its reverse complement.
    #[test]
    fn stitches_reverse_step() {
        let u0 = "TTTTACGTAC";
        let u1 = String::from_utf8(rev_comp(b"ACGTACCCCC").collect()).unwrap();
        let us = unitigs(&["u0", "u1"], &[u0, &u1]);
        let layouts = vec![layout(vec![
            LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: 10,
                overlap_len: 0,
            },
            LayoutStep {
                unitig: 1,
                strand: '-',
                q_start: 4,
                q_end: 14,
                overlap_len: 6,
            },
        ])];
        let contigs = consensus(&us, &layouts, 1).unwrap();
        assert_eq!(contigs[0].seq, b"TTTTACGTACCCCC");
    }

    /// Short layouts are dropped by the minimum contig length.
    #[test]
    fn filters_short_contigs() {
        let us = unitigs(&["u0"], &["ACGTACGT"]);
        let layouts = vec![layout(vec![LayoutStep {
            unitig: 0,
            strand: '+',
            q_start: 0,
            q_end: 8,
            overlap_len: 0,
        }])];
        assert_eq!(consensus(&us, &layouts, 9).unwrap().len(), 0);
        assert_eq!(consensus(&us, &layouts, 8).unwrap().len(), 1);
    }

    /// A disagreeing overlap is a friendly error, not a panic.
    #[test]
    fn disagreeing_overlap_errors() {
        let us = unitigs(&["u0", "u1"], &["AAAACCCC", "CCCCGGGG"]);
        let layouts = vec![layout(vec![
            LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: 8,
                overlap_len: 0,
            },
            LayoutStep {
                unitig: 1,
                strand: '+',
                q_start: 2,
                q_end: 10,
                overlap_len: 6,
            },
        ])];
        // The claimed 6 bp overlap does not match (u0 suffix "AACCCC" vs
        // u1 prefix "CCCCGG"): the stitch must fail cleanly.
        let err = consensus(&us, &layouts, 1).unwrap_err();
        assert!(err.to_string().contains("disagree"), "{err}");
    }
}
