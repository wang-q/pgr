//! Read-side trace expansion: turn decoded `.1aln` records back into
//! base-level aligned columns (and from there PAF/PSL).
//!
//! A `.1aln` record stores the alignment as trace points: the `b`-coordinate
//! of the path at every `tspace`-multiple of the `a`-axis, plus the number of
//! differences in each interval. To recover the base-level alignment we fill
//! the path between consecutive trace points with an in-box DP. pgr reuses
//! `pgi::wave::banded_edit_ops` for that box DP, so this module is a thin
//! orchestration layer (`P5` in the design doc `notes/design/1aln.md`).

use anyhow::{anyhow, bail, Result};

use super::record::{AlnRecord, Skeleton};

use crate::libs::pgi::wave::{banded_edit_ops, ops_to_columns, EditOp};

/// Expand one record's trace into aligned columns `(a_aln, b_aln)`.
///
/// `a_seq` is the full reference (`a`) contig; `b_seq` is the full query (`b`)
/// contig already placed in aligned orientation (reverse-complemented when
/// `comp`, so the alignment is forward). The returned columns are in the
/// forward orientation of the trace (a increasing, b increasing in `b_seq`).
pub fn trace_to_columns(
    rec: &AlnRecord,
    tspace: i64,
    a_seq: &[u8],
    b_seq: &[u8],
) -> Result<(Vec<u8>, Vec<u8>)> {
    let k = rec.num_points();
    if k == 0 {
        bail!("empty trace for record {rec:?}");
    }
    // a-boundaries: trace points sit at tspace multiples of `a`; the last
    // segment ends at aepos. b-boundaries come from the running tpoints.
    let mut a_pts = vec![0i64; k + 1];
    a_pts[0] = rec.abpos;
    let mut next = (rec.abpos / tspace) * tspace + tspace;
    for slot in a_pts.iter_mut().take(k).skip(1) {
        *slot = next;
        next += tspace;
    }
    a_pts[k] = rec.aepos;

    let mut b_pts = vec![0i64; k + 1];
    b_pts[0] = rec.bbpos;
    for i in 1..k {
        b_pts[i] = b_pts[i - 1] + rec.tpoints[i - 1];
    }
    b_pts[k] = rec.bepos;

    let mut ops: Vec<EditOp> = Vec::new();
    for i in 0..k {
        let (as_, ae) = (a_pts[i] as usize, a_pts[i + 1] as usize);
        let (bs, be) = (b_pts[i] as usize, b_pts[i + 1] as usize);
        let q = &a_seq[as_..ae];
        let t = &b_seq[bs..be];
        // Diagonal k = t_pos - q_pos. d0 on the segment start; cover the whole
        // box so banded_edit_ops degenerates to the exact min-edit DP.
        let d0 = bs as i64 - as_ as i64;
        banded_edit_ops(
            q,
            t,
            as_,
            bs,
            d0 - ae as i64 + as_ as i64,
            d0 + be as i64 - bs as i64,
            &mut ops,
        );
    }
    let (a_aln, b_aln, _) = ops_to_columns(
        a_seq,
        b_seq,
        rec.abpos as usize,
        rec.bbpos as usize,
        rec.aepos as usize,
        rec.bepos as usize,
        &ops,
    );
    Ok((a_aln, b_aln))
}

/// Resolve a skeleton contig to `(scaffold_name, scaffold_length, contig_seq)`.
///
/// `seqs` holds the scaffold sequences keyed by scaffold name; a contig is the
/// `[sbeg, sbeg+clen)` slice of its scaffold. Error if the scaffold (or the
/// contig slice) is not found.
pub fn contig_sequence<'a>(
    skeleton: &'a Skeleton,
    seqs: &'a [(String, Vec<u8>)],
    c: usize,
) -> Result<(&'a str, i64, &'a [u8])> {
    let contig = skeleton
        .contigs
        .get(c)
        .ok_or_else(|| anyhow!("contig index {c} out of range"))?;
    let scaff = skeleton
        .scaffolds
        .get(contig.scaf)
        .ok_or_else(|| anyhow!("scaffold index {} out of range", contig.scaf))?;
    let (_, seq) = seqs
        .iter()
        .find(|(name, _)| name == &scaff.name)
        .ok_or_else(|| anyhow!("scaffold sequence '{}' not found", scaff.name))?;
    let sbeg = contig.sbeg as usize;
    let clen = contig.clen as usize;
    let seq: &[u8] = seq;
    let sub = seq
        .get(sbeg..sbeg + clen)
        .ok_or_else(|| anyhow!("contig {c} slice out of scaffold bounds"))?;
    Ok((scaff.name.as_str(), scaff.slen, sub))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::onepack::record::AlnFile;
    use crate::libs::paf::cigar::cigar_from_alignment;
    use crate::libs::pgi::build::read_fasta;

    fn golden_path() -> String {
        env!("CARGO_MANIFEST_DIR").to_string() + "/tests/genome/mg1655-sakai.1aln"
    }

    #[allow(clippy::type_complexity)] // test helper returning the two FASTA genomes
    fn load_genomes() -> (Vec<(String, Vec<u8>)>, Vec<(String, Vec<u8>)>) {
        let root = env!("CARGO_MANIFEST_DIR").to_string() + "/tests/genome/";
        let ref_seqs = read_fasta(&(root.clone() + "mg1655.fa.gz")).unwrap();
        let qry_seqs = read_fasta(&(root + "sakai.fa.gz")).unwrap();
        (ref_seqs, qry_seqs)
    }

    #[test]
    fn expands_first_records() {
        let (ref_seqs, qry_seqs) = load_genomes();
        let mut aln = AlnFile::open(&golden_path()).unwrap();
        // The first skeleton's scaffold name must match the ref genome.
        let g0 = &aln.skeletons[0];
        let scaff0_name = g0.scaffolds[0].name.clone();
        assert!(ref_seqs.iter().any(|(n, _)| *n == scaff0_name));
        let mut count = 0;
        while let Some(rec) = aln.next_record().unwrap() {
            let (a_name, _, a_sub) =
                contig_sequence(&aln.skeletons[0], &ref_seqs, rec.aread as usize).unwrap();
            let (b_name, _, b_sub) =
                contig_sequence(&aln.skeletons[1], &qry_seqs, rec.bread as usize).unwrap();
            assert!(!a_name.is_empty());
            assert!(!b_name.is_empty());
            let (a_aln, b_aln) = if rec.comp {
                let rc: Vec<u8> = crate::libs::nt::rev_comp(b_sub).collect();
                trace_to_columns(&rec, aln.tspace, a_sub, &rc).unwrap()
            } else {
                trace_to_columns(&rec, aln.tspace, a_sub, b_sub).unwrap()
            };
            assert_eq!(a_aln.len(), b_aln.len());
            assert!(!a_aln.is_empty());
            // The alignment interval in the columns matches the record bounds
            // (count non-gap bases).
            let a_bases = a_aln.iter().filter(|&&b| b != b'-').count();
            let b_bases = b_aln.iter().filter(|&&b| b != b'-').count();
            assert_eq!(a_bases as i64, rec.aepos - rec.abpos);
            assert_eq!(b_bases as i64, rec.bepos - rec.bbpos);
            // CIGAR is well-formed over these columns.
            let ops = cigar_from_alignment(&a_aln, &b_aln).unwrap();
            assert!(!ops.is_empty());
            count += 1;
            if count >= 5 {
                break;
            }
        }
        assert_eq!(count, 5);
    }
}
