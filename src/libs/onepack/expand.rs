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

use super::record::{AlnFile, AlnRecord, Skeleton};

use crate::libs::fmt::psl::Psl;
use crate::libs::paf::cigar::{cigar_from_alignment, format_cigar};
use crate::libs::paf::record::PafRecord;
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

/// A resolved skeleton contig: its scaffold name/length and the contig slice.
///
/// `seq` is the `[sbeg, sbeg+clen)` slice of the scaffold sequence. Holding
/// `sbeg`/`clen` lets callers map a contig-relative interval onto the forward
/// scaffold (needed to reverse `b` coordinates for reverse-complemented
/// records, matching FastGA `ALNtoPAF` `contigs[sbeg] + clen`).
pub struct ContigSeq<'a> {
    /// Scaffold (chromosome) name.
    pub name: &'a str,
    /// Scaffold length (in bases).
    pub slen: i64,
    /// Contig start offset within the scaffold.
    pub sbeg: i64,
    /// Contig length.
    pub clen: i64,
    /// The contig sequence slice.
    pub seq: &'a [u8],
}

/// Resolve a skeleton contig to its scaffold name/length and sequence slice.
///
/// `seqs` holds the scaffold sequences keyed by scaffold name; a contig is the
/// `[sbeg, sbeg+clen)` slice of its scaffold. Error if the scaffold (or the
/// contig slice) is not found.
pub fn contig_sequence<'a>(
    skeleton: &'a Skeleton,
    seqs: &'a [(String, Vec<u8>)],
    c: usize,
) -> Result<ContigSeq<'a>> {
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
    let sub = seq
        .get(sbeg..sbeg + clen)
        .ok_or_else(|| anyhow!("contig {c} slice out of scaffold bounds"))?;
    Ok(ContigSeq {
        name: scaff.name.as_str(),
        slen: scaff.slen,
        sbeg: contig.sbeg,
        clen: contig.clen,
        seq: sub,
    })
}

/// Loaded source genomes for the two sides of a `.1aln` file.
pub struct Genomes {
    /// Reference (`a`) side sequences, keyed by scaffold name.
    pub ref_seqs: Vec<(String, Vec<u8>)>,
    /// Query (`b`) side sequences, keyed by scaffold name.
    pub qry_seqs: Vec<(String, Vec<u8>)>,
}

/// Expand one record into a `PafRecord`, following the design doc §7.8 step 9
/// (query = `a` side, target = `b` side; reverse `b` coordinates when `comp`).
pub fn record_to_paf(
    rec: &AlnRecord,
    tspace: i64,
    genomes: &Genomes,
    aln: &AlnFile,
    with_cigar: bool,
) -> Result<PafRecord> {
    let a = contig_sequence(&aln.skeletons[0], &genomes.ref_seqs, rec.aread as usize)?;
    let b = contig_sequence(&aln.skeletons[1], &genomes.qry_seqs, rec.bread as usize)?;
    let (a_aln, b_aln) = expand_columns(rec, tspace, a.seq, b.seq)?;

    let blocksum = a_aln.len() as i64;
    let iid = blocksum - rec.diffs;
    let identity = if blocksum > 0 {
        iid as f64 / blocksum as f64
    } else {
        0.0
    };

    let mut tags = vec![format!("dv:f:{identity:.6}"), format!("df:i:{}", rec.diffs)];
    if with_cigar {
        // PAF CIGAR is query-vs-target. The PAF query is the `a` side and the
        // target the `b` side, so orient the CIGAR as `a`-vs-`b` by passing
        // `b_aln` as the reference argument (the second arg is the "query").
        let ops = cigar_from_alignment(&b_aln, &a_aln)?;
        tags.push(format!("cg:Z:{}", format_cigar(&ops)));
    }

    // Map the contig-relative `a`/`b` intervals onto the forward scaffolds. The
    // `.1aln` stores the `b` interval on the forward source; for a
    // reverse-complemented (`comp`) record the PAF target is reported in
    // forward orientation by reversing against the `b` contig end, matching
    // FastGA `ALNtoPAF` (`boff - bepos`/`boff - bbpos`).
    let (q_start, q_end) = (a.sbeg + rec.abpos, a.sbeg + rec.aepos);
    let (t_start, t_end) = if rec.comp {
        let boff = b.sbeg + b.clen;
        (boff - rec.bepos, boff - rec.bbpos)
    } else {
        (b.sbeg + rec.bbpos, b.sbeg + rec.bepos)
    };

    Ok(PafRecord {
        query_name: a.name.to_string(),
        query_length: a.slen as u32,
        query_start: q_start as u32,
        query_end: q_end as u32,
        strand: if rec.comp { '-' } else { '+' },
        target_name: b.name.to_string(),
        target_length: b.slen as u32,
        target_start: t_start as u32,
        target_end: t_end as u32,
        matches: iid as u32,
        block_length: blocksum as u32,
        mapq: 255,
        tags,
    })
}

/// Expand one record into a `Psl` (query = `a` side, target = `b` side).
pub fn record_to_psl(
    rec: &AlnRecord,
    tspace: i64,
    genomes: &Genomes,
    aln: &AlnFile,
) -> Result<Psl> {
    let a = contig_sequence(&aln.skeletons[0], &genomes.ref_seqs, rec.aread as usize)?;
    let b = contig_sequence(&aln.skeletons[1], &genomes.qry_seqs, rec.bread as usize)?;
    let (a_aln, b_aln) = expand_columns(rec, tspace, a.seq, b.seq)?;

    let a_str = String::from_utf8_lossy(&a_aln).into_owned();
    let b_str = String::from_utf8_lossy(&b_aln).into_owned();
    // strand: query `a` is always forward; target `b` is reverse when `comp`.
    let strand = if rec.comp { "+-" } else { "+" };
    // `Psl::from_align` expects forward-strand coordinates (it reverses them
    // internally for block construction), so pass the forward-mapped intervals.
    let (q_start, q_end) = (a.sbeg + rec.abpos, a.sbeg + rec.aepos);
    let (t_start, t_end) = if rec.comp {
        let boff = b.sbeg + b.clen;
        (boff - rec.bepos, boff - rec.bbpos)
    } else {
        (b.sbeg + rec.bbpos, b.sbeg + rec.bepos)
    };
    Psl::from_align(
        a.name,
        a.slen as u32,
        q_start as i32,
        q_end as i32,
        &a_str,
        b.name,
        b.slen as u32,
        t_start as i32,
        t_end as i32,
        &b_str,
        strand,
    )
    .ok_or_else(|| anyhow!("failed to build PSL for record {}:{}", rec.aread, rec.abpos))
}

/// Expand a record's trace into aligned columns, reverse-complementing `b`
/// when `comp` so the alignment is forward.
fn expand_columns(
    rec: &AlnRecord,
    tspace: i64,
    a_sub: &[u8],
    b_sub: &[u8],
) -> Result<(Vec<u8>, Vec<u8>)> {
    if rec.comp {
        let rc: Vec<u8> = crate::libs::nt::rev_comp(b_sub).collect();
        trace_to_columns(rec, tspace, a_sub, &rc)
    } else {
        trace_to_columns(rec, tspace, a_sub, b_sub)
    }
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
            let a = contig_sequence(&aln.skeletons[0], &ref_seqs, rec.aread as usize).unwrap();
            let b = contig_sequence(&aln.skeletons[1], &qry_seqs, rec.bread as usize).unwrap();
            assert!(!a.name.is_empty());
            assert!(!b.name.is_empty());
            let (a_aln, b_aln) = if rec.comp {
                let rc: Vec<u8> = crate::libs::nt::rev_comp(b.seq).collect();
                trace_to_columns(&rec, aln.tspace, a.seq, &rc).unwrap()
            } else {
                trace_to_columns(&rec, aln.tspace, a.seq, b.seq).unwrap()
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
