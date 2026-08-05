//! Shared read-side orchestration for `pgr 1aln to-paf` / `to-psl`.
//!
//! Loads the source genomes referenced by a `.1aln` file, expands each record's
//! trace points back into base-level aligned columns (via `onepack::expand`),
//! and exposes the resulting PAF / PSL records. Orchestration lives here (the
//! thin CLI layer) rather than in `libs::onepack`, per the design doc
//! `notes/design/1aln.md` §7.8.

use anyhow::{anyhow, Context, Result};
use std::io::Write;

use pgr::libs::fmt::psl::Psl;
use pgr::libs::onepack::expand::{contig_sequence, trace_to_columns};
use pgr::libs::onepack::record::{AlnFile, AlnRecord};
use pgr::libs::paf::cigar::cigar_from_alignment;
use pgr::libs::paf::record::PafRecord;
use pgr::libs::pgi::build::read_fasta;

/// Loaded source genomes for the two sides of a `.1aln` file.
pub struct Genomes {
    /// Reference (`a`) side sequences, keyed by scaffold name.
    pub ref_seqs: Vec<(String, Vec<u8>)>,
    /// Query (`b`) side sequences, keyed by scaffold name.
    pub qry_seqs: Vec<(String, Vec<u8>)>,
}

/// Open the `.1aln` file and load the two source genomes.
pub fn open_aln(args: &clap::ArgMatches) -> Result<(AlnFile, Genomes)> {
    let infile = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let ref_seq = args
        .get_one::<String>("ref_seq")
        .context("missing required argument: --ref-seq")?;
    let qry_seq = args
        .get_one::<String>("query_seq")
        .context("missing required argument: --query-seq")?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        [infile.as_str(), ref_seq.as_str(), qry_seq.as_str()],
    )?;

    let aln =
        AlnFile::open(infile).with_context(|| format!("Failed to open .1aln file {infile}"))?;
    let genomes = Genomes {
        ref_seqs: read_fasta(ref_seq)
            .with_context(|| format!("Failed to load reference FASTA {ref_seq}"))?,
        qry_seqs: read_fasta(qry_seq)
            .with_context(|| format!("Failed to load query FASTA {qry_seq}"))?,
    };
    Ok((aln, genomes))
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
    let (a_name, a_slen, a_sub) =
        contig_sequence(&aln.skeletons[0], &genomes.ref_seqs, rec.aread as usize)?;
    let (b_name, b_slen, b_sub) =
        contig_sequence(&aln.skeletons[1], &genomes.qry_seqs, rec.bread as usize)?;
    let (a_aln, b_aln) = expand_columns(rec, tspace, a_sub, b_sub)?;

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
        tags.push(format!(
            "cg:Z:{}",
            pgr::libs::paf::cigar::format_cigar(&ops)
        ));
    }

    // Forward `b` coordinates. The `.1aln` stores the `b` interval on the
    // forward source ([bbpos, bepos], bepos > bbpos); when `comp` the PAF
    // target is still reported in forward orientation (matching FastGA
    // ALNtoPAF), with the `-` strand flag carrying the orientation.
    let (b_start, b_end) = (rec.bbpos, rec.bepos);

    Ok(PafRecord {
        query_name: a_name.to_string(),
        query_length: a_slen as u32,
        query_start: rec.abpos as u32,
        query_end: rec.aepos as u32,
        strand: if rec.comp { '-' } else { '+' },
        target_name: b_name.to_string(),
        target_length: b_slen as u32,
        target_start: b_start as u32,
        target_end: b_end as u32,
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
    let (a_name, a_slen, a_sub) =
        contig_sequence(&aln.skeletons[0], &genomes.ref_seqs, rec.aread as usize)?;
    let (b_name, b_slen, b_sub) =
        contig_sequence(&aln.skeletons[1], &genomes.qry_seqs, rec.bread as usize)?;
    let (a_aln, b_aln) = expand_columns(rec, tspace, a_sub, b_sub)?;

    let a_str = String::from_utf8_lossy(&a_aln).into_owned();
    let b_str = String::from_utf8_lossy(&b_aln).into_owned();
    // strand: query `a` is always forward; target `b` is reverse when `comp`.
    let strand = if rec.comp { "+-" } else { "+" };
    Psl::from_align(
        a_name,
        a_slen as u32,
        rec.abpos as i32,
        rec.aepos as i32,
        &a_str,
        b_name,
        b_slen as u32,
        rec.bbpos as i32,
        rec.bepos as i32,
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
        let rc: Vec<u8> = pgr::libs::nt::rev_comp(b_sub).collect();
        trace_to_columns(rec, tspace, a_sub, &rc)
    } else {
        trace_to_columns(rec, tspace, a_sub, b_sub)
    }
}

/// Write a single PAF record to `out`.
pub fn write_paf<W: Write>(out: &mut W, rec: &PafRecord) -> Result<()> {
    pgr::libs::paf::record::write_paf_record(out, rec)?;
    Ok(())
}
