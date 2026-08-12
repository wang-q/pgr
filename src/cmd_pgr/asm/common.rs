//! Shared helpers for the asm command family.

use anyhow::Context;
use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use pgr::libs::olc::layout::Layout;
use pgr::libs::olc::overlap::{Overlap, OverlapType, Unitig};
use pgr::libs::paf::record::PafRecord;
use std::collections::HashSet;
use std::io::Write;
use std::path::Path;

/// Converts an overlap record into a PAF line.
pub fn to_paf(ov: &Overlap, unitigs: &[Unitig]) -> PafRecord {
    PafRecord {
        query_name: unitigs[ov.qid].name.clone(),
        query_length: unitigs[ov.qid].seq.len() as u32,
        query_start: ov.q_start as u32,
        query_end: ov.q_end as u32,
        strand: ov.strand,
        target_name: unitigs[ov.tid].name.clone(),
        target_length: unitigs[ov.tid].seq.len() as u32,
        target_start: ov.t_start as u32,
        target_end: ov.t_end as u32,
        matches: ov.length as u32,
        block_length: ov.length as u32,
        mapq: 255,
        tags: vec![format!(
            "ov:A:{}",
            match ov.otype {
                OverlapType::Dovetail => "D",
                OverlapType::Contain => "C",
            }
        )],
    }
}

/// Writes layouts as TSV lines (`contig_N<TAB>step<TAB>name<TAB>strand
/// <TAB>q_start<TAB>q_end<TAB>overlap_len`).
pub fn write_layout_tsv<W: Write>(
    out: &mut W,
    unitigs: &[Unitig],
    layouts: &[Layout],
) -> anyhow::Result<()> {
    for (ci, layout) in layouts.iter().enumerate() {
        for (si, step) in layout.steps.iter().enumerate() {
            writeln!(
                out,
                "contig_{}\t{}\t{}\t{}\t{}\t{}\t{}",
                ci + 1,
                si,
                unitigs[step.unitig].name,
                step.strand,
                step.q_start,
                step.q_end,
                step.overlap_len
            )?;
        }
    }
    Ok(())
}

/// `ByteBuilder.append(double, 1)`: half-up fixed-point with one decimal.
pub fn format_cov(x: f64) -> String {
    if x == x.trunc() {
        return format!("{}", x as i64);
    }
    let x = x + 0.05;
    let upper = x as i64;
    let lower = ((x - upper as f64) * 10.0) as i64;
    format!("{upper}.{lower}")
}

/// Reads all unitig FASTA files, prefixing names with a unique file tag.
///
/// The tag is the sanitized file stem (`stem:name`), so identical
/// `unitig_<id>` names across k files stay unique; collisions get a `.i`
/// suffix in file order (deterministic).
pub fn read_unitigs(infiles: &[String]) -> anyhow::Result<Vec<Unitig>> {
    let mut tags = Vec::with_capacity(infiles.len());
    let mut used = HashSet::new();
    for (i, path) in infiles.iter().enumerate() {
        let mut tag = tag_for(path);
        if used.contains(&tag) {
            tag = format!("{tag}.{i}");
        }
        used.insert(tag.clone());
        tags.push(tag);
    }
    let mut unitigs = Vec::new();
    for (path, tag) in infiles.iter().zip(tags) {
        let mut reader =
            SeqReader::new(path).with_context(|| format!("failed to open input {path}"))?;
        let mut rec = SeqRecord::new();
        while reader.read_record(&mut rec)? {
            unitigs.push(Unitig {
                name: format!("{tag}:{}", rec.name()),
                seq: rec.sequence().to_vec(),
            });
        }
    }
    Ok(unitigs)
}

/// File stem sanitized to `[A-Za-z0-9_.-]` (empty -> empty string).
fn tag_for(path: &str) -> String {
    let stem = Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    let mut tag = String::with_capacity(stem.len());
    for c in stem.chars() {
        if c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-') {
            tag.push(c);
        } else {
            tag.push('_');
        }
    }
    tag
}
