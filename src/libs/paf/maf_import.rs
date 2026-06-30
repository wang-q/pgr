//! Conversion from MAF blocks (in `libs::fmt::maf`) into PAF records.
//!
//! Lives in `libs::paf` (rather than `libs::fmt::maf`) to keep the
//! `fmt -> paf` dependency edge one-directional: `paf` reads `fmt::maf`
//! structures, but `fmt::maf` never imports `paf`.

use crate::libs::fmt::maf::MafAli;
use crate::libs::paf::cigar::{
    block_identity, cigar_from_alignment, cigar_stats, format_cigar, gap_compressed_identity,
};
use crate::libs::paf::record::PafRecord;

/// Convert a two-sequence MAF block into a PAF record.
///
/// Returns `Ok(None)` if the block does not have exactly two components.
/// Multi-sequence blocks are skipped (caller should log a warning if desired).
pub fn maf_block_to_paf(block: &MafAli) -> anyhow::Result<Option<PafRecord>> {
    if block.components.len() < 2 {
        return Ok(None);
    }
    if block.components.len() > 2 {
        return Ok(None); // caller logs warning
    }

    let ref_entry = &block.components[0];
    let qry_entry = &block.components[1];

    let cigar_ops = cigar_from_alignment(ref_entry.text.as_bytes(), qry_entry.text.as_bytes())?;
    let stats = cigar_stats(&cigar_ops);
    let gi = gap_compressed_identity(&cigar_ops);
    let bi = block_identity(&cigar_ops);
    let cigar_str = format_cigar(&cigar_ops);

    let mut tags = vec![
        format!("gi:f:{gi:.6}"),
        format!("bi:f:{bi:.6}"),
        format!("cg:Z:{cigar_str}"),
    ];
    if let Some(s) = block.score {
        tags.push(format!("ms:i:{}", s as u64));
    }

    let rec = PafRecord {
        query_name: qry_entry.src.clone(),
        query_length: qry_entry.src_size as u32,
        query_start: qry_entry.start as u32,
        query_end: (qry_entry.start + qry_entry.size) as u32,
        strand: qry_entry.strand,
        target_name: ref_entry.src.clone(),
        target_length: ref_entry.src_size as u32,
        target_start: ref_entry.start as u32,
        target_end: (ref_entry.start + ref_entry.size) as u32,
        matches: stats.matches,
        block_length: crate::libs::paf::cigar::block_length(&stats),
        mapq: 255,
        tags,
    };

    Ok(Some(rec))
}
