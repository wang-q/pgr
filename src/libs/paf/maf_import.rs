//! Conversion from MAF blocks (in `libs::fmt::maf`) into PAF records.
//!
//! Lives in `libs::paf` (rather than `libs::fmt::maf`) to keep the
//! `fmt -> paf` dependency edge one-directional: `paf` reads `fmt::maf`
//! structures, but `fmt::maf` never imports `paf`.

use crate::libs::fmt::maf::MafAli;
use crate::libs::paf::cigar::{
    block_identity, cigar_from_alignment, cigar_stats, cs_from_alignment, format_cigar,
    gap_compressed_identity,
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
    let cs = cs_from_alignment(ref_entry.text.as_bytes(), qry_entry.text.as_bytes())?;

    let mut tags = vec![
        format!("gi:f:{gi:.6}"),
        format!("bi:f:{bi:.6}"),
        format!("cg:Z:{cigar_str}"),
        format!("cs:Z:{cs}"),
    ];
    if let Some(s) = block.score {
        // MAF scores are f64 and may be negative; PAF `ms:i` is a signed
        // integer tag. Casting to u64 would silently clamp negative scores
        // to 0, so we round to the nearest integer and preserve the sign.
        tags.push(format!("ms:i:{}", s.round() as i64));
    }

    // MAF `start` for '-' strand is relative to the reverse-complemented
    // source; PAF query coordinates must be forward-strand 0-based half-open.
    let (q_start, q_end) = if qry_entry.strand == '-' {
        // Guard against malformed MAF where `start + size` exceeds `src_size`
        // (would underflow the usize subtraction in `reverse_range_pair`).
        let q_end_rc = qry_entry.start.checked_add(qry_entry.size).ok_or_else(|| {
            anyhow::anyhow!(
                "MAF query {} start+size overflow on reverse strand",
                qry_entry.src
            )
        })?;
        if q_end_rc > qry_entry.src_size {
            anyhow::bail!(
                "MAF query {} reverse-strand interval [{},{}) exceeds src_size {}",
                qry_entry.src,
                qry_entry.start,
                q_end_rc,
                qry_entry.src_size
            );
        }
        crate::libs::alignment::coords::reverse_range_pair(
            qry_entry.start,
            q_end_rc,
            qry_entry.src_size,
        )
    } else {
        (qry_entry.start, qry_entry.start + qry_entry.size)
    };

    let rec = PafRecord {
        query_name: qry_entry.src.clone(),
        query_length: qry_entry.src_size as u32,
        query_start: q_start as u32,
        query_end: q_end as u32,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::fmt::maf::MafComp;

    fn comp(
        src: &str,
        start: usize,
        size: usize,
        strand: char,
        src_size: usize,
        text: &str,
    ) -> MafComp {
        MafComp {
            src: src.to_string(),
            start,
            size,
            strand,
            src_size,
            text: text.to_string(),
        }
    }

    #[test]
    fn test_minus_strand_query_coords_converted_to_forward() {
        // query on '-' strand: start=5, size=12, src_size=100.
        // MAF reverse-strand interval [5, 17); forward = [100-17, 100-5) = [83, 95).
        let block = MafAli {
            score: None,
            components: vec![
                comp("ref", 0, 12, '+', 12, "ACGTACGTACGT"),
                comp("qry", 5, 12, '-', 100, "ACGTACGTACGT"),
            ],
        };
        let rec = maf_block_to_paf(&block)
            .unwrap()
            .expect("expected Some(rec)");
        assert_eq!(rec.strand, '-');
        assert_eq!(rec.query_start, 83);
        assert_eq!(rec.query_end, 95);
        assert_eq!(rec.query_length, 100);
        assert!(
            rec.tags.iter().any(|t| t.starts_with("cs:Z:")),
            "MAF-to-PAF must emit a cs:Z tag"
        );
    }

    #[test]
    fn test_plus_strand_query_coords_unchanged() {
        let block = MafAli {
            score: None,
            components: vec![
                comp("ref", 10, 12, '+', 50, "ACGTACGTACGT"),
                comp("qry", 20, 12, '+', 100, "ACGTACGTACGT"),
            ],
        };
        let rec = maf_block_to_paf(&block)
            .unwrap()
            .expect("expected Some(rec)");
        assert_eq!(rec.strand, '+');
        assert_eq!(rec.query_start, 20);
        assert_eq!(rec.query_end, 32);
    }

    #[test]
    fn test_minus_strand_interval_exceeding_src_size_rejected() {
        // Malformed MAF: reverse-strand query interval [5, 200) exceeds
        // src_size 100. Without the guard this would underflow the usize
        // subtraction in `reverse_range_pair` (debug panic). It must instead
        // return a friendly error.
        let block = MafAli {
            score: None,
            components: vec![
                comp("ref", 0, 12, '+', 12, "ACGTACGTACGT"),
                comp("qry", 5, 195, '-', 100, "ACGTACGTACGT"),
            ],
        };
        let err = maf_block_to_paf(&block).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("exceeds src_size"),
            "expected a clear error, got: {msg}"
        );
    }

    #[test]
    fn test_negative_score_preserved_in_ms_tag() {
        // MAF alignment scores may be negative. The `ms:i` tag must keep the
        // sign (→ i64), not clamp to 0 via a u64 cast.
        let block = MafAli {
            score: Some(-42.6),
            components: vec![
                comp("ref", 0, 12, '+', 12, "ACGTACGTACGT"),
                comp("qry", 0, 12, '+', 12, "ACGTACGTACGT"),
            ],
        };
        let rec = maf_block_to_paf(&block)
            .unwrap()
            .expect("expected Some(rec)");
        assert!(
            rec.tags.iter().any(|t| t == "ms:i:-43"),
            "negative score must round to i64 and keep sign, tags={:?}",
            rec.tags
        );
    }
}
