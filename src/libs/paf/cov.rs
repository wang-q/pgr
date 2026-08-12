//! Per-target coverage depth from PAF `cg:Z` tags.
//!
//! M/=/X/D operations cover the target (deletions still mean the target
//! region is aligned); insertions do not. Records without a `cg:Z` tag
//! contribute nothing (e.g. `psl to-paf` output, which skips the tag).

use super::cigar::extract_cigar;
use super::record::PafRecord;
use anyhow::Result;
use std::collections::BTreeMap;

/// Per-target coverage segments: (start, end, depth), 0-based half-open.
pub type CoverageSegments = BTreeMap<String, Vec<(u32, u32, u32)>>;

/// Coverage segments per target: (start, end, depth), 0-based half-open,
/// maximal runs of constant depth >= `min_depth`, sorted by target name then
/// start. Targets with no covering record are absent.
pub fn coverage_segments(records: &[PafRecord], min_depth: u32) -> Result<CoverageSegments> {
    let mut events: BTreeMap<String, Vec<(u32, i64)>> = BTreeMap::new();
    let mut target_len: BTreeMap<String, u32> = BTreeMap::new();
    for rec in records {
        target_len
            .entry(rec.target_name.clone())
            .and_modify(|l| *l = (*l).max(rec.target_length))
            .or_insert(rec.target_length);
        let cigar = extract_cigar(&rec.tags)?;
        if cigar.is_empty() {
            continue;
        }
        let mut tpos = rec.target_start;
        for op in cigar {
            if op.op() != 'I' {
                let len = op.len();
                events
                    .entry(rec.target_name.clone())
                    .or_default()
                    .push((tpos, 1));
                events
                    .entry(rec.target_name.clone())
                    .or_default()
                    .push((tpos + len, -1));
            }
            tpos += op.target_delta();
        }
    }

    let mut out = BTreeMap::new();
    for (target, mut ev) in events {
        ev.sort_unstable_by_key(|(p, _)| *p);
        let end = target_len.get(&target).copied().unwrap_or(0);
        let mut depth: i64 = 0;
        let mut last_pos = 0u32;
        let mut segs: Vec<(u32, u32, u32)> = Vec::new();
        for (pos, delta) in ev {
            if pos > last_pos && depth >= min_depth as i64 {
                let d = depth as u32;
                if let Some(last) = segs.last_mut() {
                    if last.2 == d {
                        last.1 = pos;
                    } else {
                        segs.push((last_pos, pos, d));
                    }
                } else {
                    segs.push((last_pos, pos, d));
                }
            }
            depth += delta;
            last_pos = pos;
        }
        if depth >= min_depth as i64 && last_pos < end {
            segs.push((last_pos, end, depth as u32));
        }
        out.insert(target, segs);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(target: &str, tlen: u32, tstart: u32, cigar: &str) -> PafRecord {
        PafRecord {
            query_name: "q".into(),
            query_length: 0,
            query_start: 0,
            query_end: 0,
            strand: '+',
            target_name: target.into(),
            target_length: tlen,
            target_start: tstart,
            target_end: 0,
            matches: 0,
            block_length: 0,
            mapq: 255,
            tags: vec![format!("cg:Z:{cigar}")],
        }
    }

    /// M/D cover the target, I does not; adjacent equal-depth segments merge.
    #[test]
    fn sweeps_constant_depth_segments() {
        let records = vec![
            rec("t", 100, 0, "10M"),   // covers 0..10
            rec("t", 100, 0, "5M"),    // covers 0..5
            rec("t", 100, 5, "5M"),    // covers 5..10
            rec("t", 100, 10, "5I5M"), // covers 10..15 (I skips target)
            rec("t", 100, 10, "5D5M"), // D covers 10..15, M covers 15..20
        ];
        let segs = coverage_segments(&records, 1).unwrap();
        let t = segs.get("t").unwrap();
        // 0..10 depth 2 (rec1+rec2/rec3, merged), 10..15 depth 2
        // (rec4 + rec5's D, merged with the left run), 15..20 depth 1.
        assert_eq!(*t, vec![(0, 15, 2), (15, 20, 1)]);
    }

    /// Records without `cg:Z` are ignored.
    #[test]
    fn ignores_records_without_cigar() {
        let records = vec![PafRecord {
            query_name: "q".into(),
            query_length: 0,
            query_start: 0,
            query_end: 0,
            strand: '+',
            target_name: "t".into(),
            target_length: 10,
            target_start: 0,
            target_end: 0,
            matches: 0,
            block_length: 0,
            mapq: 255,
            tags: vec![],
        }];
        let segs = coverage_segments(&records, 1).unwrap();
        assert!(segs.is_empty());
    }

    /// `min_depth` filters shallow segments.
    #[test]
    fn filters_by_min_depth() {
        let records = vec![rec("t", 100, 0, "10M"), rec("t", 100, 5, "10M")];
        let segs = coverage_segments(&records, 2).unwrap();
        assert_eq!(*segs.get("t").unwrap(), vec![(5, 10, 2)]);
    }
}
