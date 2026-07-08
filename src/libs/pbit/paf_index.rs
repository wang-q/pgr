//! pbit-local PAF query-side index for CIGAR-driven delta encoding.
//!
//! Builds an independent `BasicCOITree<PafAlign, u32>` keyed by `query_id`,
//! interval = query coordinates. Stores the ORIGINAL PAF metadata (no role
//! swap, no CIGAR reversal), unlike `PafIndex.reverse_trees` which only holds
//! `+` strand mirror entries with swapped query/target roles.
//!
//! Used by `Compressor::append_sample_with_paf` to locate alignments covering
//! each sample segment.

use std::collections::HashMap;
use std::io::BufRead;

use anyhow::{Context, Result};
use coitrees::{BasicCOITree, Interval, IntervalTree};
use indexmap::IndexMap;

use crate::libs::paf::cigar::{extract_cigar, CigarOp};
use crate::libs::paf::parser::parse_paf_line;

/// One PAF alignment stored in the query-side index. Holds the ORIGINAL
/// (non-swapped) fields needed by pbit's CIGAR encoding path. `cigar` is
/// already `extract_cigar`-parsed; empty Vec means the record was skipped
/// during indexing.
#[derive(Debug, Clone)]
pub struct PafAlign {
    pub query_id: u32,
    pub query_start: i32,
    pub query_end: i32,
    pub target_name: String,
    pub target_start: i32,
    pub target_end: i32,
    pub strand: char,
    pub cigar: Vec<CigarOp>,
}

/// Query-side PAF index: name→id map + per-query interval trees.
pub struct PafQueryIndex {
    pub names: IndexMap<String, u32>,
    pub trees: HashMap<u32, BasicCOITree<PafAlign, u32>>,
}

impl PafQueryIndex {
    /// Build from a PAF reader. Lines are parsed one-by-one: a single
    /// malformed record is skipped with `log::warn` (design §11 decision 8),
    /// and records without `cg:Z:` CIGAR are skipped silently (decision 7).
    /// If every non-empty line fails to parse (non-PAF format), bail.
    /// An empty file yields an empty index (caller falls back to LZ-diff).
    pub fn build<R: BufRead>(reader: R) -> Result<Self> {
        let mut names: IndexMap<String, u32> = IndexMap::new();
        let mut by_query: HashMap<u32, Vec<Interval<PafAlign>>> = HashMap::new();
        let mut total_lines = 0usize;
        let mut failed_count = 0usize;

        for line in reader.lines() {
            let line = line.context("failed to read PAF line")?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            total_lines += 1;
            let rec = match parse_paf_line(&line) {
                Ok(r) => r,
                Err(e) => {
                    log::warn!("skipping invalid PAF line: {}: {}", line, e);
                    failed_count += 1;
                    continue;
                }
            };
            let cigar = match extract_cigar(&rec.tags) {
                Ok(c) => c,
                Err(e) => {
                    log::warn!(
                        "skipping PAF record with malformed CIGAR tag: {}: {}",
                        line,
                        e
                    );
                    failed_count += 1;
                    continue;
                }
            };
            if cigar.is_empty() {
                // No CIGAR tag: skip this record (design §11 decision 7). Not a failure.
                continue;
            }

            // Allocate query_id by query_name (insertion order).
            let next_id = names.len() as u32;
            let query_id = *names.entry(rec.query_name.clone()).or_insert(next_id);

            let meta = PafAlign {
                query_id,
                query_start: rec.query_start as i32,
                query_end: rec.query_end as i32,
                target_name: rec.target_name.clone(),
                target_start: rec.target_start as i32,
                target_end: rec.target_end as i32,
                strand: rec.strand,
                cigar,
            };
            by_query.entry(query_id).or_default().push(Interval::new(
                rec.query_start as i32,
                rec.query_end as i32,
                meta,
            ));
        }

        // All non-empty lines failed to parse → treat as non-PAF format (decision 8).
        if total_lines > 0 && failed_count == total_lines {
            anyhow::bail!(
                "PAF file contains no valid records (all {} lines failed to parse)",
                total_lines
            );
        }

        // Build sorted interval trees (one per query_id).
        let mut trees: HashMap<u32, BasicCOITree<PafAlign, u32>> = HashMap::new();
        for (qid, mut intervals) in by_query {
            intervals.sort_by_key(|iv| iv.first);
            trees.insert(qid, BasicCOITree::new(&intervals));
        }

        Ok(Self { names, trees })
    }

    /// Build from a PAF file path (supports `stdin` and `.gz`).
    pub fn build_from_path(path: &str) -> Result<Self> {
        let reader = crate::libs::io::reader(path)
            .with_context(|| format!("failed to open PAF file: {}", path))?;
        Self::build(reader)
    }

    /// Look up query_id by query name.
    pub fn query_id(&self, name: &str) -> Option<u32> {
        self.names.get(name).copied()
    }

    /// Return true if the index has no alignments (all records were skipped
    /// or the PAF file was empty). Used by the compressor to decide whether
    /// to fall back to pure LZ-diff.
    pub fn is_empty(&self) -> bool {
        self.trees.is_empty()
    }

    /// Query all alignments overlapping `[start, end)` on the given query.
    /// Returns owned clones of the matching `PafAlign` entries (the
    /// `coitrees::IntervalTree::query` closure borrows nodes with a lifetime
    /// scoped to the closure body, so references cannot escape).
    pub fn query(&self, query_id: u32, start: i32, end: i32) -> Vec<PafAlign> {
        let mut out: Vec<PafAlign> = Vec::new();
        if let Some(tree) = self.trees.get(&query_id) {
            tree.query(start, end, |iv| {
                out.push(iv.metadata.clone());
            });
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[allow(clippy::too_many_arguments)]
    fn paf_line(
        qname: &str,
        qs: u32,
        qe: u32,
        strand: &str,
        tname: &str,
        ts: u32,
        te: u32,
        cigar: &str,
    ) -> String {
        // query_len=1000, target_len=2000, matches=0, block=0, mapq=255 (placeholders).
        format!(
            "{}\t1000\t{}\t{}\t{}\t{}\t2000\t{}\t{}\t0\t0\t255\tcg:Z:{}",
            qname, qs, qe, strand, tname, ts, te, cigar
        )
    }

    #[test]
    fn test_build_empty() -> Result<()> {
        let idx = PafQueryIndex::build(Cursor::new(""))?;
        assert!(idx.is_empty());
        assert_eq!(idx.names.len(), 0);
        Ok(())
    }

    #[test]
    fn test_build_single_record_plus_strand() -> Result<()> {
        let paf = paf_line("qry1", 0, 100, "+", "ref1", 50, 150, "100=");
        let idx = PafQueryIndex::build(Cursor::new(paf))?;
        assert!(!idx.is_empty());
        assert_eq!(idx.names.len(), 1);
        let qid = idx.query_id("qry1").unwrap();
        let hits = idx.query(qid, 0, 100);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].strand, '+');
        assert_eq!(hits[0].target_name, "ref1");
        assert_eq!(hits[0].target_start, 50);
        assert_eq!(hits[0].target_end, 150);
        assert_eq!(hits[0].cigar.len(), 1);
        assert_eq!(hits[0].cigar[0].op(), '=');
        Ok(())
    }

    #[test]
    fn test_build_minus_strand_record_inserted() -> Result<()> {
        // Minus-strand records must be inserted (unlike PafIndex.reverse_trees
        // which only holds '+' strand mirror entries).
        let paf = paf_line("qry1", 0, 100, "-", "ref1", 50, 150, "100=");
        let idx = PafQueryIndex::build(Cursor::new(paf))?;
        assert!(!idx.is_empty());
        let qid = idx.query_id("qry1").unwrap();
        let hits = idx.query(qid, 0, 100);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].strand, '-');
        Ok(())
    }

    #[test]
    fn test_build_skips_records_without_cigar() -> Result<()> {
        // Record without cg:Z: tag (only a generic tag).
        let line = "qry1\t1000\t0\t100\t+\tref1\t2000\t50\t150\t0\t0\t255\tgi:f:0.9";
        let idx = PafQueryIndex::build(Cursor::new(line))?;
        assert!(idx.is_empty()); // No CIGAR → skipped → empty.
        Ok(())
    }

    #[test]
    fn test_build_multiple_records_multiple_queries() -> Result<()> {
        let paf = format!(
            "{}\n{}\n",
            paf_line("qry1", 0, 100, "+", "ref1", 0, 100, "100="),
            paf_line("qry2", 0, 50, "+", "ref2", 0, 50, "50=")
        );
        let idx = PafQueryIndex::build(Cursor::new(paf))?;
        assert_eq!(idx.names.len(), 2);
        assert!(idx.query_id("qry1").is_some());
        assert!(idx.query_id("qry2").is_some());
        assert!(idx.query_id("qry3").is_none());
        Ok(())
    }

    #[test]
    fn test_query_multiple_alignments_overlap() -> Result<()> {
        // Two alignments on qry1: [0,100) and [50,200). Query [60,80) hits both.
        let paf = format!(
            "{}\n{}\n",
            paf_line("qry1", 0, 100, "+", "ref1", 0, 100, "100="),
            paf_line("qry1", 50, 200, "+", "ref2", 0, 150, "150=")
        );
        let idx = PafQueryIndex::build(Cursor::new(paf))?;
        let qid = idx.query_id("qry1").unwrap();
        let hits = idx.query(qid, 60, 80);
        assert_eq!(hits.len(), 2);
        Ok(())
    }

    #[test]
    fn test_query_non_overlapping() -> Result<()> {
        let paf = paf_line("qry1", 0, 100, "+", "ref1", 0, 100, "100=");
        let idx = PafQueryIndex::build(Cursor::new(paf))?;
        let qid = idx.query_id("qry1").unwrap();
        let hits = idx.query(qid, 200, 300);
        assert!(hits.is_empty());
        Ok(())
    }

    #[test]
    fn test_query_unknown_query_id() {
        let idx = PafQueryIndex {
            names: IndexMap::new(),
            trees: HashMap::new(),
        };
        let hits = idx.query(999, 0, 100);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_build_skips_malformed_lines() -> Result<()> {
        // 1 malformed line + 2 valid lines → 2 records indexed, no error.
        let valid1 = paf_line("qry1", 0, 100, "+", "ref1", 0, 100, "100=");
        let valid2 = paf_line("qry1", 100, 200, "+", "ref1", 100, 200, "100=");
        let malformed = "not_a_paf_line_at_all";
        let paf = format!("{}\n{}\n{}\n", valid1, malformed, valid2);
        let idx = PafQueryIndex::build(Cursor::new(paf))?;
        assert_eq!(idx.names.len(), 1);
        let qid = idx.query_id("qry1").unwrap();
        let hits = idx.query(qid, 0, 200);
        assert_eq!(hits.len(), 2);
        Ok(())
    }

    #[test]
    fn test_build_all_malformed_bails() {
        // All lines malformed → bail (non-PAF format).
        let paf = "garbage_line_one\ngarbage_line_two\n";
        match PafQueryIndex::build(Cursor::new(paf)) {
            Ok(_) => panic!("expected error for all-malformed PAF"),
            Err(e) => assert!(
                e.to_string().contains("no valid records"),
                "unexpected error: {}",
                e
            ),
        }
    }
}
