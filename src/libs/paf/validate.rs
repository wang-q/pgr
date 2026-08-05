//! PAF self-consistency validation from the `cg:Z:` CIGAR tag.
//!
//! For each record, the expected query/target end position is reconstructed
//! from its CIGAR (matches + mismatches + insertion/deletion bases) and
//! compared with the declared coordinate. Disagreements are flagged as
//! invalid. Records without a usable `cg:Z:` tag are counted and skipped.

use std::io::Write;

use super::cigar::{cigar_stats, extract_cigar};
use super::record::PafRecord;

/// Report of PAF records whose declared end coordinates disagree with the
/// span implied by their `cg:Z:` CIGAR tag.
#[derive(Debug, Clone, Default)]
pub struct ValidationReport {
    /// Total number of records examined.
    pub total: usize,
    /// Records with a query end not explained by the CIGAR.
    pub query_invalid: usize,
    /// Records with a target end not explained by the CIGAR.
    pub target_invalid: usize,
    /// Records lacking a usable `cg:Z:` tag (coordinate check skipped).
    pub missing_cigar: usize,
    /// Records with a malformed `cg:Z:` tag (coordinate check skipped).
    pub malformed_cigar: usize,
    /// Identifiers of records with a bad query end (`name:start-end`).
    query_idents: Vec<String>,
    /// Identifiers of records with a bad target end (`name:start-end`).
    target_idents: Vec<String>,
}

impl ValidationReport {
    /// Validate one record against its CIGAR-derived expected ends.
    ///
    /// Records without a usable `cg:Z:` tag are counted under `missing_cigar`
    /// / `malformed_cigar` and skipped; this never fails.
    pub fn validate(&mut self, rec: &PafRecord) -> anyhow::Result<()> {
        self.total += 1;
        let ops = match extract_cigar(&rec.tags) {
            Ok(ops) if !ops.is_empty() => ops,
            Ok(_) => {
                self.missing_cigar += 1;
                return Ok(());
            }
            Err(_) => {
                self.malformed_cigar += 1;
                return Ok(());
            }
        };
        let s = cigar_stats(&ops);

        let exp_query_end =
            rec.query_start as u64 + s.matches as u64 + s.mismatches as u64 + s.ins_bp as u64;
        if exp_query_end != rec.query_end as u64 {
            self.query_invalid += 1;
            self.query_idents.push(format!(
                "{}:{}-{}",
                rec.query_name, rec.query_start, rec.query_end
            ));
        }

        let exp_target_end =
            rec.target_start as u64 + s.matches as u64 + s.mismatches as u64 + s.del_bp as u64;
        if exp_target_end != rec.target_end as u64 {
            self.target_invalid += 1;
            self.target_idents.push(format!(
                "{}:{}-{}",
                rec.target_name, rec.target_start, rec.target_end
            ));
        }

        Ok(())
    }

    /// Write the plain-text validation report to `writer`.
    pub fn write_report<W: Write>(&self, mut writer: W) -> anyhow::Result<()> {
        writeln!(writer, "Total records: {}", self.total)?;
        writeln!(writer, "Query invalid records: {}", self.query_invalid)?;
        writeln!(writer, "Target invalid records: {}", self.target_invalid)?;
        writeln!(writer, "Records without cg:Z tag: {}", self.missing_cigar)?;
        writeln!(
            writer,
            "Records with malformed cg:Z tag: {}",
            self.malformed_cigar
        )?;
        writeln!(writer, "Query invalid list:")?;
        for id in &self.query_idents {
            writeln!(writer, "{}", id)?;
        }
        writeln!(writer, "Target invalid list:")?;
        for id in &self.target_idents {
            writeln!(writer, "{}", id)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(query_end: u32, target_end: u32, cigar: &str) -> PafRecord {
        PafRecord {
            query_name: "q".into(),
            query_length: 1000,
            query_start: 10,
            query_end,
            strand: '+',
            target_name: "t".into(),
            target_length: 2000,
            target_start: 20,
            target_end,
            matches: 0,
            block_length: 0,
            mapq: 255,
            tags: if cigar.is_empty() {
                vec![]
            } else {
                vec![format!("cg:Z:{cigar}")]
            },
        }
    }

    #[test]
    fn valid_record_is_not_flagged() {
        // 10= then 5I on query, 5D on target.
        let r = rec(10 + 10 + 5, 20 + 10 + 5, "10=5I5D");
        let mut rep = ValidationReport::default();
        rep.validate(&r).unwrap();
        assert_eq!(rep.total, 1);
        assert_eq!(rep.query_invalid, 0);
        assert_eq!(rep.target_invalid, 0);
    }

    #[test]
    fn bad_query_end_is_flagged() {
        // Declared query_end is 5 too short vs the CIGAR span.
        let r = rec(10 + 5, 20 + 10 + 5, "10=5D");
        let mut rep = ValidationReport::default();
        rep.validate(&r).unwrap();
        assert_eq!(rep.query_invalid, 1);
        assert_eq!(rep.target_invalid, 0);
        assert_eq!(rep.query_idents, vec!["q:10-15".to_string()]);
    }

    #[test]
    fn bad_target_end_is_flagged() {
        // Declared target_end is 3 too long vs the CIGAR span.
        let r = rec(10 + 10, 20 + 10 + 8, "10=3D");
        let mut rep = ValidationReport::default();
        rep.validate(&r).unwrap();
        assert_eq!(rep.query_invalid, 0);
        assert_eq!(rep.target_invalid, 1);
        assert_eq!(rep.target_idents, vec!["t:20-38".to_string()]);
    }

    #[test]
    fn missing_cigar_is_counted_not_flagged() {
        let r = rec(10, 20, "");
        let mut rep = ValidationReport::default();
        rep.validate(&r).unwrap();
        assert_eq!(rep.total, 1);
        assert_eq!(rep.missing_cigar, 1);
        assert_eq!(rep.query_invalid, 0);
        assert_eq!(rep.target_invalid, 0);
    }

    #[test]
    fn malformed_cigar_invalid_op_is_counted() {
        // 'N' is not a valid PAF CIGAR op -> parse_cigar errors.
        let r = rec(10, 20, "10N");
        let mut rep = ValidationReport::default();
        rep.validate(&r).unwrap();
        assert_eq!(rep.total, 1);
        assert_eq!(rep.malformed_cigar, 1);
        assert_eq!(rep.query_invalid, 0);
        assert_eq!(rep.target_invalid, 0);
    }

    #[test]
    fn malformed_cigar_trailing_digits_is_counted() {
        // "25M25" ends in digits with no following op -> parse_cigar errors.
        let r = rec(10, 20, "10M5");
        let mut rep = ValidationReport::default();
        rep.validate(&r).unwrap();
        assert_eq!(rep.total, 1);
        assert_eq!(rep.malformed_cigar, 1);
        assert_eq!(rep.query_invalid, 0);
        assert_eq!(rep.target_invalid, 0);
    }

    #[test]
    fn missing_and_malformed_are_counted_separately() {
        let mut rep = ValidationReport::default();
        rep.validate(&rec(10, 20, "")).unwrap(); // missing
        rep.validate(&rec(10, 20, "10Q")).unwrap(); // malformed
        rep.validate(&rec(10 + 5, 20 + 5, "5I5D")).unwrap(); // valid
        assert_eq!(rep.total, 3);
        assert_eq!(rep.missing_cigar, 1);
        assert_eq!(rep.malformed_cigar, 1);
        assert_eq!(rep.query_invalid, 0);
        assert_eq!(rep.target_invalid, 0);
    }

    #[test]
    fn x_and_m_contribute_to_both_axes() {
        // 5= 2X 3I 4D -> query span 5+2+3=10, target span 5+2+4=11.
        let r = rec(10 + 10, 20 + 11, "5=2X3I4D");
        let mut rep = ValidationReport::default();
        rep.validate(&r).unwrap();
        assert_eq!(rep.query_invalid, 0);
        assert_eq!(rep.target_invalid, 0);
    }

    #[test]
    fn report_renders_counts_and_lists() {
        let mut rep = ValidationReport::default();
        rep.validate(&rec(10, 20, "")).unwrap();
        rep.validate(&rec(10 + 5, 20 + 10, "10=5X")).unwrap();
        let mut buf = Vec::new();
        rep.write_report(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("Total records: 2"));
        assert!(text.contains("Records without cg:Z tag: 1"));
        assert!(text.contains("Query invalid list:"));
        assert!(text.contains("q:10-15"));
    }
}
