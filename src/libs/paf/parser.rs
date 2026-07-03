/// PAF format parser.
///
/// Uses `csv::Reader` in flexible mode to handle variable column counts
/// (12 mandatory + optional tags), following wgatools' approach.
use super::record::PafRecord;
use anyhow::Context;
use std::io::BufRead;

/// Parse a complete PAF file into a vector of `PafRecord`s.
///
/// Skips empty lines and lines starting with `#`.
pub fn parse_paf<R: BufRead>(reader: R) -> anyhow::Result<Vec<PafRecord>> {
    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let rec = parse_paf_line(&line).with_context(|| format!("invalid PAF line: {line}"))?;
        records.push(rec);
    }
    Ok(records)
}

/// Parse a single PAF line (not starting with `#`, not empty).
pub fn parse_paf_line(line: &str) -> anyhow::Result<PafRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        anyhow::bail!(
            "need at least 12 tab-separated fields, got {}",
            fields.len()
        );
    }

    let query_name = fields[0].to_string();
    let query_length = fields[1]
        .parse::<u32>()
        .map_err(|_| anyhow::anyhow!("invalid query_length: {}", fields[1]))?;
    let query_start = fields[2]
        .parse::<u32>()
        .map_err(|_| anyhow::anyhow!("invalid query_start: {}", fields[2]))?;
    let query_end = fields[3]
        .parse::<u32>()
        .map_err(|_| anyhow::anyhow!("invalid query_end: {}", fields[3]))?;
    let strand = match fields[4] {
        "+" => '+',
        "-" => '-',
        _ => anyhow::bail!("invalid strand: {}", fields[4]),
    };
    let target_name = fields[5].to_string();
    let target_length = fields[6]
        .parse::<u32>()
        .map_err(|_| anyhow::anyhow!("invalid target_length: {}", fields[6]))?;
    let target_start = fields[7]
        .parse::<u32>()
        .map_err(|_| anyhow::anyhow!("invalid target_start: {}", fields[7]))?;
    let target_end = fields[8]
        .parse::<u32>()
        .map_err(|_| anyhow::anyhow!("invalid target_end: {}", fields[8]))?;
    let matches = fields[9]
        .parse::<u32>()
        .map_err(|_| anyhow::anyhow!("invalid matches: {}", fields[9]))?;
    let block_length = fields[10]
        .parse::<u32>()
        .map_err(|_| anyhow::anyhow!("invalid block_length: {}", fields[10]))?;
    let mapq = fields[11]
        .parse::<u8>()
        .map_err(|_| anyhow::anyhow!("invalid mapq: {}", fields[11]))?;

    let tags: Vec<String> = fields[12..].iter().map(|s| s.to_string()).collect();

    Ok(PafRecord {
        query_name,
        query_length,
        query_start,
        query_end,
        strand,
        target_name,
        target_length,
        target_start,
        target_end,
        matches,
        block_length,
        mapq,
        tags,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    // ── Line parsing ──────────────────────────────────────────

    #[test]
    fn test_parse_paf_line_valid() {
        let line = "seq1\t100\t0\t100\t+\tseq2\t200\t50\t150\t80\t100\t255";
        let rec = parse_paf_line(line).unwrap();
        assert_eq!(rec.query_name, "seq1");
        assert_eq!(rec.query_length, 100);
        assert_eq!(rec.query_start, 0);
        assert_eq!(rec.query_end, 100);
        assert_eq!(rec.strand, '+');
        assert_eq!(rec.target_name, "seq2");
        assert_eq!(rec.target_length, 200);
        assert_eq!(rec.target_start, 50);
        assert_eq!(rec.target_end, 150);
        assert_eq!(rec.matches, 80);
        assert_eq!(rec.block_length, 100);
        assert_eq!(rec.mapq, 255);
        assert!(rec.tags.is_empty());
    }

    #[test]
    fn test_parse_paf_line_reverse_strand() {
        let line = "qry\t500\t10\t60\t-\tref\t1000\t100\t150\t45\t50\t255";
        let rec = parse_paf_line(line).unwrap();
        assert_eq!(rec.strand, '-');
    }

    #[test]
    fn test_parse_paf_line_with_tags() {
        let line = "q\t100\t0\t50\t+\tt\t200\t0\t50\t45\t50\t255\tcg:Z:50M\tgi:f:0.9";
        let rec = parse_paf_line(line).unwrap();
        assert_eq!(rec.tags.len(), 2);
        assert_eq!(rec.tags[0], "cg:Z:50M");
        assert_eq!(rec.tags[1], "gi:f:0.9");
    }

    #[test]
    fn test_parse_paf_line_too_few_fields() {
        assert!(parse_paf_line("a\tb\tc").is_err());
    }

    #[test]
    fn test_parse_paf_line_invalid_number() {
        assert!(parse_paf_line("q\txxx\t0\t100\t+\tt\t200\t0\t100\t80\t100\t255").is_err());
    }

    #[test]
    fn test_parse_paf_line_invalid_strand() {
        assert!(parse_paf_line("q\t100\t0\t100\t*\tt\t200\t0\t100\t80\t100\t255").is_err());
    }

    #[test]
    fn test_parse_paf_line_strand_rejects_multi_char() {
        // Strand field must be exactly "+" or "-"; leading-correct but
        // multi-char values (e.g. "+foo", "+-") must be rejected.
        assert!(parse_paf_line("q\t100\t0\t100\t+foo\tt\t200\t0\t100\t80\t100\t255").is_err());
        assert!(parse_paf_line("q\t100\t0\t100\t+-\tt\t200\t0\t100\t80\t100\t255").is_err());
    }

    // ── Full file parsing ─────────────────────────────────────

    #[test]
    fn test_parse_paf_basic() {
        let input = "q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255\n\
                     q2\t300\t10\t60\t-\tt2\t400\t10\t60\t45\t50\t255\n";
        let reader = BufReader::new(input.as_bytes());
        let records = parse_paf(reader).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].query_name, "q1");
        assert_eq!(records[1].query_name, "q2");
    }

    #[test]
    fn test_parse_paf_skip_comments_and_empty() {
        let input = "# header\n\
                     \n\
                     q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255\n\
                     # comment\n\
                     q2\t300\t10\t60\t-\tt2\t400\t10\t60\t45\t50\t255\n";
        let reader = BufReader::new(input.as_bytes());
        let records = parse_paf(reader).unwrap();
        assert_eq!(records.len(), 2);
    }

    #[test]
    fn test_parse_paf_invalid_line_reports_error() {
        let input = "q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255\nbad\n";
        let reader = BufReader::new(input.as_bytes());
        let err = parse_paf(reader).unwrap_err();
        assert!(err.to_string().contains("bad"));
    }

    // ── Roundtrip: parse → write → parse ──────────────────────

    #[test]
    fn test_parse_write_roundtrip() {
        let rec = PafRecord {
            query_name: "q".into(),
            query_length: 100,
            query_start: 0,
            query_end: 50,
            strand: '+',
            target_name: "t".into(),
            target_length: 200,
            target_start: 10,
            target_end: 60,
            matches: 45,
            block_length: 50,
            mapq: 255,
            tags: vec!["cg:Z:50M".into(), "gi:f:0.9".into()],
        };

        // Write to buffer
        let mut buf = Vec::new();
        super::super::record::write_paf_record(&mut buf, &rec).unwrap();

        // Parse back
        let line = String::from_utf8(buf).unwrap();
        let rec2 = parse_paf_line(line.trim()).unwrap();

        assert_eq!(rec.query_name, rec2.query_name);
        assert_eq!(rec.query_start, rec2.query_start);
        assert_eq!(rec.strand, rec2.strand);
        assert_eq!(rec.tags, rec2.tags);
    }
}
