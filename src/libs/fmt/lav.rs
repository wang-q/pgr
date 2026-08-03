use std::io::{self, BufRead, Write};

use crate::libs::fmt::psl::Psl;

#[derive(Debug, Clone, Default)]
pub struct Block {
    pub score: i32,
    pub t_start: i64,
    pub t_end: i64,
    pub q_start: i64,
    pub q_end: i64,
    pub percent_id: i32,
}

#[derive(Debug, Clone)]
pub enum LavStanza {
    Sizes {
        t_size: i64,
        q_size: i64,
    },
    Header {
        t_name: String,
        q_name: String,
        is_rc: bool,
    },
    Data {
        lines: Vec<String>,
    },
    Alignment {
        blocks: Vec<Block>,
    },
    /// The `m { ... }` mask stanza (lastz writes it; ignored by converters).
    Mask,
    Unknown(String),
}

pub struct LavReader<R: BufRead> {
    lines: std::iter::Peekable<std::io::Lines<R>>,
}

impl<R: BufRead> LavReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            lines: reader.lines().peekable(),
        }
    }

    pub fn next_stanza(&mut self) -> io::Result<Option<LavStanza>> {
        while let Some(line_res) = self.lines.next() {
            let line = line_res?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if line.starts_with("s {") {
                return Ok(Some(self.parse_s()?));
            } else if line.starts_with("h {") {
                return Ok(Some(self.parse_h()?));
            } else if line.starts_with("d {") {
                return Ok(Some(self.parse_d()?));
            } else if line.starts_with("a {") {
                return Ok(Some(self.parse_a()?));
            } else if line.starts_with("m {") {
                self.skip_stanza()?;
                return Ok(Some(LavStanza::Mask));
            } else if line.ends_with('{') {
                self.skip_stanza()?;
                return Ok(Some(LavStanza::Unknown(line.to_string())));
            }
        }
        Ok(None)
    }

    fn parse_s(&mut self) -> io::Result<LavStanza> {
        let t_size = self.read_size_line()?;
        let q_size = self.read_size_line()?;
        self.skip_until_brace()?;
        Ok(LavStanza::Sizes { t_size, q_size })
    }

    fn read_size_line(&mut self) -> io::Result<i64> {
        loop {
            let line = self.read_line()?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 3 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid s stanza line (expected >= 3 words): {}", line),
                ));
            }
            return parts[2]
                .parse::<i64>()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e));
        }
    }

    fn parse_h(&mut self) -> io::Result<LavStanza> {
        let mut t_name = String::new();
        let mut q_name = String::new();
        let mut is_rc = false;

        let mut i = 0;
        loop {
            let line = self.read_line()?;
            let line_trim = line.trim();
            if line_trim == "}" {
                break;
            }
            if line_trim.starts_with('#') || line_trim.is_empty() {
                continue;
            }

            let parts: Vec<&str> = line_trim.split_whitespace().collect();
            if parts.is_empty() {
                continue;
            }

            let word = parts[0];
            let content = parse_header_word(word);

            if i == 0 {
                t_name = content;
            } else if i == 1 {
                q_name = content;
            }

            if line.contains("(reverse") {
                is_rc = true;
            }

            i += 1;
        }

        Ok(LavStanza::Header {
            t_name,
            q_name,
            is_rc,
        })
    }

    fn parse_d(&mut self) -> io::Result<LavStanza> {
        let mut lines = Vec::new();
        loop {
            let line = self.read_line()?;
            if line.trim() == "}" {
                break;
            }
            lines.push(line);
        }
        Ok(LavStanza::Data { lines })
    }

    fn parse_a(&mut self) -> io::Result<LavStanza> {
        let mut blocks = Vec::new();
        let mut current_score = 0;

        loop {
            let line = self.read_line()?;
            if line.trim() == "}" {
                break;
            }
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if line.starts_with('s') {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() < 2 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("malformed s line (need >=2 fields): {}", line),
                    ));
                }
                let s: i32 = parts[1].parse().map_err(|e| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid score in s line: {}", e),
                    )
                })?;
                // C code: score = lineFileNeedNum(lf, words, 1) - 1;
                current_score = s - 1;
            } else if line.starts_with('l') {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() < 6 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("malformed l line (need >=6 fields): {}", line),
                    ));
                }
                let t_start = parts[1]
                    .parse::<i64>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
                    - 1;
                let q_start = parts[2]
                    .parse::<i64>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
                    - 1;
                let t_end = parts[3]
                    .parse::<i64>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let q_end = parts[4]
                    .parse::<i64>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let percent_id = parts[5]
                    .parse::<i32>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                if (q_end - q_start) != (t_end - t_start) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Block size mismatch",
                    ));
                }

                if q_end == q_start && t_end == t_start {
                    continue;
                }

                blocks.push(Block {
                    score: current_score,
                    t_start,
                    t_end,
                    q_start,
                    q_end,
                    percent_id,
                });
            }
        }

        blocks = remove_frayed_ends(blocks);

        Ok(LavStanza::Alignment { blocks })
    }

    fn read_line(&mut self) -> io::Result<String> {
        if let Some(res) = self.lines.next() {
            res
        } else {
            Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Unexpected EOF",
            ))
        }
    }

    fn skip_until_brace(&mut self) -> io::Result<()> {
        loop {
            let line = self.read_line()?;
            if line.trim() == "}" {
                break;
            }
        }
        Ok(())
    }

    fn skip_stanza(&mut self) -> io::Result<()> {
        self.skip_until_brace()
    }
}

fn parse_header_word(word: &str) -> String {
    let mut s = word;
    if s.starts_with('"') {
        s = &s[1..];
    }
    if s.starts_with('>') {
        s = &s[1..];
    }

    if let Some(idx) = s.find('"') {
        s = &s[..idx];
    }

    // Remove range specifiers (e.g., :start-end)
    // if let Some(idx) = s.find(':') {
    //     s = &s[..idx];
    // }

    // Extract filename from path
    if let Some(idx) = s.rfind('/') {
        s = &s[idx + 1..];
    }

    // Remove common extensions
    if let Some(stripped) = s.strip_suffix(".nib") {
        s = stripped;
    } else if let Some(stripped) = s.strip_suffix(".fa") {
        s = stripped;
    } else if let Some(stripped) = s.strip_suffix(".fasta") {
        s = stripped;
    } else if let Some(stripped) = s.strip_suffix(".2bit") {
        s = stripped;
    }

    s.to_string()
}

fn remove_frayed_ends(mut blocks: Vec<Block>) -> Vec<Block> {
    while !blocks.is_empty() && blocks[0].q_start == blocks[0].q_end {
        blocks.remove(0);
    }
    while !blocks.is_empty() && blocks[blocks.len() - 1].q_start == blocks[blocks.len() - 1].q_end {
        blocks.pop();
    }
    blocks
}

/// Parse a LAV `d` stanza into UCSC-style `##` metadata comment lines.
///
/// Mirrors UCSC `lavToPsl`'s `parseD` + `axtScoreSchemeReadLf` +
/// `axtScoreSchemeDnaWrite`: when the first line mentions `lastz`, emit
/// `##aligner`, `##matrix`, `##gapPenalties`, and `##blastzParms` lines.
/// Returns an empty vec for non-lastz d stanzas (matching UCSC, which only
/// processes lastz-style stanzas).
pub fn parse_d_stanza_to_comments(lines: &[String]) -> anyhow::Result<Vec<String>> {
    if lines.is_empty() {
        return Ok(Vec::new());
    }

    // UCSC parseD: stripChar(line, '"') then chopLine. Only proceeds if
    // stringIn("lastz", line).
    let first_stripped: String = lines[0].chars().filter(|c| *c != '"').collect();
    if !first_stripped.contains("lastz") {
        return Ok(Vec::new());
    }

    let words: Vec<&str> = first_stripped.split_whitespace().collect();
    if words.is_empty() {
        return Ok(Vec::new());
    }
    let aligner = words[0];

    // ##aligner=<words[0]> + " <words[3]> " for each remaining param. The
    // surrounding spaces match UCSC's `fprintf(f, " %s ", words[i])` format,
    // producing double spaces between params and a trailing space.
    let mut aligner_line = format!("##aligner={}", aligner);
    for w in &words[3..] {
        aligner_line.push(' ');
        aligner_line.push_str(w);
        aligner_line.push(' ');
    }
    let mut comments = vec![aligner_line];

    // Locate the matrix header line ("A C G T"). UCSC axtScoreSchemeReadLf
    // scans for the row where row[0..4] == ['A','C','G','T'].
    let mut header_idx = 1;
    while header_idx < lines.len() {
        let parts: Vec<&str> = lines[header_idx].split_whitespace().collect();
        if parts.len() >= 4
            && parts[0] == "A"
            && parts[1] == "C"
            && parts[2] == "G"
            && parts[3] == "T"
        {
            break;
        }
        header_idx += 1;
    }
    if header_idx + 5 > lines.len() {
        anyhow::bail!(
            "lastz d stanza missing matrix rows or params line (header at {}, need {} lines, have {})",
            header_idx,
            header_idx + 6,
            lines.len()
        );
    }

    // Read 4 matrix rows. UCSC skips the first column when wordCount == 5
    // (row-labeled form like "A 91 -114 ...").
    let mut matrix_vals: Vec<i32> = Vec::with_capacity(16);
    for row_i in 0..4 {
        let parts: Vec<&str> = lines[header_idx + 1 + row_i].split_whitespace().collect();
        let start = if parts.len() == 5 { 1 } else { 0 };
        if start + 4 > parts.len() {
            anyhow::bail!(
                "matrix row {} has too few columns: {:?}",
                row_i,
                lines[header_idx + 1 + row_i]
            );
        }
        for part in parts.iter().skip(start).take(4) {
            let v: i32 = part
                .parse()
                .map_err(|e| anyhow::anyhow!("invalid matrix value {:?}: {}", part, e))?;
            matrix_vals.push(v);
        }
    }

    // Params line: "O = 400, E = 30, K = 3000, L = 3000, M = 0"
    let params_line = &lines[header_idx + 5];
    let param_tokens: Vec<&str> = params_line
        .split([' ', '=', ',', '\t'])
        .filter(|s| !s.is_empty())
        .collect();
    let mut gap_open: i32 = 400;
    let mut gap_extend: i32 = 30;
    let mut k = 0;
    while k + 1 < param_tokens.len() {
        let val = param_tokens[k + 1].trim_end_matches('"');
        match param_tokens[k] {
            "O" => {
                gap_open = val.parse().unwrap_or(gap_open);
            }
            "E" => {
                gap_extend = val.parse().unwrap_or(gap_extend);
            }
            _ => {}
        }
        k += 2;
    }

    // ##matrix=<aligner> 16 v1,v2,...,v16
    let matrix_str: Vec<String> = matrix_vals.iter().map(|v| v.to_string()).collect();
    comments.push(format!("##matrix={} 16 {}", aligner, matrix_str.join(",")));

    // ##gapPenalties=<aligner> O=<open> E=<extend>
    comments.push(format!(
        "##gapPenalties={} O={} E={}",
        aligner, gap_open, gap_extend
    ));

    // ##blastzParms=<extra> — UCSC strips spaces and quotes from the params
    // line (the `extra` field accumulates the line, then stripChar ' ' and '"').
    let extra: String = params_line
        .chars()
        .filter(|c| *c != ' ' && *c != '"')
        .collect();
    comments.push(format!("##blastzParms={}", extra));

    Ok(comments)
}

/// Convert LAV alignment blocks into a Psl record.
pub fn blocks_to_psl(
    blocks: &[Block],
    t_size: u32,
    q_size: u32,
    t_name: &str,
    q_name: &str,
    strand: &str,
) -> Psl {
    let mut psl = Psl::new();
    psl.t_size = t_size;
    psl.q_size = q_size;
    psl.t_name = t_name.to_string();
    psl.q_name = q_name.to_string();
    psl.strand = strand.to_string();

    // Calculate overall range and stats
    let mut q_min = i64::MAX;
    let mut q_max = i64::MIN;
    let mut t_min = i64::MAX;
    let mut t_max = i64::MIN;

    for block in blocks {
        let len = (block.t_end - block.t_start) as u32;
        // UCSC lavToPsl calculation: match = (width * identity + 50)/100
        let match_cnt = (len.saturating_mul(block.percent_id as u32) + 50) / 100;
        let match_cnt = match_cnt.min(len);
        let mismatch_cnt = len.saturating_sub(match_cnt);

        psl.match_count += match_cnt;
        psl.mismatch_count += mismatch_cnt;

        psl.block_count += 1;
        psl.block_sizes.push(len);
        psl.q_starts.push(block.q_start as u32);
        psl.t_starts.push(block.t_start as u32);

        if block.q_start < q_min {
            q_min = block.q_start;
        }
        if block.q_end > q_max {
            q_max = block.q_end;
        }
        if block.t_start < t_min {
            t_min = block.t_start;
        }
        if block.t_end > t_max {
            t_max = block.t_end;
        }
    }

    if !blocks.is_empty() {
        if strand == "-" {
            psl.q_start = (q_size as i64 - q_max) as i32;
            psl.q_end = (q_size as i64 - q_min) as i32;
        } else {
            psl.q_start = q_min as i32;
            psl.q_end = q_max as i32;
        }
        psl.t_start = t_min as i32;
        psl.t_end = t_max as i32;
    }

    // Gaps (inserts)
    for i in 0..blocks.len().saturating_sub(1) {
        let curr = &blocks[i];
        let next = &blocks[i + 1];

        // Assumption: blocks are sorted by T. LAV usually implies this.
        // If not, gap calculation might be weird (negative).
        // Let's assume non-negative gaps for now, or clamp to 0.

        let q_gap = next.q_start - curr.q_end;
        let t_gap = next.t_start - curr.t_end;

        if q_gap > 0 {
            psl.q_num_insert += 1;
            psl.q_base_insert += q_gap as i32;
        }

        if t_gap > 0 {
            psl.t_num_insert += 1;
            psl.t_base_insert += t_gap as i32;
        }
    }

    psl
}

/// Convert a LAV stream to PSL with optional target strand annotation.
///
/// Iterates LAV stanzas, accumulating sizes/header state, and emits one PSL
/// per Alignment stanza. Unknown stanzas trigger a `log::warn!` unless
/// `strict` is set, in which case they bail.
pub fn lav_to_psl<R: BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    target_strand: Option<&str>,
    strict: bool,
) -> anyhow::Result<()> {
    let mut lav_reader = LavReader::new(reader);

    let mut t_size: Option<u32> = None;
    let mut q_size: Option<u32> = None;
    let mut t_name: Option<String> = None;
    let mut q_name: Option<String> = None;
    let mut strand: Option<String> = None;

    while let Some(stanza) = lav_reader.next_stanza()? {
        match stanza {
            LavStanza::Sizes {
                t_size: t,
                q_size: q,
            } => {
                t_size =
                    Some(u32::try_from(t).map_err(|_| anyhow::anyhow!("invalid t_size: {}", t))?);
                q_size =
                    Some(u32::try_from(q).map_err(|_| anyhow::anyhow!("invalid q_size: {}", q))?);
            }
            LavStanza::Header {
                t_name: t,
                q_name: q,
                is_rc,
            } => {
                t_name = Some(t);
                q_name = Some(q);
                strand = Some(if is_rc {
                    "-".to_string()
                } else {
                    "+".to_string()
                });
            }
            LavStanza::Alignment { blocks } => {
                if blocks.is_empty() {
                    continue;
                }

                let t_size = t_size.ok_or_else(|| {
                    anyhow::anyhow!("Alignment stanza encountered before Sizes stanza")
                })?;
                let q_size = q_size.ok_or_else(|| {
                    anyhow::anyhow!("Alignment stanza encountered before Sizes stanza")
                })?;
                let t_name = t_name.as_deref().ok_or_else(|| {
                    anyhow::anyhow!("Alignment stanza encountered before Header stanza")
                })?;
                let q_name = q_name.as_deref().ok_or_else(|| {
                    anyhow::anyhow!("Alignment stanza encountered before Header stanza")
                })?;
                let strand = strand.as_deref().ok_or_else(|| {
                    anyhow::anyhow!("Alignment stanza encountered before Header stanza")
                })?;

                let mut psl = blocks_to_psl(&blocks, t_size, q_size, t_name, q_name, strand);

                if let Some(ts) = target_strand {
                    // Append target strand if provided
                    if psl.strand.len() == 1 {
                        let ts_char = ts
                            .chars()
                            .next()
                            .ok_or_else(|| anyhow::anyhow!("--target-strand cannot be empty"))?;
                        psl.strand.push(ts_char);
                    }
                }

                psl.write_to(writer)?;
            }
            LavStanza::Data { lines } => {
                // UCSC lavToPsl parseD: emit ##aligner/##matrix/##gapPenalties/
                // ##blastzParms metadata lines from the d stanza.
                let comments = parse_d_stanza_to_comments(&lines)?;
                for comment in &comments {
                    writeln!(writer, "{}", comment)?;
                }
            }
            LavStanza::Mask => {}
            other => {
                if strict {
                    anyhow::bail!("unknown lav stanza: {:?}", other);
                }
                log::warn!("skipping unknown lav stanza: {:?}", other);
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_lav_simple() {
        let data = r#"#:lav
s {
    "/path/target.fa" 1 1000
    "/path/query.fa" 1 500
}
h {
    ">target.fa"
    ">query.fa"
}
a {
    s 100
    l 1 1 10 10 95
}
m {
    x y z
}
"#;
        let mut reader = LavReader::new(Cursor::new(data));

        // 1. s stanza
        let stanza = reader.next_stanza().unwrap().unwrap();
        match stanza {
            LavStanza::Sizes { t_size, q_size } => {
                assert_eq!(t_size, 1000);
                assert_eq!(q_size, 500);
            }
            _ => panic!("Expected Sizes stanza, got {:?}", stanza),
        }

        // 2. h stanza
        let stanza = reader.next_stanza().unwrap().unwrap();
        match stanza {
            LavStanza::Header {
                t_name,
                q_name,
                is_rc,
            } => {
                assert_eq!(t_name, "target");
                assert_eq!(q_name, "query");
                assert!(!is_rc);
            }
            _ => panic!("Expected Header stanza, got {:?}", stanza),
        }

        // 3. a stanza
        let stanza = reader.next_stanza().unwrap().unwrap();
        match stanza {
            LavStanza::Alignment { blocks } => {
                assert_eq!(blocks.len(), 1);
                assert_eq!(blocks[0].score, 99); // 100 - 1
                assert_eq!(blocks[0].t_start, 0); // 1 - 1
                assert_eq!(blocks[0].q_start, 0); // 1 - 1
                assert_eq!(blocks[0].t_end, 10);
                assert_eq!(blocks[0].q_end, 10);
            }
            _ => panic!("Expected Alignment stanza, got {:?}", stanza),
        }

        // 4. m stanza (mask, ignored)
        let stanza = reader.next_stanza().unwrap().unwrap();
        match stanza {
            LavStanza::Mask => {}
            _ => panic!("Expected Mask stanza, got {:?}", stanza),
        }

        // End
        assert!(reader.next_stanza().unwrap().is_none());
    }

    #[test]
    fn test_parse_lav_rc() {
        let data = r#"
h {
    ">target"
    ">query" (reverse)
}
"#;
        let mut reader = LavReader::new(Cursor::new(data));
        let stanza = reader.next_stanza().unwrap().unwrap();
        match stanza {
            LavStanza::Header { is_rc, .. } => {
                assert!(is_rc);
            }
            _ => panic!("Expected Header stanza, got {:?}", stanza),
        }
    }

    #[test]
    fn test_parse_d_stanza_to_comments_lastz() {
        // Mirrors tests/pgr/cmp/shared/lastz.lav d stanza (UCSC reference).
        let lines = vec![
            r#"  "lastz.v1.04.41 tests/pgr/pseudocat.fa tests/pgr/pseudopig.fa "#.to_string(),
            "     A    C    G    T".to_string(),
            "    91 -114  -31 -123".to_string(),
            "  -114  100 -125  -31".to_string(),
            "   -31 -125  100 -114".to_string(),
            "  -123  -31 -114   91".to_string(),
            r#"  O = 400, E = 30, K = 3000, L = 3000, M = 0""#.to_string(),
        ];
        let comments = parse_d_stanza_to_comments(&lines).unwrap();
        assert_eq!(comments.len(), 4);
        // No params beyond word[2] -> ##aligner has no trailing params.
        assert_eq!(comments[0], "##aligner=lastz.v1.04.41");
        assert_eq!(
            comments[1],
            "##matrix=lastz.v1.04.41 16 91,-114,-31,-123,-114,100,-125,-31,-31,-125,100,-114,-123,-31,-114,91"
        );
        assert_eq!(comments[2], "##gapPenalties=lastz.v1.04.41 O=400 E=30");
        assert_eq!(comments[3], "##blastzParms=O=400,E=30,K=3000,L=3000,M=0");
    }

    #[test]
    fn test_parse_d_stanza_to_comments_with_params() {
        // Mirrors tests/lav/newStyleLastz.lav d stanza with extra params.
        let lines = vec![
            r#"  "lastz.v1.03.46 hg19.chrM.fa susScr3.chrM.fa M=50 T=2 O=400 E=30 Q=hg19.susScr3.chrM.Q.txt --output=hg19.susScr3.chrM.lav "#.to_string(),
            "     A    C    G    T".to_string(),
            "    79  -84  -55 -128".to_string(),
            "   -84  100 -174  -55".to_string(),
            "   -55 -174  100  -84".to_string(),
            "  -128  -55  -84   79".to_string(),
            r#"  O = 400, E = 30, K = 3000, L = 3000, M = 50""#.to_string(),
        ];
        let comments = parse_d_stanza_to_comments(&lines).unwrap();
        assert_eq!(comments.len(), 4);
        // Words[3..] each get surrounding spaces -> double space between params.
        assert_eq!(
            comments[0],
            "##aligner=lastz.v1.03.46 M=50  T=2  O=400  E=30  Q=hg19.susScr3.chrM.Q.txt  --output=hg19.susScr3.chrM.lav "
        );
        assert_eq!(
            comments[1],
            "##matrix=lastz.v1.03.46 16 79,-84,-55,-128,-84,100,-174,-55,-55,-174,100,-84,-128,-55,-84,79"
        );
        assert_eq!(comments[2], "##gapPenalties=lastz.v1.03.46 O=400 E=30");
        assert_eq!(comments[3], "##blastzParms=O=400,E=30,K=3000,L=3000,M=50");
    }

    #[test]
    fn test_parse_d_stanza_to_comments_non_lastz() {
        // UCSC parseD only proceeds when stringIn("lastz", line). An aligner
        // name without the "lastz" substring yields no comments.
        let lines = vec!["  \"blurz.v1 target.fa query.fa".to_string()];
        let comments = parse_d_stanza_to_comments(&lines).unwrap();
        assert!(comments.is_empty());
    }

    #[test]
    fn test_parse_d_stanza_to_comments_blastz_v7() {
        // "blastz" contains "lastz" as a substring (b-lastz), so UCSC's
        // stringIn("lastz", ...) matches and comments are emitted.
        let lines = vec![
            r#"  "blastz.v7 hg19.chrM.fa susScr3.chrM.fa M=50 T=2 O=400 E=30 Q=hg19.susScr3.blastz.q"#.to_string(),
            "     A    C    G    T".to_string(),
            "    79  -84  -55 -128".to_string(),
            "   -84  100 -174  -55".to_string(),
            "   -55 -174  100  -84".to_string(),
            "  -128  -55  -84   79".to_string(),
            r#"  O = 400, E = 30, K = 3000, L = 3000, M = 50""#.to_string(),
        ];
        let comments = parse_d_stanza_to_comments(&lines).unwrap();
        assert_eq!(comments.len(), 4);
        assert_eq!(
            comments[0],
            "##aligner=blastz.v7 M=50  T=2  O=400  E=30  Q=hg19.susScr3.blastz.q "
        );
        assert_eq!(
            comments[1],
            "##matrix=blastz.v7 16 79,-84,-55,-128,-84,100,-174,-55,-55,-174,100,-84,-128,-55,-84,79"
        );
        assert_eq!(comments[2], "##gapPenalties=blastz.v7 O=400 E=30");
        assert_eq!(comments[3], "##blastzParms=O=400,E=30,K=3000,L=3000,M=50");
    }
}
