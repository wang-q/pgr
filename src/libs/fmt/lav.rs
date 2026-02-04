use crate::libs::psl::Psl;
use std::io::{self, BufRead};

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
                if parts.len() >= 2 {
                    if let Ok(s) = parts[1].parse::<i32>() {
                        // C code: score = lineFileNeedNum(lf, words, 1) - 1;
                        current_score = s - 1;
                    }
                }
            } else if line.starts_with('l') {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 6 {
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
    if let Some(idx) = s.find(':') {
        s = &s[..idx];
    }

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
        let match_cnt = (len * block.percent_id as u32 + 50) / 100;
        let mismatch_cnt = len - match_cnt;

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
    for i in 0..blocks.len() - 1 {
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
                assert_eq!(is_rc, false);
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

        // 4. m stanza (Unknown)
        let stanza = reader.next_stanza().unwrap().unwrap();
        match stanza {
            LavStanza::Unknown(line) => {
                assert!(line.contains("m {"));
            }
            _ => panic!("Expected Unknown stanza, got {:?}", stanza),
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
}
