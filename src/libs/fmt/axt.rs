use std::io::BufRead;

#[derive(Debug, Clone, Default)]
pub struct Axt {
    pub id: u64, // The first number in the header line
    pub t_name: String,
    pub t_start: usize, // 0-based
    pub t_end: usize,   // 0-based, half-open
    pub t_strand: char,
    pub q_name: String,
    pub q_start: usize, // 0-based
    pub q_end: usize,   // 0-based, half-open
    pub q_strand: char,
    pub score: Option<i32>,
    pub t_sym: String,
    pub q_sym: String,
}

impl Axt {
    pub fn new() -> Self {
        Self::default()
    }
}

pub struct AxtReader<R> {
    reader: std::io::BufReader<R>,
    line_buf: String,
}

impl<R: std::io::Read> AxtReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            reader: std::io::BufReader::new(inner),
            line_buf: String::new(),
        }
    }

    fn read_line(&mut self) -> std::io::Result<usize> {
        self.line_buf.clear();
        self.reader.read_line(&mut self.line_buf)
    }
}

impl<R: std::io::Read> Iterator for AxtReader<R> {
    type Item = anyhow::Result<Axt>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Read first line (header)
            match self.read_line() {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    let line = self.line_buf.trim();
                    if line.is_empty() || line.starts_with('#') {
                        continue;
                    }

                    // Parse header
                    // Format: id tName tStart tEnd tStrand? qName qStart qEnd qStrand score?
                    // Example: 0 chr19 3001012 3001075 chr11 70568380 70568443 - 3500
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() < 8 {
                        return Some(Err(anyhow::anyhow!("Invalid AXT header: {}", line)));
                    }

                    let id: u64 = match parts[0].parse() {
                        Ok(v) => v,
                        Err(_) => {
                            return Some(Err(anyhow::anyhow!("Invalid AXT id: {}", parts[0])))
                        }
                    };

                    let t_name = parts[1].to_string();
                    let t_start: usize = match parts[2].parse::<usize>() {
                        Ok(v) => {
                            if v > 0 {
                                v - 1
                            } else {
                                0
                            }
                        } // 1-based to 0-based
                        Err(_) => {
                            return Some(Err(anyhow::anyhow!("Invalid tStart: {}", parts[2])))
                        }
                    };
                    let t_end: usize = match parts[3].parse::<usize>() {
                        Ok(v) => v,
                        Err(_) => return Some(Err(anyhow::anyhow!("Invalid tEnd: {}", parts[3]))),
                    };

                    let q_name = parts[4].to_string();
                    let q_start: usize = match parts[5].parse::<usize>() {
                        Ok(v) => {
                            if v > 0 {
                                v - 1
                            } else {
                                0
                            }
                        } // 1-based to 0-based
                        Err(_) => {
                            return Some(Err(anyhow::anyhow!("Invalid qStart: {}", parts[5])))
                        }
                    };
                    let q_end: usize = match parts[6].parse::<usize>() {
                        Ok(v) => v,
                        Err(_) => return Some(Err(anyhow::anyhow!("Invalid qEnd: {}", parts[6]))),
                    };

                    let q_strand = parts[7].chars().next().unwrap_or('?');

                    let score = if parts.len() > 8 {
                        parts[8].parse().ok()
                    } else {
                        None
                    };

                    let mut axt = Axt {
                        id,
                        t_name,
                        t_start,
                        t_end,
                        t_strand: '+',
                        q_name,
                        q_start,
                        q_end,
                        q_strand,
                        score,
                        t_sym: String::new(),
                        q_sym: String::new(),
                    };

                    // Read tSym
                    match self.read_line() {
                        Ok(0) => return Some(Err(anyhow::anyhow!("Unexpected EOF after header"))),
                        Ok(_) => axt.t_sym = self.line_buf.trim().to_string(),
                        Err(e) => return Some(Err(anyhow::Error::new(e))),
                    }

                    // Read qSym
                    match self.read_line() {
                        Ok(0) => return Some(Err(anyhow::anyhow!("Unexpected EOF after tSym"))),
                        Ok(_) => axt.q_sym = self.line_buf.trim().to_string(),
                        Err(e) => return Some(Err(anyhow::Error::new(e))),
                    }

                    if axt.t_sym.len() != axt.q_sym.len() {
                        return Some(Err(anyhow::anyhow!(
                            "Alignment lengths differ: {} vs {}",
                            axt.t_sym.len(),
                            axt.q_sym.len()
                        )));
                    }

                    return Some(Ok(axt));
                }
                Err(e) => return Some(Err(anyhow::Error::new(e))),
            }
        }
    }
}

pub fn write_axt<W: std::io::Write>(writer: &mut W, axt: &Axt) -> std::io::Result<()> {
    // Write header
    // Convert 0-based to 1-based
    let t_start = axt.t_start + 1;
    let t_end = axt.t_end;
    let q_start = axt.q_start + 1;
    let q_end = axt.q_end;

    let score_str = match axt.score {
        Some(s) => format!(" {}", s),
        None => String::new(),
    };

    writeln!(
        writer,
        "{} {} {} {} {} {} {} {} {}",
        axt.id,
        axt.t_name,
        t_start,
        t_end,
        axt.q_name,
        q_start,
        q_end,
        axt.q_strand,
        score_str // Placeholder
    )?;

    writeln!(writer, "{}", axt.t_sym)?;
    writeln!(writer, "{}", axt.q_sym)?;
    writeln!(writer)?; // Blank line

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_axt() {
        let input = "\
0 chr19 3001012 3001075 chr11 70568380 70568443 - 3500
TCAGCTCTGGACGCAATACTCGCTCCCAGTCCAGATTCCTTCCTGATACTCGTCATGTGAGGA
TCTGTTCGTTGCACAT---TCGCTCCCAGTCCAGATTCCTTCCTGATACTCGTCATGTGAGGA

";
        let reader = AxtReader::new(input.as_bytes());
        let axts: Vec<Axt> = reader.collect::<Result<Vec<_>, _>>().unwrap();

        assert_eq!(axts.len(), 1);
        let a = &axts[0];
        assert_eq!(a.id, 0);
        assert_eq!(a.t_name, "chr19");
        assert_eq!(a.t_start, 3001011); // 3001012 - 1
        assert_eq!(a.t_end, 3001075);
        assert_eq!(a.q_name, "chr11");
        assert_eq!(a.q_start, 70568379); // 70568380 - 1
        assert_eq!(a.q_end, 70568443);
        assert_eq!(a.q_strand, '-');
        assert_eq!(a.score, Some(3500));
        assert_eq!(a.t_sym.len(), 63);
        assert_eq!(a.q_sym.len(), 63);
    }
}
