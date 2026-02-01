use std::io::BufRead;
use std::str::FromStr;

#[derive(Debug, Clone, Default)]
pub struct ChainHeader {
    pub score: f64,
    pub t_name: String,
    pub t_size: u64,
    pub t_strand: char,
    pub t_start: u64,
    pub t_end: u64,
    pub q_name: String,
    pub q_size: u64,
    pub q_strand: char,
    pub q_start: u64,
    pub q_end: u64,
    pub id: u64,
}

#[derive(Debug, Clone, Default)]
pub struct ChainData {
    pub size: u64,
    pub dt: u64,
    pub dq: u64,
}

#[derive(Debug, Clone, Default)]
pub struct Chain {
    pub header: ChainHeader,
    pub data: Vec<ChainData>,
}

impl Chain {
    pub fn write<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writeln!(writer, "chain {} {} {} {} {} {} {} {} {} {} {} {}",
            self.header.score,
            self.header.t_name,
            self.header.t_size,
            self.header.t_strand,
            self.header.t_start,
            self.header.t_end,
            self.header.q_name,
            self.header.q_size,
            self.header.q_strand,
            self.header.q_start,
            self.header.q_end,
            self.header.id
        )?;

        let len = self.data.len();
        for (i, d) in self.data.iter().enumerate() {
            if i == len - 1 {
                writeln!(writer, "{}", d.size)?;
            } else {
                writeln!(writer, "{} {} {}", d.size, d.dt, d.dq)?;
            }
        }
        writeln!(writer)?;
        Ok(())
    }
}

impl FromStr for ChainHeader {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s.split_whitespace().collect();
        if parts.len() < 13 || parts[0] != "chain" {
            return Err(anyhow::anyhow!("Invalid chain header line"));
        }

        Ok(ChainHeader {
            score: parts[1].parse()?,
            t_name: parts[2].to_string(),
            t_size: parts[3].parse()?,
            t_strand: parts[4].chars().next().unwrap(),
            t_start: parts[5].parse()?,
            t_end: parts[6].parse()?,
            q_name: parts[7].to_string(),
            q_size: parts[8].parse()?,
            q_strand: parts[9].chars().next().unwrap(),
            q_start: parts[10].parse()?,
            q_end: parts[11].parse()?,
            id: parts[12].parse()?,
        })
    }
}

pub struct ChainReader<R> {
    reader: std::io::BufReader<R>,
    next_line: Option<String>,
}

impl<R: std::io::Read> ChainReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            reader: std::io::BufReader::new(inner),
            next_line: None,
        }
    }

    fn read_line(&mut self) -> std::io::Result<Option<String>> {
        if let Some(line) = self.next_line.take() {
            return Ok(Some(line));
        }
        let mut buf = String::new();
        let n = self.reader.read_line(&mut buf)?;
        if n == 0 {
            Ok(None)
        } else {
            Ok(Some(buf))
        }
    }

    fn push_back(&mut self, line: String) {
        self.next_line = Some(line);
    }
}

impl<R: std::io::Read> Iterator for ChainReader<R> {
    type Item = anyhow::Result<Chain>;

    fn next(&mut self) -> Option<Self::Item> {
        // Find next chain header
        let header_line = loop {
            match self.read_line() {
                Ok(Some(line)) => {
                    let trimmed = line.trim();
                    if trimmed.is_empty() || trimmed.starts_with('#') {
                        continue;
                    }
                    if trimmed.starts_with("chain") {
                        break trimmed.to_string();
                    }
                    // If we find non-chain data while looking for header, it's garbage or error, ignore
                }
                Ok(None) => return None, // EOF
                Err(e) => return Some(Err(anyhow::Error::new(e))),
            }
        };

        let header = match ChainHeader::from_str(&header_line) {
            Ok(h) => h,
            Err(e) => return Some(Err(e)),
        };

        let mut data = Vec::new();

        // Read data lines
        loop {
            match self.read_line() {
                Ok(Some(line)) => {
                    let trimmed = line.trim();
                    if trimmed.is_empty() {
                        continue; 
                    }
                    if trimmed.starts_with("chain") {
                        self.push_back(line); // Push back for next iteration
                        break;
                    }

                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() == 1 {
                        // Last line of block: size
                        if let Ok(size) = parts[0].parse() {
                            data.push(ChainData { size, dt: 0, dq: 0 });
                        }
                    } else if parts.len() == 3 {
                        // size dt dq
                        if let (Ok(size), Ok(dt), Ok(dq)) = (
                            parts[0].parse(),
                            parts[1].parse(),
                            parts[2].parse(),
                        ) {
                            data.push(ChainData { size, dt, dq });
                        }
                    } else {
                        // Invalid data line
                    }
                }
                Ok(None) => break, // EOF ends the current chain
                Err(e) => return Some(Err(anyhow::Error::new(e))),
            }
        }

        Some(Ok(Chain { header, data }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_chain() {
        let input = "\
chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
16 0 4
60 4 0
10
";
        let reader = ChainReader::new(input.as_bytes());
        let chains: Vec<Chain> = reader.collect::<Result<Vec<_>, _>>().unwrap();

        assert_eq!(chains.len(), 1);
        let c = &chains[0];
        assert_eq!(c.header.score, 4900.0);
        assert_eq!(c.header.t_name, "chrY");
        assert_eq!(c.header.t_strand, '+');
        assert_eq!(c.header.q_strand, '-');
        assert_eq!(c.data.len(), 3);
        assert_eq!(c.data[0].size, 16);
        assert_eq!(c.data[0].dt, 0);
        assert_eq!(c.data[0].dq, 4);
        assert_eq!(c.data[2].size, 10);
        assert_eq!(c.data[2].dt, 0);
        assert_eq!(c.data[2].dq, 0);
    }
}
