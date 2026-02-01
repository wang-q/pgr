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
pub struct Block {
    pub t_start: u64,
    pub t_end: u64,
    pub q_start: u64,
    pub q_end: u64,
}

#[derive(Debug, Clone, Default)]
pub struct Chain {
    pub header: ChainHeader,
    pub data: Vec<ChainData>,
}

impl Chain {
    /// Convert chain data (relative coordinates) to blocks (absolute coordinates).
    /// Note: Coordinates are 0-based, half-open [start, end).
    /// This handles strand logic:
    /// - t_strand is always '+'.
    /// - q_strand can be '-' (which means q_start/q_end are on the reverse strand coordinates,
    ///   but the Block struct stores them as increasing numbers on that strand).
    ///   Wait, in UCSC chain format, if qStrand is '-', qStart and qEnd are already in reverse strand coordinates
    ///   relative to the qSize. Specifically, qStart = qSize - real_end, qEnd = qSize - real_start.
    ///   However, the chain format defines the alignment blocks.
    ///   Let's verify the logic:
    ///   t_current starts at header.t_start
    ///   q_current starts at header.q_start
    ///   For each line:
    ///     block size: both advance by size
    ///     dt: t advances by dt
    ///     dq: q advances by dq
    pub fn to_blocks(&self) -> Vec<Block> {
        let mut blocks = Vec::with_capacity(self.data.len());
        let mut t_curr = self.header.t_start;
        let mut q_curr = self.header.q_start;

        for d in &self.data {
            blocks.push(Block {
                t_start: t_curr,
                t_end: t_curr + d.size,
                q_start: q_curr,
                q_end: q_curr + d.size,
            });

            t_curr += d.size + d.dt;
            q_curr += d.size + d.dq;
        }

        blocks
    }

    /// Reconstruct chain data from blocks.
    /// Assumes blocks are sorted by t_start and consistent with the header.
    /// Will update header.t_start, t_end, q_start, q_end based on the blocks.
    pub fn from_blocks(header: &mut ChainHeader, blocks: &[Block]) -> Vec<ChainData> {
        if blocks.is_empty() {
            return Vec::new();
        }

        // Update header range
        header.t_start = blocks.first().unwrap().t_start;
        header.t_end = blocks.last().unwrap().t_end;
        header.q_start = blocks.first().unwrap().q_start;
        header.q_end = blocks.last().unwrap().q_end;

        let mut data = Vec::with_capacity(blocks.len());
        for i in 0..blocks.len() {
            let curr = &blocks[i];
            let size = curr.t_end - curr.t_start;
            
            // Sanity check
            // assert_eq!(size, curr.q_end - curr.q_start);

            let (dt, dq) = if i < blocks.len() - 1 {
                let next = &blocks[i + 1];
                (
                    next.t_start - curr.t_end,
                    next.q_start - curr.q_end
                )
            } else {
                (0, 0)
            };

            data.push(ChainData {
                size,
                dt,
                dq,
            });
        }
        data
    }

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

    pub fn subset(&self, t_start: u64, t_end: u64) -> Option<Chain> {
        let blocks = self.to_blocks();
        let mut new_blocks = Vec::new();

        for b in blocks {
            // Check for overlap
            let start = std::cmp::max(b.t_start, t_start);
            let end = std::cmp::min(b.t_end, t_end);

            if start < end {
                let offset = start - b.t_start;
                let len = end - start;
                new_blocks.push(Block {
                    t_start: start,
                    t_end: end,
                    q_start: b.q_start + offset,
                    q_end: b.q_start + offset + len,
                });
            }
        }

        if new_blocks.is_empty() {
            None
        } else {
            let mut header = self.header.clone();
            let data = Chain::from_blocks(&mut header, &new_blocks);
            Some(Chain { header, data })
        }
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

pub fn read_chains<R: std::io::Read>(reader: R) -> anyhow::Result<Vec<Chain>> {
    let mut reader = ChainReader::new(reader);
    let mut chains = Vec::new();

    while let Some(line) = reader.read_line()? {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with("chain") {
            let header = ChainHeader::from_str(line)?;
            let mut data = Vec::new();

            while let Some(inner_line) = reader.read_line()? {
                let inner_line_trim = inner_line.trim();
                if inner_line_trim.is_empty() {
                    break;
                }
                // Check if next line is a new chain (shouldn't happen if properly formatted with blank lines, but just in case)
                if inner_line_trim.starts_with("chain") {
                    reader.push_back(inner_line);
                    break;
                }

                let parts: Vec<&str> = inner_line_trim.split_whitespace().collect();
                if parts.len() == 1 {
                    data.push(ChainData {
                        size: parts[0].parse()?,
                        dt: 0,
                        dq: 0,
                    });
                } else if parts.len() == 3 {
                    data.push(ChainData {
                        size: parts[0].parse()?,
                        dt: parts[1].parse()?,
                        dq: parts[2].parse()?,
                    });
                } else {
                    return Err(anyhow::anyhow!("Invalid chain data line: {}", inner_line_trim));
                }
            }
            chains.push(Chain { header, data });
        }
    }

    Ok(chains)
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
