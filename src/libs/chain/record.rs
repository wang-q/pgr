use std::io::BufRead;
use std::str::FromStr;

/// Header information for a chain, containing score and sequence details.
#[derive(Debug, Clone, Default)]
pub struct ChainHeader {
    /// Chain score.
    pub score: f64,
    /// Target sequence name.
    pub t_name: String,
    /// Target sequence size.
    pub t_size: u64,
    /// Target strand ('+' or '-').
    pub t_strand: char,
    /// Target start coordinate (0-based).
    pub t_start: u64,
    /// Target end coordinate (0-based, exclusive).
    pub t_end: u64,
    /// Query sequence name.
    pub q_name: String,
    /// Query sequence size.
    pub q_size: u64,
    /// Query strand ('+' or '-').
    pub q_strand: char,
    /// Query start coordinate (0-based).
    pub q_start: u64,
    /// Query end coordinate (0-based, exclusive).
    pub q_end: u64,
    /// Chain ID.
    pub id: u64,
}

/// Data for a single block in a chain (size, and gap to next block).
#[derive(Debug, Clone, Default)]
pub struct ChainData {
    /// Size of the alignment block.
    pub size: u64,
    /// Gap in target sequence to the next block.
    pub dt: u64,
    /// Gap in query sequence to the next block.
    pub dq: u64,
}

/// A simplified representation of an alignment block with absolute coordinates.
#[derive(Debug, Clone, Default)]
pub struct Block {
    /// Target start (0-based).
    pub t_start: u64,
    /// Target end (0-based, exclusive).
    pub t_end: u64,
    /// Query start (0-based).
    pub q_start: u64,
    /// Query end (0-based, exclusive).
    pub q_end: u64,
}

/// Represents a complete chain with header and data blocks.
#[derive(Debug, Clone, Default)]
pub struct Chain {
    /// Chain header containing score, sequence names, coordinates, and ID.
    pub header: ChainHeader,
    /// Alignment blocks (size, dt, dq) constituting the chain.
    pub data: Vec<ChainData>,
}

impl Chain {
    /// Convert chain data (relative coordinates) to blocks (absolute coordinates).
    ///
    /// Note: Coordinates are 0-based, half-open [start, end).
    /// This handles strand logic:
    /// - t_strand is always '+'.
    /// - q_strand can be '-' (which means q_start/q_end are on the reverse strand coordinates,
    ///   but the Block struct stores them as increasing numbers on that strand).
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
    ///
    /// Assumes blocks are sorted by t_start, non-overlapping, and consistent
    /// with the header. Will update header.t_start, t_end, q_start, q_end
    /// based on the blocks. Overlapping blocks are treated as having zero gap.
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
                    next.t_start.saturating_sub(curr.t_end),
                    next.q_start.saturating_sub(curr.q_end),
                )
            } else {
                (0, 0)
            };

            data.push(ChainData { size, dt, dq });
        }
        data
    }

    /// Write the chain in UCSC Chain format.
    ///
    /// Header fields are space-separated; gap data lines are tab-separated
    /// (`size\tdt\tdq`), matching the UCSC convention. The score is formatted
    /// as an integer (`{:.0}`) since chain scores are conventionally whole
    /// numbers.
    pub fn write<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writeln!(
            writer,
            "chain {:.0} {} {} {} {} {} {} {} {} {} {} {}",
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
                writeln!(writer, "{}\t{}\t{}", d.size, d.dt, d.dq)?;
            }
        }
        writeln!(writer)?;
        Ok(())
    }

    /// Extract a subset of the chain overlapping with the given target range.
    ///
    /// Returns `None` if no overlap is found.
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
            t_strand: parts[4]
                .chars()
                .next()
                .ok_or_else(|| anyhow::anyhow!("empty t_strand field"))?,
            t_start: parts[5].parse()?,
            t_end: parts[6].parse()?,
            q_name: parts[7].to_string(),
            q_size: parts[8].parse()?,
            q_strand: parts[9]
                .chars()
                .next()
                .ok_or_else(|| anyhow::anyhow!("empty q_strand field"))?,
            q_start: parts[10].parse()?,
            q_end: parts[11].parse()?,
            id: parts[12].parse()?,
        })
    }
}

/// A buffered reader for UCSC Chain format files.
///
/// Non-chain, non-comment lines encountered while scanning for a header are silently ignored.
pub struct ChainReader<R> {
    reader: std::io::BufReader<R>,
    next_line: Option<String>,
    /// Header/comments lines (starting with `#`) collected before the first chain header.
    pub header_comments: Vec<String>,
}

impl<R: std::io::Read> ChainReader<R> {
    /// Creates a new `ChainReader` from any `Read` source.
    pub fn new(inner: R) -> Self {
        Self {
            reader: std::io::BufReader::new(inner),
            next_line: None,
            header_comments: Vec::new(),
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

/// Reads all chains from a reader into a vector.
pub fn read_chains<R: std::io::Read>(reader: R) -> anyhow::Result<Vec<Chain>> {
    let chain_reader = ChainReader::new(reader);
    chain_reader.collect()
}

impl<R: std::io::Read> Iterator for ChainReader<R> {
    type Item = anyhow::Result<Chain>;

    fn next(&mut self) -> Option<Self::Item> {
        // Find next chain header
        let header_line = loop {
            match self.read_line() {
                Ok(Some(line)) => {
                    let trimmed = line.trim();
                    if trimmed.is_empty() {
                        continue;
                    }
                    if trimmed.starts_with('#') {
                        self.header_comments.push(line);
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
                        match parts[0].parse() {
                            Ok(size) => data.push(ChainData { size, dt: 0, dq: 0 }),
                            Err(e) => {
                                return Some(Err(anyhow::anyhow!(
                                    "Invalid chain data line '{}': failed to parse size: {}",
                                    trimmed,
                                    e
                                )))
                            }
                        }
                    } else if parts.len() == 3 {
                        // size dt dq
                        let size = match parts[0].parse() {
                            Ok(v) => v,
                            Err(e) => {
                                return Some(Err(anyhow::anyhow!(
                                    "Invalid chain data line '{}': failed to parse size: {}",
                                    trimmed,
                                    e
                                )))
                            }
                        };
                        let dt = match parts[1].parse() {
                            Ok(v) => v,
                            Err(e) => {
                                return Some(Err(anyhow::anyhow!(
                                    "Invalid chain data line '{}': failed to parse dt: {}",
                                    trimmed,
                                    e
                                )))
                            }
                        };
                        let dq = match parts[2].parse() {
                            Ok(v) => v,
                            Err(e) => {
                                return Some(Err(anyhow::anyhow!(
                                    "Invalid chain data line '{}': failed to parse dq: {}",
                                    trimmed,
                                    e
                                )))
                            }
                        };
                        data.push(ChainData { size, dt, dq });
                    } else {
                        return Some(Err(anyhow::anyhow!("Invalid chain data line: {}", trimmed)));
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

    #[test]
    fn test_chain_block_conversion() {
        let mut header = ChainHeader {
            t_start: 100,
            q_start: 200,
            ..Default::default()
        };
        let data = vec![
            ChainData {
                size: 10,
                dt: 5,
                dq: 5,
            },
            ChainData {
                size: 20,
                dt: 0,
                dq: 0,
            },
        ];
        let chain = Chain {
            header: header.clone(),
            data: data.clone(),
        };

        // To blocks
        let blocks = chain.to_blocks();
        assert_eq!(blocks.len(), 2);

        // Block 1: start at 100/200, size 10
        assert_eq!(blocks[0].t_start, 100);
        assert_eq!(blocks[0].t_end, 110);
        assert_eq!(blocks[0].q_start, 200);
        assert_eq!(blocks[0].q_end, 210);

        // Gap: 5, 5. Next start: 110+5=115, 210+5=215
        // Block 2: start at 115/215, size 20
        assert_eq!(blocks[1].t_start, 115);
        assert_eq!(blocks[1].t_end, 135);
        assert_eq!(blocks[1].q_start, 215);
        assert_eq!(blocks[1].q_end, 235);

        // From blocks
        let new_data = Chain::from_blocks(&mut header, &blocks);
        assert_eq!(new_data.len(), 2);
        assert_eq!(new_data[0].size, 10);
        assert_eq!(new_data[0].dt, 5);
        assert_eq!(new_data[0].dq, 5);
        assert_eq!(new_data[1].size, 20);
        assert_eq!(new_data[1].dt, 0);
        assert_eq!(new_data[1].dq, 0);
    }

    #[test]
    fn test_iterator_rejects_malformed_data_line() {
        let input = "\
chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
16 0 4
malformed
";
        let reader = ChainReader::new(input.as_bytes());
        let result: Result<Vec<_>, _> = reader.collect();
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("Invalid chain data line"));
    }

    #[test]
    fn test_from_blocks_overlapping_does_not_panic() {
        let mut header = ChainHeader::default();
        let blocks = vec![
            Block {
                t_start: 10,
                t_end: 20,
                q_start: 100,
                q_end: 110,
            },
            Block {
                t_start: 15,
                t_end: 25,
                q_start: 105,
                q_end: 115,
            },
        ];
        let data = Chain::from_blocks(&mut header, &blocks);
        assert_eq!(data.len(), 2);
        // Overlap produces zero gap via saturating_sub.
        assert_eq!(data[0].dt, 0);
        assert_eq!(data[0].dq, 0);
    }

    #[test]
    fn test_subset_reverse_strand() {
        let chain = Chain {
            header: ChainHeader {
                score: 100.0,
                t_name: "chrT".to_string(),
                t_size: 1000,
                t_strand: '+',
                t_start: 100,
                t_end: 350,
                q_name: "chrQ".to_string(),
                q_size: 1000,
                q_strand: '-',
                q_start: 600,
                q_end: 850,
                id: 1,
            },
            data: vec![
                ChainData {
                    size: 100,
                    dt: 50,
                    dq: 50,
                },
                ChainData {
                    size: 100,
                    dt: 0,
                    dq: 0,
                },
            ],
        };

        let sub = chain.subset(120, 280).expect("subset should overlap");

        assert_eq!(sub.header.t_strand, '+');
        assert_eq!(sub.header.q_strand, '-');
        assert_eq!(sub.header.t_start, 120);
        assert_eq!(sub.header.t_end, 280);
        // Query coordinates remain in reverse-strand space.
        assert_eq!(sub.header.q_start, 620);
        assert_eq!(sub.header.q_end, 780);

        assert_eq!(sub.data.len(), 2);
        assert_eq!(sub.data[0].size, 80);
        assert_eq!(sub.data[0].dt, 50);
        assert_eq!(sub.data[0].dq, 50);
        assert_eq!(sub.data[1].size, 30);
        assert_eq!(sub.data[1].dt, 0);
        assert_eq!(sub.data[1].dq, 0);
    }
}
