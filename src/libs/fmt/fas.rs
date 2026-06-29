use intspan::Range;
use std::collections::VecDeque;
use std::{fmt, io, str};

use crate::libs::io::LinesRef;

#[derive(Default, Clone)]
pub struct FasEntry {
    range: Range,
    seq: Vec<u8>,
}

impl FasEntry {
    // Immutable accessors
    pub fn range(&self) -> &Range {
        &self.range
    }
    pub fn seq(&self) -> &Vec<u8> {
        &self.seq
    }

    pub fn new() -> Self {
        Self {
            range: Range::new(),
            seq: vec![],
        }
    }

    /// Constructed from range and seq
    ///
    /// ```
    /// # use intspan::Range;
    /// # use pgr::libs::fmt::fas::FasEntry;
    /// let range = Range::from("I", 1, 10);
    /// let seq = "ACAGCTGA-AA".as_bytes().to_vec();
    /// let entry = FasEntry::from(&range, &seq);
    /// # assert_eq!(*entry.range().chr(), "I");
    /// # assert_eq!(*entry.range().start(), 1);
    /// # assert_eq!(*entry.range().end(), 10);
    /// # assert_eq!(std::str::from_utf8(entry.seq()).unwrap(), "ACAGCTGA-AA".to_string());
    /// ```
    pub fn from(range: &Range, seq: &[u8]) -> Self {
        Self {
            range: range.clone(),
            seq: seq.to_owned(),
        }
    }
}

/// To string
///
/// ```
/// # use intspan::Range;
/// # use pgr::libs::fmt::fas::FasEntry;
/// let range = Range::from("I", 1, 10);
/// let seq = "ACAGCTGA-AA".as_bytes().to_vec();
/// let entry = FasEntry::from(&range, &seq);
/// assert_eq!(entry.to_string(), ">I:1-10\nACAGCTGA-AA\n");
/// ```
impl fmt::Display for FasEntry {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            ">{}\n{}\n",
            self.range(),
            str::from_utf8(self.seq()).unwrap()
        )?;
        Ok(())
    }
}

/// A Fas alignment block.
pub struct FasBlock {
    pub entries: Vec<FasEntry>,
    pub names: Vec<String>,
    pub headers: Vec<String>,
}

/// Get the next FasBlock out of the input.
pub fn next_fas_block<T: io::BufRead + ?Sized>(mut input: &mut T) -> Result<FasBlock, io::Error> {
    let mut header: Option<String> = None;
    {
        let lines = LinesRef { buf: &mut input };
        for line_res in lines {
            let line: String = line_res?;
            if line.trim().is_empty() {
                // Blank line
                continue;
            }
            if line.starts_with('#') {
                // Fas comment
                continue;
            } else if line.starts_with('>') {
                // Start of a block
                header = Some(line);
                break;
            } else {
                // Shouldn't see this.
                return Err(io::Error::other("Unexpected line"));
            }
        }
    }
    let block = parse_fas_block(
        header.ok_or(io::Error::other("EOF"))?,
        LinesRef { buf: &mut input },
    )?;
    Ok(block)
}

pub fn parse_fas_block(
    header: String,
    iter: impl Iterator<Item = Result<String, io::Error>>,
) -> Result<FasBlock, io::Error> {
    let mut block_lines: VecDeque<String> = VecDeque::new();
    block_lines.push_back(header);

    for line_res in iter {
        let line: String = line_res?;
        if line.is_empty() {
            // Blank lines terminate the "paragraph".
            break;
        }
        block_lines.push_back(line);
    }
    let mut block_entries: Vec<FasEntry> = vec![];
    let mut block_names: Vec<String> = vec![];
    let mut block_headers: Vec<String> = vec![];

    while let Some(h) = block_lines.pop_front() {
        let header = match h.starts_with('>') {
            true => &h[1..],
            false => h.as_str(),
        };
        let range = Range::from_str(header);
        let seq = block_lines
            .pop_front()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "FAS block missing sequence line",
                )
            })?
            .as_bytes()
            .to_vec();

        let entry = FasEntry::from(&range, &seq);
        block_entries.push(entry);

        let name = if let Some(idx) = header.find("|species=") {
            let species = &header[idx + "|species=".len()..];
            let species = species
                .split(['|', ' ', '\t'])
                .next()
                .unwrap_or("")
                .to_string();
            species
        } else {
            range.name().to_string()
        };
        block_names.push(name);
        block_headers.push(header.to_string());
    }

    Ok(FasBlock {
        entries: block_entries,
        names: block_names,
        headers: block_headers,
    })
}

#[cfg(test)]
mod fas_tests {
    use std::io::BufReader;

    #[test]
    fn parse_fas_block_range() {
        let str = ">S288c.I(+):13267-13287|species=S288c
TCGTCAGTTGGTTGACCATTA
>YJM789.gi_151941327(-):5668-5688|species=YJM789
TCGTCAGTTGGTTGACCATTA
>RM11.gi_61385832(-):5590-5610|species=RM11
TCGTCAGTTGGTTGACCATTA
>Spar.gi_29362400(+):2477-2497|species=Spar
TCATCAGTTGGCAAACCGTTA

>S288c.I(+):185273-185334|species=S288c
GCATATAATATGAACCAATATCTA-TTCATGAAGAGACTATGGTATACCCGGTACTATTTCTA
>YJM789.gi_151941327(+):156665-156726|species=YJM789
GCGTATAATATGAACCAGTATCTTTTTCATGAAG-GGCTATGGTATACTCCATATTACTTCTA
>RM11.gi_61385833(-):3668-3730|species=RM11
GCATATAATATGAACCAATATCTATTTCATGGAGAGACTATGATAT-CCCCGTACTATTTCTA
>Spar.gi_29362478(-):2102-2161|species=Spar
GC-TAAAATATGAA-CGATATTTA-CCTGTAGAGGGACTATGGGAT-CCCCATACTACTTT--
";
        let mut reader = BufReader::new(str.as_bytes());

        let block = crate::libs::fmt::fas::next_fas_block(&mut reader).unwrap();
        assert_eq!(
            block.entries.first().unwrap().range.to_string(),
            "S288c.I(+):13267-13287".to_string()
        );
        assert_eq!(
            block.entries.get(2).unwrap().range.to_string(),
            "RM11.gi_61385832(-):5590-5610".to_string()
        );

        let block = crate::libs::fmt::fas::next_fas_block(&mut reader).unwrap();
        assert_eq!(
            String::from_utf8(block.entries.get(1).unwrap().seq.clone()).unwrap(),
            "GCGTATAATATGAACCAGTATCTTTTTCATGAAG-GGCTATGGTATACTCCATATTACTTCTA".to_string()
        );
    }
}
