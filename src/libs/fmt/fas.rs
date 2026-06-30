use intspan::Range;
use std::collections::VecDeque;
use std::io::Write;
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

/// Crossbeam parallel pipeline: 1 reader → N workers → 1 writer.
///
/// Reads FasBlock records from each path in `infiles`, calls `proc_block` on
/// each (in parallel across `parallel` workers), and writes the resulting
/// string chunks to `writer`. Output order may differ from input order.
pub fn run_parallel<W, F>(
    infiles: &[String],
    parallel: usize,
    writer: &mut W,
    proc_block: &F,
) -> anyhow::Result<()>
where
    W: Write,
    F: Fn(&FasBlock) -> anyhow::Result<String> + Sync,
{
    let (snd1, rcv1) = crossbeam::channel::bounded::<FasBlock>(10);
    let (snd2, rcv2) = crossbeam::channel::bounded::<String>(10);

    crossbeam::scope(|s| {
        // Reader thread.
        s.spawn(|_| {
            for infile in infiles {
                let mut reader = match crate::reader(infile) {
                    Ok(r) => r,
                    Err(_) => break,
                };
                while let Ok(block) = next_fas_block(&mut reader) {
                    if snd1.send(block).is_err() {
                        break;
                    }
                }
            }
            drop(snd1);
        });

        // Worker threads.
        for _ in 0..parallel {
            let (sendr, recvr) = (snd2.clone(), rcv1.clone());
            s.spawn(move |_| {
                for block in recvr.iter() {
                    if let Ok(out_string) = proc_block(&block) {
                        if sendr.send(out_string).is_err() {
                            break;
                        }
                    }
                }
            });
        }
        drop(snd2);

        // Writer thread (runs on this thread).
        for out_string in rcv2.iter() {
            if writer.write_all(out_string.as_ref()).is_err() {
                break;
            }
        }
    })
    .map_err(|_| anyhow::anyhow!("parallel pipeline failed (worker panic)"))?;

    Ok(())
}

/// Check a FasEntry's sequence against the reference genome.
pub fn check_entry_against_ref(
    entry: &FasEntry,
    reader: &mut crate::libs::loc::Input,
    loc_of: &indexmap::IndexMap<String, (u64, usize)>,
) -> anyhow::Result<String> {
    let range = entry.range();
    let seq = entry.seq().to_vec();
    let seq = std::str::from_utf8(&seq)?
        .to_string()
        .to_ascii_uppercase()
        .replace('-', "");

    let gseq = if loc_of.contains_key(range.chr()) {
        crate::libs::loc::fetch_range_seq(reader, loc_of, range)?.to_ascii_uppercase()
    } else {
        String::new()
    };

    let status = if seq == gseq { "OK" } else { "FAILED" };
    Ok(status.to_string())
}

/// Process fas blocks from reader, concatenating sequences for needed names.
pub fn concat_blocks_into<R: io::BufRead>(
    reader: &mut R,
    needed: &[String],
    seq_of: &mut std::collections::BTreeMap<String, String>,
) -> anyhow::Result<()> {
    while let Ok(block) = next_fas_block(reader) {
        let block_names = block.names;
        let length = block.entries.first().unwrap().seq().len();

        for name in needed {
            if block_names.contains(name) {
                for entry in &block.entries {
                    let entry_name = entry.range().name();
                    if entry_name == name {
                        let seq = std::str::from_utf8(entry.seq())?;
                        seq_of.entry(name.to_string()).and_modify(|e| *e += seq);
                    }
                }
            } else {
                seq_of
                    .entry(name.to_string())
                    .and_modify(|e| *e += "-".repeat(length).as_str());
            }
        }
    }
    Ok(())
}

/// Process fas blocks from reader, aggregating coverage into res_of.
pub fn aggregate_coverage_into<R: io::BufRead>(
    reader: &mut R,
    res_of: &mut std::collections::BTreeMap<String, std::collections::BTreeMap<String, intspan::IntSpan>>,
    name_filter: &str,
    trim: i32,
) -> anyhow::Result<()> {
    while let Ok(block) = next_fas_block(reader) {
        let block_names = block.names;

        if !name_filter.is_empty() {
            if !res_of.contains_key(name_filter) {
                res_of.insert(name_filter.to_string(), std::collections::BTreeMap::new());
            }
        } else {
            for name in &block_names {
                if !res_of.contains_key(name) {
                    res_of.insert(name.to_string(), std::collections::BTreeMap::new());
                }
            }
        }

        for entry in &block.entries {
            let range = entry.range();
            if !range.is_valid() {
                continue;
            }

            if !name_filter.is_empty() && name_filter != range.name() {
                continue;
            }

            let res = res_of.get_mut(range.name()).unwrap();

            if !res.contains_key(range.chr()) {
                res.insert(range.chr().to_string(), intspan::IntSpan::new());
            }

            let intspan = range.intspan().clone().trim(trim);
            res.get_mut(range.chr()).unwrap().merge(&intspan);
        }
    }
    Ok(())
}

/// Find best-to-best bilateral pairs based on sequence distance.
pub fn find_best_pairs(entries: &[FasEntry]) -> Vec<(usize, usize)> {
    let n = entries.len();
    let mut best_pair: Vec<(usize, usize)> = vec![];
    for i in 0..n {
        let mut dist_idx: (f32, usize) = (1.0, n - 1);
        for j in 0..n {
            if i == j {
                continue;
            }
            let dist = crate::libs::alignment::pair_d(entries[i].seq(), entries[j].seq());
            if dist < dist_idx.0 {
                dist_idx = (dist, j);
            }
        }
        if i < dist_idx.1 {
            best_pair.push((i, dist_idx.1));
        } else {
            best_pair.push((dist_idx.1, i));
        }
    }
    // Deduplicate pairs preserving first-occurrence order
    let mut seen = std::collections::HashSet::new();
    best_pair.into_iter().filter(|p| seen.insert(*p)).collect()
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
