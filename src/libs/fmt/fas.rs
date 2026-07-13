use intspan::Range;
use std::collections::VecDeque;
use std::io::Write;
use std::{fmt, io, str};

use crate::libs::io::LinesRef;

/// A single sequence entry in a block FA file, with its genomic range and aligned sequence.
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
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    /// Creates an empty FasEntry.
    pub fn new() -> Self {
        Self {
            range: Range::new(),
            seq: vec![],
        }
    }

    /// Constructed from range and seq
    ///
    /// ```ignore
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
/// ```ignore
/// # use intspan::Range;
/// # use pgr::libs::fmt::fas::FasEntry;
/// let range = Range::from("I", 1, 10);
/// let seq = "ACAGCTGA-AA".as_bytes().to_vec();
/// let entry = FasEntry::from(&range, &seq);
/// assert_eq!(entry.to_string(), ">I:1-10\nACAGCTGA-AA\n");
/// ```
impl fmt::Display for FasEntry {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let seq = str::from_utf8(self.seq()).map_err(|_| fmt::Error)?;
        write!(f, ">{}\n{}\n", self.range(), seq)?;
        Ok(())
    }
}

/// A block FA alignment block, containing one or more aligned sequence entries.
pub struct FasBlock {
    /// Aligned sequence entries in this block.
    pub entries: Vec<FasEntry>,
    /// Species/genome name for each entry.
    pub names: Vec<String>,
    /// Header strings (range descriptions) for each entry.
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
        header.ok_or(io::Error::new(io::ErrorKind::UnexpectedEof, "EOF"))?,
        LinesRef { buf: &mut input },
    )?;
    Ok(block)
}

/// Iterator over FasBlock records from a reader.
///
/// Wraps [`next_fas_block`], treating `UnexpectedEof` as the end of iteration
/// and propagating other errors as `anyhow::Error`.
pub struct FasBlockIter<'a, R: io::BufRead + ?Sized> {
    reader: &'a mut R,
}

impl<'a, R: io::BufRead + ?Sized> Iterator for FasBlockIter<'a, R> {
    type Item = anyhow::Result<FasBlock>;

    fn next(&mut self) -> Option<Self::Item> {
        match next_fas_block(self.reader) {
            Ok(block) => Some(Ok(block)),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => None,
            Err(e) => Some(Err(anyhow::Error::from(e))),
        }
    }
}

/// Create a FasBlock iterator from a reader.
pub fn iter_fas_blocks<R: io::BufRead + ?Sized>(reader: &mut R) -> FasBlockIter<'_, R> {
    FasBlockIter { reader }
}

/// Parse a single FasBlock from its header line and the following non-empty lines.
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
        if !h.starts_with('>') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Expected FAS header line starting with '>'",
            ));
        }
        let header = &h[1..];
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

        block_names.push(range.name().to_string());
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
    use std::sync::{Arc, Mutex};

    let (snd1, rcv1) = crossbeam::channel::bounded::<FasBlock>(10);
    let (snd2, rcv2) = crossbeam::channel::bounded::<String>(10);
    let errors: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));

    let write_result = crossbeam::scope(|s| {
        // Reader thread.
        let reader_errors = Arc::clone(&errors);
        s.spawn(move |_| {
            for infile in infiles {
                let mut reader = match crate::reader(infile) {
                    Ok(r) => r,
                    Err(e) => {
                        reader_errors
                            .lock()
                            .unwrap()
                            .push(format!("failed to open reader for {}: {}", infile, e));
                        continue;
                    }
                };
                for block_result in iter_fas_blocks(&mut reader) {
                    match block_result {
                        Ok(block) => {
                            if snd1.send(block).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            let _ = writeln!(
                                std::io::stderr(),
                                "pgr: warning: skipping malformed fas block: {}",
                                e
                            );
                        }
                    }
                }
            }
            drop(snd1);
        });

        // Worker threads.
        for _ in 0..parallel {
            let (sendr, recvr) = (snd2.clone(), rcv1.clone());
            let errors = Arc::clone(&errors);
            s.spawn(move |_| {
                for block in recvr.iter() {
                    match proc_block(&block) {
                        Ok(out_string) => {
                            if sendr.send(out_string).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            errors
                                .lock()
                                .unwrap()
                                .push(format!("fas block processing failed: {}", e));
                        }
                    }
                }
            });
        }
        drop(snd2);

        // Writer thread (runs on this thread).
        let mut result = Ok(());
        for out_string in rcv2.iter() {
            if let Err(e) = writer.write_all(out_string.as_ref()) {
                result = Err(e);
                break;
            }
        }
        result
    })
    .map_err(|_| anyhow::anyhow!("parallel pipeline failed (worker panic)"))?;

    if let Err(e) = write_result {
        return Err(anyhow::Error::from(e));
    }

    let errors = Arc::try_unwrap(errors)
        .map_err(|_| anyhow::anyhow!("errors Arc still shared"))?
        .into_inner()
        .map_err(|e| anyhow::anyhow!("errors mutex poisoned: {}", e))?;
    if !errors.is_empty() {
        anyhow::bail!(
            "{} block(s) failed during parallel processing:\n{}",
            errors.len(),
            errors.join("\n")
        );
    }

    Ok(())
}

/// Process FasBlock files either single-threaded or in parallel.
///
/// For each block read from `infiles`, calls `proc_block` to produce a string
/// chunk, and writes all chunks to `writer`. When `parallel > 1`, delegates to
/// [`run_parallel`] with `parallel` worker threads (output order may differ
/// from input order). Flushes `writer` before returning.
pub fn run_pipeline<W, F>(
    writer: &mut W,
    infiles: &[String],
    parallel: usize,
    proc_block: F,
) -> anyhow::Result<()>
where
    W: Write,
    F: Fn(&FasBlock) -> anyhow::Result<String> + Sync,
{
    if parallel <= 1 {
        for infile in infiles {
            let mut reader = crate::reader(infile)?;
            for block_result in iter_fas_blocks(&mut reader) {
                match block_result {
                    Ok(block) => {
                        let out_string = proc_block(&block)?;
                        writer.write_all(out_string.as_ref())?;
                    }
                    Err(e) => {
                        let _ = writeln!(
                            std::io::stderr(),
                            "pgr: warning: skipping malformed fas block: {}",
                            e
                        );
                    }
                }
            }
        }
    } else {
        run_parallel(infiles, parallel, writer, &proc_block)?;
    }

    writer.flush()?;
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
    for block_result in iter_fas_blocks(reader) {
        let block = block_result?;
        let first_entry = block
            .entries
            .first()
            .ok_or_else(|| anyhow::anyhow!("empty fas block encountered while concatenating"))?;
        let length = first_entry.seq().len();

        for name in needed {
            if let Some(idx) = block.names.iter().position(|n| n == name) {
                let seq = std::str::from_utf8(block.entries[idx].seq())?;
                seq_of.entry(name.to_string()).and_modify(|e| *e += seq);
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
    res_of: &mut std::collections::BTreeMap<
        String,
        std::collections::BTreeMap<String, intspan::IntSpan>,
    >,
    name_filter: &str,
    trim: i32,
) -> anyhow::Result<()> {
    for block_result in iter_fas_blocks(reader) {
        let block = block_result?;
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

        for (idx, entry) in block.entries.iter().enumerate() {
            let range = entry.range();
            if !range.is_valid() {
                continue;
            }

            let name = &block_names[idx];
            if !name_filter.is_empty() && name_filter != name {
                continue;
            }

            let res = res_of
                .get_mut(name)
                .ok_or_else(|| anyhow::anyhow!("name not found in res_of: {}", name))?;

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
pub fn find_best_pairs(entries: &[FasEntry]) -> anyhow::Result<Vec<(usize, usize)>> {
    let n = entries.len();
    if n < 2 {
        return Ok(vec![]);
    }
    let mut best_pair: Vec<(usize, usize)> = vec![];
    for i in 0..n {
        let mut dist_idx: Option<(f32, usize)> = None;
        for j in 0..n {
            if i == j {
                continue;
            }
            let dist = crate::libs::alignment::pair_d(entries[i].seq(), entries[j].seq())?;
            if dist_idx.map(|d| dist < d.0).unwrap_or(true) {
                dist_idx = Some((dist, j));
            }
        }
        if let Some((_, j)) = dist_idx {
            if i < j {
                best_pair.push((i, j));
            } else {
                best_pair.push((j, i));
            }
        }
    }
    // Deduplicate pairs preserving first-occurrence order
    let mut seen = std::collections::HashSet::new();
    Ok(best_pair.into_iter().filter(|p| seen.insert(*p)).collect())
}

/// Add entries from a block to the join map, keyed by the target entry's range.
pub fn join_block_entries(
    block: &FasBlock,
    name: &str,
    block_of: &mut std::collections::BTreeMap<String, Vec<FasEntry>>,
) -> anyhow::Result<()> {
    let idx = match block.names.iter().position(|x| x == name) {
        Some(i) => i,
        None => return Ok(()),
    };
    let header = block.entries[idx].range().to_string();

    if !block_of.contains_key(&header) {
        block_of.insert(header.clone(), vec![]);
        block_of
            .get_mut(&header)
            .ok_or_else(|| anyhow::anyhow!("inserted header missing"))?
            .push(block.entries[idx].clone());
    }

    for entry in &block.entries {
        if entry.range().name() != name {
            block_of
                .get_mut(&header)
                .ok_or_else(|| anyhow::anyhow!("header missing in block_of"))?
                .push(entry.clone());
        }
    }
    Ok(())
}

/// Concatenate FasEntry records into a single block string without a trailing newline.
fn block_to_string(entries: &[FasEntry]) -> String {
    let mut s = String::new();
    for entry in entries {
        s.push_str(&entry.to_string());
    }
    if s.ends_with('\n') {
        s.pop();
    }
    s
}

/// Generate output blocks (each a complete string) with header replacement applied.
pub fn replace_block_lines(
    block: &FasBlock,
    replace_of: &std::collections::BTreeMap<String, Vec<String>>,
) -> anyhow::Result<Vec<String>> {
    let matched: Vec<&String> = replace_of
        .keys()
        .filter(|e| block.headers.contains(*e))
        .collect();

    let mut blocks = Vec::new();

    if matched.len() != 1 {
        if matched.len() > 1 {
            log::warn!("Doesn't support replacing multiple records in one block");
        }
        blocks.push(block_to_string(&block.entries));
    } else {
        let original = matched[0];
        let occ = block.headers.iter().filter(|h| *h == original).count();
        if occ != 1 {
            log::warn!(
                "Header '{}' appears {} times in one block; keeping block unchanged",
                original,
                occ
            );
            blocks.push(block_to_string(&block.entries));
        } else {
            let idx = block
                .headers
                .iter()
                .position(|e| e == original)
                .ok_or_else(|| anyhow::anyhow!("matched header not found"))?;
            for new in &replace_of[original] {
                let mut s = String::new();
                for (i, entry) in block.entries.iter().enumerate() {
                    if i == idx {
                        s.push_str(&format!(
                            ">{}\n{}\n",
                            new,
                            String::from_utf8(entry.seq().to_vec())?
                        ));
                    } else {
                        s.push_str(&entry.to_string());
                    }
                }
                if s.ends_with('\n') {
                    s.pop();
                }
                blocks.push(s);
            }
        }
    }
    Ok(blocks)
}

/// Create block FA content from a links-of-ranges TSV reader. For each line,
/// splits on tab, parses each field as an intspan::Range, optionally overrides
/// the species name, fetches the sequence via `get_seq_loc`, and writes
/// `>{range}\n{seq}\n` per range with a blank line separating blocks.
pub fn create_from_links<R: io::BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    genome: &str,
    name: &str,
) -> anyhow::Result<()> {
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        let mut wrote_entry = false;
        for part in &parts {
            let mut range = Range::from_str(part);
            if !range.is_valid() {
                log::warn!("skipping invalid range: {}", part);
                continue;
            }
            if !name.is_empty() {
                *range.name_mut() = name.to_string();
            }
            let seq = crate::libs::loc::get_seq_loc(genome, &range.to_string())?;
            if seq.is_empty() {
                log::warn!("skipping range with no sequence: {}", range);
                continue;
            }
            writer.write_all(format!(">{}\n{}\n", range, seq).as_ref())?;
            wrote_entry = true;
        }
        if wrote_entry {
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

// ============================================================================
// consensus_block
// ============================================================================

/// Options for [`consensus_block`].
pub struct ConsensusOptions {
    /// Consensus sequence name written to the output header.
    pub cname: String,
    /// Whether the last entry of each block is an outgroup to be preserved.
    pub has_outgroup: bool,
    /// POA engine selector: `"builtin"` or `"spoa"`.
    pub engine: String,
    /// POA scoring parameters.
    pub params: crate::libs::poa::AlignmentParams,
    /// Alignment mode code: 0=local, 1=global, 2=semi_global.
    pub algo_code: i32,
}

/// Build consensus for one [`FasBlock`] and return a fas-formatted string.
pub fn consensus_block(block: &FasBlock, opts: &ConsensusOptions) -> anyhow::Result<String> {
    use std::fmt::Write;
    let outgroup = if opts.has_outgroup {
        block.entries.last()
    } else {
        None
    };

    let mut seqs: Vec<&[u8]> = Vec::with_capacity(block.entries.len());
    for entry in &block.entries {
        seqs.push(entry.seq());
    }
    if outgroup.is_some() {
        seqs.pop(); // Remove the outgroup sequence
    }

    // Generate consensus sequence
    let mut cons = match opts.engine.as_str() {
        "spoa" => crate::libs::alignment::get_consensus_poa_external(
            &seqs,
            opts.params.match_score,
            opts.params.mismatch_score,
            opts.params.gap_open,
            opts.params.gap_extend,
            opts.algo_code,
        )?,
        _ => crate::libs::alignment::get_consensus_poa_builtin(
            &seqs,
            opts.params.match_score,
            opts.params.mismatch_score,
            opts.params.gap_open,
            opts.params.gap_extend,
            opts.algo_code,
        )?,
    };
    cons = cons.replace('-', "");

    let mut range = match block.entries.first() {
        Some(e) => e.range().clone(),
        None => anyhow::bail!("empty block"),
    };

    let mut out_string = String::new();
    if range.is_valid() {
        *range.name_mut() = opts.cname.clone();
        writeln!(out_string, ">{}\n{}", range, cons)?;
    } else {
        writeln!(out_string, ">{}\n{}", opts.cname, cons)?;
    }
    if let Some(og) = outgroup {
        out_string.push_str(&og.to_string());
    }

    // end of a block
    out_string.push('\n');

    Ok(out_string)
}

// ============================================================================
// refine_block
// ============================================================================

/// Options for [`refine_block`].
pub struct RefineOptions<'a> {
    /// MSA engine selector: `"builtin"`, `"clustalw"`, `"mafft"`, `"muscle"`, `"spoa"`, or `"none"`.
    pub engine: &'a str,
    /// Whether the last entry of each block is an outgroup.
    pub has_outgroup: bool,
    /// Chop head and tail indels of this length (0 disables).
    pub chop: usize,
    /// Quick mode: only align indel-adjacent regions.
    pub is_quick: bool,
    /// In quick mode, enlarge indel regions by this padding.
    pub pad: usize,
    /// In quick mode, fill holes between indels up to this distance.
    pub fill: usize,
}

/// Realign and trim one [`FasBlock`], return a fas-formatted string.
pub fn refine_block(block: &FasBlock, opts: &RefineOptions) -> anyhow::Result<String> {
    use std::fmt::Write;

    let n = block.entries.len();
    let mut seqs: Vec<String> = Vec::with_capacity(n);
    let mut ranges = Vec::with_capacity(n);
    for entry in &block.entries {
        seqs.push(String::from_utf8(entry.seq().to_vec())?);
        ranges.push(entry.range().clone());
    }

    let mut aligned = vec![];
    if opts.engine == "none" {
        aligned = seqs;
    } else if opts.is_quick {
        let pad_i32 = i32::try_from(opts.pad)
            .map_err(|_| anyhow::anyhow!("--indel-pad {} exceeds i32 range", opts.pad))?;
        let fill_i32 = i32::try_from(opts.fill)
            .map_err(|_| anyhow::anyhow!("--fill {} exceeds i32 range", opts.fill))?;
        aligned = crate::libs::alignment::align_seqs_quick(&seqs, opts.engine, pad_i32, fill_i32)?;
    } else {
        aligned = crate::libs::alignment::align_seqs(&seqs, opts.engine)?;
    }

    crate::libs::alignment::trim_pure_dash(&mut aligned);
    if opts.has_outgroup {
        crate::libs::alignment::trim_outgroup(&mut aligned)?;
        crate::libs::alignment::trim_complex_indel(&mut aligned)?;
    }

    if opts.chop > 0 {
        crate::libs::alignment::trim_head_tail(&mut aligned, &mut ranges, opts.chop);
    }

    let mut out_string = String::new();
    for (range, seq) in ranges.iter().zip(aligned) {
        writeln!(out_string, ">{}\n{}", range, seq)?;
    }

    // end of a block
    out_string.push('\n');

    Ok(out_string)
}

#[cfg(test)]
mod fas_tests {
    use std::io::BufReader;

    #[test]
    fn parse_fas_block_range() {
        let str = ">S288c.I(+):13267-13287
TCGTCAGTTGGTTGACCATTA
>YJM789.gi_151941327(-):5668-5688
TCGTCAGTTGGTTGACCATTA
>RM11.gi_61385832(-):5590-5610
TCGTCAGTTGGTTGACCATTA
>Spar.gi_29362400(+):2477-2497
TCATCAGTTGGCAAACCGTTA

>S288c.I(+):185273-185334
GCATATAATATGAACCAATATCTA-TTCATGAAGAGACTATGGTATACCCGGTACTATTTCTA
>YJM789.gi_151941327(+):156665-156726
GCGTATAATATGAACCAGTATCTTTTTCATGAAG-GGCTATGGTATACTCCATATTACTTCTA
>RM11.gi_61385833(-):3668-3730
GCATATAATATGAACCAATATCTATTTCATGGAGAGACTATGATAT-CCCCGTACTATTTCTA
>Spar.gi_29362478(-):2102-2161
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

    #[test]
    fn parse_fas_block_rejects_missing_header() {
        let str = ">S288c.I(+):13267-13287
TCGTCAGTTGGTTGACCATTA
ACGT\n";
        let mut reader = BufReader::new(str.as_bytes());
        let result = crate::libs::fmt::fas::next_fas_block(&mut reader);
        assert!(
            result.is_err(),
            "non-header sequence line should be rejected"
        );
    }
}
