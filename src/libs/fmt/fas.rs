use anyhow::Context;
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

    /// Creates an entry from a range and sequence.
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
                let preview: String = line.chars().take(80).collect();
                return Err(io::Error::other(format!(
                    "Unexpected line in block FA: {}",
                    preview
                )));
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
        if line.trim().is_empty() {
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
    let seq = std::str::from_utf8(entry.seq())?
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
            res_of.entry(name_filter.to_string()).or_default();
        } else {
            for name in &block_names {
                res_of.entry(name.clone()).or_default();
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
            let intspan = range.intspan().clone().trim(trim);
            res.entry(range.chr().to_string())
                .or_default()
                .merge(&intspan);
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

    let entries = block_of.entry(header).or_default();
    if entries.is_empty() {
        entries.push(block.entries[idx].clone());
    }

    for entry in &block.entries {
        if entry.range().name() != name {
            entries.push(entry.clone());
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

    if matched.len() != 1 {
        if matched.len() > 1 {
            log::warn!("Doesn't support replacing multiple records in one block");
        }
        return Ok(vec![block_to_string(&block.entries)]);
    }

    let original = matched[0];
    let occ = block.headers.iter().filter(|h| *h == original).count();
    if occ != 1 {
        log::warn!(
            "Header '{}' appears {} times in one block; keeping block unchanged",
            original,
            occ
        );
        return Ok(vec![block_to_string(&block.entries)]);
    }

    let idx = block
        .headers
        .iter()
        .position(|e| e == original)
        .ok_or_else(|| anyhow::anyhow!("matched header not found"))?;

    let mut blocks = Vec::with_capacity(replace_of[original].len());
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
    Ok(blocks)
}

/// Format a sequence byte slice according to `--upper` and `--dash` flags.
///
/// When `is_dash` is true, gap characters (`-`) are removed. When `is_upper`
/// is true, the result is converted to ASCII uppercase.
pub fn format_sequence(seq: &[u8], is_dash: bool, is_upper: bool) -> String {
    let mut out = String::with_capacity(seq.len());
    for &nt in seq {
        if is_dash && nt == b'-' {
            continue;
        }
        let c = char::from(nt);
        out.push(if is_upper { c.to_ascii_uppercase() } else { c });
    }
    out
}

/// Filter and format one FasBlock.
///
/// Returns `Ok(None)` when the block should be skipped (missing species or
/// length out of range). Otherwise returns the formatted block string.
pub fn filter_block(
    block: &FasBlock,
    opt_name: &str,
    opt_min: Option<usize>,
    opt_max: Option<usize>,
    is_upper: bool,
    is_dash: bool,
) -> anyhow::Result<Option<String>> {
    if block.entries.is_empty() {
        return Ok(None);
    }

    let idx = if !opt_name.is_empty() {
        match block.names.iter().position(|x| x == opt_name) {
            Some(i) => i,
            None => return Ok(None),
        }
    } else {
        0
    };

    let idx_seq = block.entries[idx].seq();
    if let Some(min) = opt_min {
        if idx_seq.len() < min {
            return Ok(None);
        }
    }
    if let Some(max) = opt_max {
        if idx_seq.len() > max {
            return Ok(None);
        }
    }

    let mut out = String::new();
    for entry in &block.entries {
        let formatted = format_sequence(entry.seq(), is_dash, is_upper);
        let out_entry = FasEntry::from(entry.range(), formatted.as_bytes());
        out.push_str(&out_entry.to_string());
    }
    out.push('\n');
    Ok(Some(out))
}

/// Statistics for one FasBlock.
#[derive(Debug)]
pub struct BlockStat {
    pub target: String,
    pub length: usize,
    pub comparable: i32,
    pub difference: i32,
    pub gap: i32,
    pub ambiguous: i32,
    pub mean_d: f32,
    pub indel_span: i32,
}

impl BlockStat {
    /// Format the statistic as a TSV line without a trailing newline.
    pub fn to_tsv(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.target,
            self.length,
            self.comparable,
            self.difference,
            self.gap,
            self.ambiguous,
            self.mean_d,
            self.indel_span,
        )
    }
}

/// Compute statistics for a FasBlock.
///
/// When `has_outgroup` is true, the last entry is excluded from all
/// calculations except `length`.
pub fn compute_block_stat(block: &FasBlock, has_outgroup: bool) -> anyhow::Result<BlockStat> {
    if block.entries.is_empty() {
        anyhow::bail!("empty fas block");
    }

    let target = block.entries[0].range().to_string();
    let full_length = block.entries[0].seq().len();

    let mut seqs: Vec<&[u8]> = block.entries.iter().map(|e| e.seq()).collect();
    if has_outgroup {
        if seqs.len() < 2 {
            anyhow::bail!(
                "block has only {} entries, cannot apply --outgroup",
                seqs.len()
            );
        }
        seqs.pop();
    }

    let (_, comparable, difference, gap, ambiguous, mean_d) =
        crate::libs::alignment::alignment_stat(&seqs)?;

    let mut indel_ints = intspan::IntSpan::new();
    for seq in &seqs {
        indel_ints.merge(&crate::libs::alignment::indel_intspan(seq));
    }

    Ok(BlockStat {
        target,
        length: full_length,
        comparable,
        difference,
        gap,
        ambiguous,
        mean_d,
        indel_span: indel_ints.span_size() as i32,
    })
}

/// Write variations (substitutions) from a FasBlock to a writer.
///
/// `has_outgroup` treats the last entry as the outgroup and polarizes
/// substitutions against it.
pub fn write_variations<W: Write>(
    block: &FasBlock,
    has_outgroup: bool,
    writer: &mut W,
) -> anyhow::Result<()> {
    if block.entries.is_empty() {
        return Ok(());
    }

    let first = &block.entries[0];
    let trange = first.range();
    let t_ints_seq = crate::libs::alignment::seq_intspan(first.seq());

    let seqs: Vec<&[u8]> = block.entries.iter().map(|e| e.seq()).collect();
    let seq_count = seqs.len();
    if has_outgroup && seq_count < 2 {
        anyhow::bail!(
            "outgroup mode requires at least 2 sequences per block, got {}",
            seq_count
        );
    }

    let subs = if has_outgroup {
        let mut unpolarized = crate::libs::alignment::get_subs(&seqs[..(seq_count - 1)])?;
        crate::libs::alignment::polarize_subs(&mut unpolarized, seqs[seq_count - 1])?;
        unpolarized
    } else {
        crate::libs::alignment::get_subs(&seqs)?
    };

    for s in subs {
        let chr = trange.chr();
        let chr_pos = crate::libs::alignment::align_to_chr(
            &t_ints_seq,
            s.pos,
            trange.start,
            trange.strand(),
        )?;
        let var_rg = format!("{}:{}", chr, chr_pos);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            trange, chr, chr_pos, var_rg, s
        )?;
    }
    Ok(())
}

/// Write VCF rows for a single FasBlock.
///
/// `block_idx` is used only for error messages.
pub fn write_vcf_block<W: Write>(
    block: &FasBlock,
    block_idx: usize,
    writer: &mut W,
) -> anyhow::Result<()> {
    if block.entries.is_empty() {
        return Ok(());
    }

    let seqs: Vec<&[u8]> = block.entries.iter().map(|e| e.seq()).collect();
    let target_entry = &block.entries[0];
    let trange = target_entry.range();
    let t_ints_seq = crate::libs::alignment::seq_intspan(target_entry.seq());

    let subs = crate::libs::alignment::get_subs(&seqs)?;

    for s in subs {
        let chr = trange.chr();
        let chr_pos =
            crate::libs::alignment::align_to_chr(&t_ints_seq, s.pos, trange.start, trange.strand())
                .with_context(|| format!("align_to_chr at pos {} in block {}", s.pos, block_idx))?;

        let pos_idx = usize::try_from(s.pos).map_err(|_| {
            anyhow::anyhow!("invalid substitution pos {} in block {}", s.pos, block_idx)
        })?;
        let pos_idx = pos_idx.checked_sub(1).ok_or_else(|| {
            anyhow::anyhow!("invalid substitution pos {} in block {}", s.pos, block_idx)
        })?;
        if pos_idx >= seqs[0].len() {
            anyhow::bail!(
                "substitution pos {} out of range (seq len {}) in block {}",
                s.pos,
                seqs[0].len(),
                block_idx
            );
        }
        let ref_base = char::from(seqs[0][pos_idx]).to_ascii_uppercase();
        let alt_bases = crate::libs::alignment::vcf_alt_bases(&s);
        let sample_bases: Vec<u8> = seqs
            .iter()
            .map(|seq| seq.get(pos_idx).copied().unwrap_or(b'-'))
            .collect();

        crate::libs::fmt::vcf::write_snp_row(
            writer,
            chr,
            chr_pos,
            ref_base,
            &alt_bases,
            &sample_bases,
        )?;
    }
    Ok(())
}

/// Concatenate accumulated sequences and write them in FASTA or relaxed PHYLIP format.
pub fn write_concat_output<W: Write>(
    writer: &mut W,
    needed: &[String],
    seq_of: &std::collections::BTreeMap<String, String>,
    is_phylip: bool,
) -> anyhow::Result<()> {
    if needed.is_empty() {
        anyhow::bail!("no species specified for concat output");
    }
    if is_phylip {
        let length = seq_of.get(&needed[0]).map(|s| s.len()).unwrap_or(0);
        if length == 0 {
            anyhow::bail!(
                "PHYLIP output requires non-empty sequences, but all sequences are empty (check --required list and input blocks)"
            );
        }
        for name in needed {
            let v = seq_of
                .get(name)
                .ok_or_else(|| anyhow::anyhow!("name not found in concat results: {}", name))?;
            if v.len() != length {
                anyhow::bail!(
                    "PHYLIP requires equal-length sequences, but {} has length {} (expected {})",
                    name,
                    v.len(),
                    length
                );
            }
        }
        writeln!(writer, "{} {}", needed.len(), length)?;
        for name in needed {
            let v = seq_of
                .get(name)
                .ok_or_else(|| anyhow::anyhow!("name not found in concat results: {}", name))?;
            writeln!(writer, "{} {}", name, v)?;
        }
    } else {
        for name in needed {
            let v = seq_of
                .get(name)
                .ok_or_else(|| anyhow::anyhow!("name not found in concat results: {}", name))?;
            writeln!(writer, ">{}\n{}", name, v)?;
        }
    }
    Ok(())
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

/// Returns the file key used to group a FasBlock when splitting.
///
/// When `is_chr` is true, the key is `{name}.{chr}` using the first entry's
/// species name and chromosome. Otherwise it is the first entry's full range
/// string. Returns `None` for an empty block.
pub fn split_block_key(block: &FasBlock, is_chr: bool) -> Option<String> {
    let first = block.entries.first()?;
    let first_name = block.names.first()?;
    let key = if is_chr {
        format!("{}.{}", first_name, first.range().chr())
    } else {
        first.range().to_string()
    };
    Some(key)
}

/// Format one FasBlock for the `split` command.
///
/// Each entry is written as `>{header}\n{seq}\n`. When `is_simple` is true,
/// the header is reduced to the species name; otherwise the full range is
/// used. The returned string does not include the trailing blank line.
pub fn format_split_block(block: &FasBlock, is_simple: bool) -> anyhow::Result<String> {
    use std::fmt::Write;
    let mut out = String::new();
    for (idx, entry) in block.entries.iter().enumerate() {
        let header_owned = if is_simple {
            block.names.get(idx).cloned().unwrap_or_default()
        } else {
            entry.range().to_string()
        };
        let seq = std::str::from_utf8(entry.seq())?;
        writeln!(out, ">{}\n{}", header_owned, seq)?;
    }
    Ok(out)
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

    #[test]
    fn parse_fas_block_handles_whitespace_separator() {
        let input = ">S288c.I(+):13267-13287
TCGTCAGTTGGTTGACCATTA
   \n>S288c.I(+):185273-185334
GCATATAATATGAACCAATATCTA\n";
        let mut reader = BufReader::new(input.as_bytes());

        let block = crate::libs::fmt::fas::next_fas_block(&mut reader).unwrap();
        assert_eq!(block.entries.len(), 1, "first block should have one entry");

        let block = crate::libs::fmt::fas::next_fas_block(&mut reader).unwrap();
        assert_eq!(block.entries.len(), 1, "second block should have one entry");
    }
}
