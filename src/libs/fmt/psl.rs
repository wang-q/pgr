use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::io::{self, BufRead, Write};

#[derive(Debug, Clone, Default)]
pub struct Psl {
    pub match_count: u32,
    pub mismatch_count: u32,
    pub rep_match: u32,
    pub n_count: u32,
    pub q_num_insert: u32,
    pub q_base_insert: i32,
    pub t_num_insert: u32,
    pub t_base_insert: i32,
    pub strand: String, // "+", "-", "++", "+-"
    pub q_name: String,
    pub q_size: u32,
    pub q_start: i32,
    pub q_end: i32,
    pub t_name: String,
    pub t_size: u32,
    pub t_start: i32,
    pub t_end: i32,
    pub block_count: u32,
    pub block_sizes: Vec<u32>,
    pub q_starts: Vec<u32>,
    pub t_starts: Vec<u32>,
}

impl Psl {
    pub fn new() -> Self {
        Default::default()
    }

    #[allow(clippy::too_many_arguments)]
    pub fn from_align(
        q_name: &str,
        q_size: u32,
        mut q_start: i32,
        mut q_end: i32,
        q_string: &str,
        t_name: &str,
        t_size: u32,
        mut t_start: i32,
        mut t_end: i32,
        t_string: &str,
        strand: &str,
    ) -> Option<Self> {
        let q_bytes = q_string.as_bytes();
        let t_bytes = t_string.as_bytes();
        let ali_len = q_bytes.len();

        if t_bytes.len() != ali_len {
            return None;
        }

        // trimAlignment logic
        let mut start_idx = 0;
        let mut end_idx = ali_len;

        // skip leading indels
        while start_idx < end_idx {
            let q = q_bytes[start_idx];
            let t = t_bytes[start_idx];
            if Self::is_del_char(q) || Self::is_del_char(t) {
                if !Self::is_del_char(q) {
                    q_start += 1;
                } else if !Self::is_del_char(t) {
                    t_start += 1;
                }
                start_idx += 1;
            } else {
                break;
            }
        }

        // skip trailing indels
        while end_idx > start_idx {
            let q = q_bytes[end_idx - 1];
            let t = t_bytes[end_idx - 1];
            if Self::is_del_char(q) || Self::is_del_char(t) {
                if !Self::is_del_char(q) {
                    q_end -= 1;
                } else if !Self::is_del_char(t) {
                    t_end -= 1;
                }
                end_idx -= 1;
            } else {
                break;
            }
        }

        if q_start == q_end || t_start == t_end {
            return None;
        }

        let mut psl = Psl {
            q_name: q_name.to_string(),
            q_size,
            q_start,
            q_end,
            t_name: t_name.to_string(),
            t_size,
            t_start,
            t_end,
            strand: strand.to_string(),
            ..Default::default()
        };

        let mut qs = psl.q_start;
        let mut qe = psl.q_end;
        if strand.starts_with('-') {
            crate::libs::alignment::coords::reverse_range(&mut qs, &mut qe, psl.q_size as i32);
        }

        let mut ts = psl.t_start;
        let mut te = psl.t_end;
        let t_strand_rev = if strand.len() >= 2 {
            strand.chars().nth(1) == Some('-')
        } else {
            false
        };

        if t_strand_rev {
            crate::libs::alignment::coords::reverse_range(&mut ts, &mut te, psl.t_size as i32);
        }

        let mut either_insert = false;
        // qe/te track current block END. qs/ts track current block START.
        // In C: qe = qs; te = ts;
        qe = qs;
        te = ts;

        let mut prev_q = 0;
        let mut prev_t = 0;

        for i in start_idx..end_idx {
            let q = q_bytes[i];
            let t = t_bytes[i];

            if Self::is_del_char(q) && Self::is_del_char(t) {
                continue;
            } else if Self::is_del_char(q) || Self::is_del_char(t) {
                if !either_insert {
                    psl.add_block(qs, qe, ts, te);
                    either_insert = true;
                }
                if !Self::is_del_char(q) {
                    qe += 1;
                }
                if !Self::is_del_char(t) {
                    te += 1;
                }
            } else {
                if either_insert {
                    qs = qe;
                    ts = te;
                    either_insert = false;
                }
                qe += 1;
                te += 1;
            }
            psl.accum_counts(prev_q, prev_t, q, t);
            prev_q = q;
            prev_t = t;
        }
        psl.add_block(qs, qe, ts, te);

        Some(psl)
    }

    fn is_del_char(c: u8) -> bool {
        matches!(c, b'-' | b'.' | b'=' | b'_')
    }

    fn add_block(&mut self, qs: i32, qe: i32, ts: i32, _te: i32) {
        let size = qe - qs;
        if size > 0 {
            self.block_count += 1;
            self.block_sizes.push(size as u32);
            self.q_starts.push(qs as u32);
            self.t_starts.push(ts as u32);
        }
    }

    fn accum_counts(&mut self, prev_q: u8, prev_t: u8, q: u8, t: u8) {
        if !Self::is_del_char(q) && !Self::is_del_char(t) {
            let qu = q.to_ascii_uppercase();
            let tu = t.to_ascii_uppercase();
            if q == b'N' || t == b'N' {
                // Strict 'N' check as in C
                self.n_count += 1;
            } else if qu == tu {
                if qu != q || tu != t {
                    self.rep_match += 1;
                } else {
                    self.match_count += 1;
                }
            } else {
                self.mismatch_count += 1;
            }
        } else if Self::is_del_char(q) && !Self::is_del_char(t) {
            self.t_base_insert += 1;
            if !Self::is_del_char(prev_q) {
                self.t_num_insert += 1;
            }
        } else if Self::is_del_char(t) && !Self::is_del_char(q) {
            self.q_base_insert += 1;
            if !Self::is_del_char(prev_t) {
                self.q_num_insert += 1;
            }
        }
    }

    pub fn swap(&mut self, no_rc: bool) {
        // Swap simple fields
        std::mem::swap(&mut self.q_base_insert, &mut self.t_base_insert);
        std::mem::swap(&mut self.t_num_insert, &mut self.q_num_insert);
        std::mem::swap(&mut self.q_name, &mut self.t_name);
        std::mem::swap(&mut self.q_size, &mut self.t_size);
        std::mem::swap(&mut self.q_start, &mut self.t_start);
        std::mem::swap(&mut self.q_end, &mut self.t_end);

        // Handle strand and blocks
        let q_strand = self.strand.chars().nth(0).unwrap_or('+');
        let t_strand = self.strand.chars().nth(1);

        if let Some(ts) = t_strand {
            // Translated
            self.strand = format!("{}{}", ts, q_strand);
            self.swap_blocks();
        } else if no_rc {
            // Untranslated with no reverse complement
            // psl->strand[1] = psl->strand[0];
            // psl->strand[0] = '+';
            self.strand = format!("+{}", q_strand);
            self.swap_blocks();
        } else {
            // Untranslated
            if q_strand == '+' {
                self.swap_blocks();
            } else {
                self.swap_rc_blocks();
                self.strand = "-".to_string();
            }
        }
    }

    fn swap_blocks(&mut self) {
        for i in 0..self.block_count as usize {
            std::mem::swap(&mut self.q_starts[i], &mut self.t_starts[i]);
        }
    }

    fn swap_rc_blocks(&mut self) {
        // Reverse arrays
        self.t_starts.reverse();
        self.q_starts.reverse();
        self.block_sizes.reverse();

        // Swap starts
        std::mem::swap(&mut self.t_starts, &mut self.q_starts);

        // Recalculate coordinates
        // qSize and tSize have already been swapped
        for i in 0..self.block_count as usize {
            self.q_starts[i] = self
                .q_size
                .saturating_sub(self.q_starts[i].saturating_add(self.block_sizes[i]));
            self.t_starts[i] = self
                .t_size
                .saturating_sub(self.t_starts[i].saturating_add(self.block_sizes[i]));
        }
    }

    pub fn is_protein(&self) -> bool {
        if self.block_count == 0 {
            return false;
        }
        let last = (self.block_count as usize) - 1;
        let t_strand = self.strand.chars().nth(1).unwrap_or('+');

        let t_end = u32::try_from(self.t_end).unwrap_or(0);
        let t_start = u32::try_from(self.t_start).unwrap_or(0);
        let t_size = self.t_size;
        let t_start_last = self.t_starts[last];
        let block_size_last = self.block_sizes[last];

        if t_strand == '+' {
            t_end == t_start_last.saturating_add(3u32.saturating_mul(block_size_last))
        } else if t_strand == '-' {
            t_start
                == t_size.saturating_sub(
                    t_start_last.saturating_add(3u32.saturating_mul(block_size_last)),
                )
        } else {
            false
        }
    }

    /// Reverse-complement a PSL alignment. This makes the target strand explicit.
    pub fn rc(&mut self) {
        let is_prot = self.is_protein();
        let mult: u32 = if is_prot { 3 } else { 1 };

        // swap strand, forcing target to have an explict strand
        let q_s = self.strand.chars().nth(0).unwrap_or('+');
        let t_s = self.strand.chars().nth(1).unwrap_or('+');

        let flip = |c| if c == '-' { '+' } else { '-' };
        let new_q_s = flip(q_s);
        let new_t_s = flip(t_s);
        self.strand = format!("{}{}", new_q_s, new_t_s);

        let t_size = self.t_size;
        let q_size = self.q_size;

        for i in 0..self.block_count as usize {
            self.t_starts[i] = t_size.saturating_sub(
                self.t_starts[i].saturating_add(mult.saturating_mul(self.block_sizes[i])),
            );
            self.q_starts[i] =
                q_size.saturating_sub(self.q_starts[i].saturating_add(self.block_sizes[i]));
        }

        self.t_starts.reverse();
        self.q_starts.reverse();
        self.block_sizes.reverse();
    }

    pub fn score(&self) -> i32 {
        let is_prot = self.is_protein();
        let size_mul: u64 = if is_prot { 3 } else { 1 };
        let raw = (size_mul * (self.match_count as u64 + (self.rep_match as u64 >> 1))) as i64
            - (size_mul * self.mismatch_count as u64) as i64
            - self.q_num_insert as i64
            - self.t_num_insert as i64;
        raw.try_into().unwrap_or(i32::MAX)
    }

    pub fn calc_aligned(&self) -> u32 {
        self.match_count + self.mismatch_count + self.rep_match + self.n_count
    }

    pub fn calc_match(&self) -> u32 {
        self.match_count + self.rep_match
    }

    pub fn calc_ident(&self) -> f32 {
        let aligned = self.calc_aligned();
        if aligned == 0 {
            0.0
        } else {
            self.calc_match() as f32 / aligned as f32
        }
    }

    pub fn calc_q_cover(&self) -> f32 {
        if self.q_size == 0 {
            0.0
        } else {
            self.calc_aligned() as f32 / self.q_size as f32
        }
    }

    pub fn calc_t_cover(&self) -> f32 {
        if self.t_size == 0 {
            0.0
        } else {
            self.calc_aligned() as f32 / self.t_size as f32
        }
    }

    pub fn calc_rep_match(&self) -> f32 {
        let aligned = self.calc_aligned();
        if aligned == 0 {
            0.0
        } else {
            self.rep_match as f32 / aligned as f32
        }
    }

    /// Fraction of query covered by this alignment (aligned bases / q_size).
    pub fn cover(&self) -> f32 {
        let aligned = self.match_count + self.mismatch_count + self.rep_match;
        if aligned == 0 {
            0.0
        } else {
            aligned as f32 / self.q_size as f32
        }
    }

    /// Fraction identity (matches+rep_matches / aligned bases).
    pub fn ident(&self) -> f32 {
        let aligned = self.match_count + self.mismatch_count + self.rep_match;
        if aligned == 0 {
            0.0
        } else {
            (self.match_count + self.rep_match) as f32 / aligned as f32
        }
    }
}

impl std::str::FromStr for Psl {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let fields: Vec<&str> = s.split('\t').collect();
        if fields.len() < 21 {
            return Err(anyhow::anyhow!("Invalid PSL line: fewer than 21 columns"));
        }

        let parse_u32 = |s: &str| {
            s.parse::<u32>()
                .map_err(|_| anyhow::anyhow!("Invalid u32: {}", s))
        };
        let parse_i32 = |s: &str| {
            s.parse::<i32>()
                .map_err(|_| anyhow::anyhow!("Invalid i32: {}", s))
        };
        let parse_vec = |s: &str| -> Result<Vec<u32>, anyhow::Error> {
            s.split(',')
                .filter(|v| !v.is_empty())
                .map(|v| {
                    v.parse::<u32>()
                        .map_err(|_| anyhow::anyhow!("Invalid array val: {}", v))
                })
                .collect()
        };

        let mut psl = Psl {
            match_count: parse_u32(fields[0])?,
            mismatch_count: parse_u32(fields[1])?,
            rep_match: parse_u32(fields[2])?,
            n_count: parse_u32(fields[3])?,
            q_num_insert: parse_u32(fields[4])?,
            q_base_insert: parse_i32(fields[5])?,
            t_num_insert: parse_u32(fields[6])?,
            t_base_insert: parse_i32(fields[7])?,
            strand: fields[8].to_string(),
            q_name: fields[9].to_string(),
            q_size: parse_u32(fields[10])?,
            q_start: parse_i32(fields[11])?,
            q_end: parse_i32(fields[12])?,
            t_name: fields[13].to_string(),
            t_size: parse_u32(fields[14])?,
            t_start: parse_i32(fields[15])?,
            t_end: parse_i32(fields[16])?,
            block_count: parse_u32(fields[17])?,
            block_sizes: parse_vec(fields[18])?,
            q_starts: parse_vec(fields[19])?,
            t_starts: parse_vec(fields[20])?,
        };

        // Ensure consistency between block_count and vector lengths
        let min_len = psl
            .block_sizes
            .len()
            .min(psl.q_starts.len())
            .min(psl.t_starts.len());

        if (psl.block_count as usize) != min_len {
            psl.block_count = min_len as u32;
            psl.block_sizes.truncate(min_len);
            psl.q_starts.truncate(min_len);
            psl.t_starts.truncate(min_len);
        }

        Ok(psl)
    }
}

impl Psl {
    pub fn write_to<W: io::Write>(&self, w: &mut W) -> io::Result<()> {
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            self.match_count,
            self.mismatch_count,
            self.rep_match,
            self.n_count,
            self.q_num_insert,
            self.q_base_insert,
            self.t_num_insert,
            self.t_base_insert,
            self.strand,
            self.q_name,
            self.q_size,
            self.q_start,
            self.q_end,
            self.t_name,
            self.t_size,
            self.t_start,
            self.t_end,
            self.block_count
        )?;

        for s in &self.block_sizes {
            write!(w, "{},", s)?;
        }
        write!(w, "\t")?;
        for s in &self.q_starts {
            write!(w, "{},", s)?;
        }
        write!(w, "\t")?;
        for s in &self.t_starts {
            write!(w, "{},", s)?;
        }

        writeln!(w)?;
        Ok(())
    }

    /// Write this PSL record to a writer in UCSC Chain format.
    pub fn write_chain<W: io::Write>(&self, writer: &mut W, chain_id: u64) -> io::Result<()> {
        let score = self.score();
        let q_strand_char = self.strand.chars().next().unwrap_or('+');

        // Chain format: tStrand is always +, qStrand can be + or -.
        // If qStrand is -, qStart/qEnd are relative to reverse end.
        let (q_start, q_end) = if q_strand_char == '-' {
            crate::libs::alignment::reverse_range_pair(self.q_start, self.q_end, self.q_size as i32)
        } else {
            (self.q_start, self.q_end)
        };

        writeln!(
            writer,
            "chain {} {} {} + {} {} {} {} {} {} {} {}",
            score,
            self.t_name,
            self.t_size,
            self.t_start,
            self.t_end,
            self.q_name,
            self.q_size,
            q_strand_char,
            q_start,
            q_end,
            chain_id
        )?;

        // Write blocks
        for i in 0..self.block_count as usize {
            let size = self.block_sizes[i];
            write!(writer, "{}", size)?;

            if i < (self.block_count as usize) - 1 {
                let dt = self.t_starts[i + 1].saturating_sub(self.t_starts[i].saturating_add(size));
                let dq = self.q_starts[i + 1].saturating_sub(self.q_starts[i].saturating_add(size));
                write!(writer, "\t{}\t{}", dt, dq)?;
            }
            writeln!(writer)?;
        }

        Ok(())
    }
}

impl fmt::Display for Psl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut buf = Vec::new();
        self.write_to(&mut buf).map_err(|_| fmt::Error)?;
        let s = String::from_utf8_lossy(&buf);
        write!(f, "{}", s.trim_end())
    }
}

/// Accumulated statistics for PSL alignments.
#[derive(Debug, Clone, Default)]
pub struct SumStats {
    pub q_name: String,
    pub query_cnt: u32,
    pub min_q_size: u32,
    pub max_q_size: u32,
    pub total_q_size: u64,
    pub aln_cnt: u32,
    pub total_align: u64,
    pub total_match: u64,
    pub total_rep_match: u64,
    pub min_ident: f32,
    pub max_ident: f32,
    pub min_q_cover: f32,
    pub max_q_cover: f32,
    pub min_t_cover: f32,
    pub max_t_cover: f32,
    pub min_rep_match: f32,
    pub max_rep_match: f32,
}

impl SumStats {
    pub fn new(q_name: &str, q_size: u32) -> Self {
        Self {
            q_name: q_name.to_string(),
            query_cnt: 1,
            min_q_size: q_size,
            max_q_size: q_size,
            total_q_size: 0,
            aln_cnt: 0,
            ..Default::default()
        }
    }

    pub fn accumulate(&mut self, psl: &Psl) {
        let ident = psl.calc_ident();
        let q_cover = psl.calc_q_cover();
        let t_cover = psl.calc_t_cover();
        let rep_match = psl.calc_rep_match();

        self.total_q_size += psl.q_size as u64;

        if self.aln_cnt == 0 {
            self.min_ident = ident;
            self.max_ident = ident;
            self.min_q_cover = q_cover;
            self.max_q_cover = q_cover;
            self.min_t_cover = t_cover;
            self.max_t_cover = t_cover;
            self.min_rep_match = rep_match;
            self.max_rep_match = rep_match;

            self.min_q_size = self.min_q_size.min(psl.q_size);
            self.max_q_size = self.max_q_size.max(psl.q_size);
        } else {
            self.min_q_size = self.min_q_size.min(psl.q_size);
            self.max_q_size = self.max_q_size.max(psl.q_size);

            self.min_ident = self.min_ident.min(ident);
            self.max_ident = self.max_ident.max(ident);

            self.min_q_cover = self.min_q_cover.min(q_cover);
            self.max_q_cover = self.max_q_cover.max(q_cover);

            self.min_t_cover = self.min_t_cover.min(t_cover);
            self.max_t_cover = self.max_t_cover.max(t_cover);

            self.min_rep_match = self.min_rep_match.min(rep_match);
            self.max_rep_match = self.max_rep_match.max(rep_match);
        }

        self.total_align += psl.calc_aligned() as u64;
        self.total_match += psl.calc_match() as u64;
        self.total_rep_match += psl.rep_match as u64;
        self.aln_cnt += 1;
    }

    pub fn merge(&mut self, other: &SumStats) {
        if self.aln_cnt == 0 {
            self.min_q_size = other.min_q_size;
            self.max_q_size = other.max_q_size;
            self.min_ident = other.min_ident;
            self.max_ident = other.max_ident;
            self.min_q_cover = other.min_q_cover;
            self.max_q_cover = other.max_q_cover;
            self.min_t_cover = other.min_t_cover;
            self.max_t_cover = other.max_t_cover;
            self.min_rep_match = other.min_rep_match;
            self.max_rep_match = other.max_rep_match;
        } else if other.aln_cnt > 0 {
            self.min_q_size = self.min_q_size.min(other.min_q_size);
            self.max_q_size = self.max_q_size.max(other.max_q_size);
            self.min_ident = self.min_ident.min(other.min_ident);
            self.max_ident = self.max_ident.max(other.max_ident);
            self.min_q_cover = self.min_q_cover.min(other.min_q_cover);
            self.max_q_cover = self.max_q_cover.max(other.max_q_cover);
            self.min_t_cover = self.min_t_cover.min(other.min_t_cover);
            self.max_t_cover = self.max_t_cover.max(other.max_t_cover);
            self.min_rep_match = self.min_rep_match.min(other.min_rep_match);
            self.max_rep_match = self.max_rep_match.max(other.max_rep_match);
        }

        self.query_cnt += other.query_cnt;
        self.total_q_size += other.total_q_size;
        self.total_align += other.total_align;
        self.total_match += other.total_match;
        self.total_rep_match += other.total_rep_match;
        self.aln_cnt += other.aln_cnt;
    }

    pub fn mean_ident(&self) -> f32 {
        if self.total_align == 0 {
            0.0
        } else {
            self.total_match as f32 / self.total_align as f32
        }
    }

    pub fn mean_q_size(&self) -> u32 {
        if self.query_cnt == 0 {
            0
        } else {
            (self.total_q_size / self.query_cnt as u64) as u32
        }
    }

    pub fn mean_q_cover(&self) -> f32 {
        if self.total_q_size == 0 {
            0.0
        } else {
            self.total_align as f32 / self.total_q_size as f32
        }
    }

    pub fn mean_rep_match(&self) -> f32 {
        if self.total_align == 0 {
            0.0
        } else {
            self.total_rep_match as f32 / self.total_align as f32
        }
    }
}

/// Statistics output mode for `run_stats`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PslStatsMode {
    PerAlignment,
    PerQuery,
    Overall,
}

/// Options for `run_stats`.
#[derive(Debug, Clone)]
pub struct PslStatsOptions {
    pub mode: PslStatsMode,
    pub tsv: bool,
}

impl Default for PslStatsOptions {
    fn default() -> Self {
        Self {
            mode: PslStatsMode::PerAlignment,
            tsv: false,
        }
    }
}

/// Read a queries TSV (q_name<TAB>q_size) into a map of pre-initialized SumStats.
pub fn read_queries<R: BufRead>(reader: R) -> anyhow::Result<HashMap<String, SumStats>> {
    let mut tbl: HashMap<String, SumStats> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            log::warn!("skipping malformed queries line: {}", line);
            continue;
        }
        let q_name = parts[0].to_string();
        let q_size: u32 = parts[1].parse()?;
        tbl.insert(q_name.clone(), SumStats::new(&q_name, q_size));
    }
    Ok(tbl)
}

/// Run PSL statistics collection and write formatted output.
///
/// If `queries` is provided, only queries present in the map are counted
/// (and queries with zero alignments are still emitted in per-query/overall
/// modes). If `queries` is `None`, all records in the reader are counted.
pub fn run_stats<R: BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    opts: &PslStatsOptions,
    queries: Option<HashMap<String, SumStats>>,
) -> anyhow::Result<()> {
    let mut query_stats_tbl: HashMap<String, SumStats> = queries.unwrap_or_default();

    match opts.mode {
        PslStatsMode::PerQuery | PslStatsMode::Overall => {
            let has_queries = !query_stats_tbl.is_empty();
            for psl in iter_psl(reader) {
                let psl = psl?;
                if has_queries {
                    if let Some(entry) = query_stats_tbl.get_mut(&psl.q_name) {
                        entry.accumulate(&psl);
                    }
                } else {
                    let entry = query_stats_tbl
                        .entry(psl.q_name.clone())
                        .or_insert_with(|| SumStats::new(&psl.q_name, psl.q_size));
                    entry.accumulate(&psl);
                }
            }

            if opts.mode == PslStatsMode::PerQuery {
                if !opts.tsv {
                    write!(writer, "#")?;
                }
                writeln!(writer, "qName\tqSize\talnCnt\tminIdent\tmaxIdent\tmeanIdent\tminQCover\tmaxQCover\tmeanQCover\tminRepMatch\tmaxRepMatch\tmeanRepMatch\tminTCover\tmaxTCover")?;

                let mut keys: Vec<String> = query_stats_tbl.keys().cloned().collect();
                keys.sort();

                for k in keys {
                    let s = &query_stats_tbl[&k];
                    writeln!(writer, "{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                        s.q_name, s.min_q_size, s.aln_cnt,
                        s.min_ident, s.max_ident, s.mean_ident(),
                        s.min_q_cover, s.max_q_cover, s.mean_q_cover(),
                        s.min_rep_match, s.max_rep_match, s.mean_rep_match(),
                        s.min_t_cover, s.max_t_cover
                    )?;
                }
            } else {
                // Overall mode
                let mut os = SumStats::default();
                let mut aligned1 = 0;
                let mut aligned_n = 0;

                for s in query_stats_tbl.values() {
                    os.merge(s);

                    if s.aln_cnt == 1 {
                        aligned1 += 1;
                    } else if s.aln_cnt > 1 {
                        aligned_n += 1;
                    }
                }

                if !opts.tsv {
                    write!(writer, "#")?;
                }
                writeln!(writer, "queryCnt\tminQSize\tmaxQSize\tmeanQSize\talnCnt\tminIdent\tmaxIdent\tmeanIdent\tminQCover\tmaxQCover\tmeanQCover\tminRepMatch\tmaxRepMatch\tmeanRepMatch\tminTCover\tmaxTCover\taligned\taligned1\talignedN\ttotalAlignedSize")?;

                writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}",
                    os.query_cnt, os.min_q_size, os.max_q_size, os.mean_q_size(),
                    os.aln_cnt,
                    os.min_ident, os.max_ident, os.mean_ident(),
                    os.min_q_cover, os.max_q_cover, os.mean_q_cover(),
                    os.min_rep_match, os.max_rep_match, os.mean_rep_match(),
                    os.min_t_cover, os.max_t_cover,
                    aligned1 + aligned_n, aligned1, aligned_n,
                    os.total_align
                )?;
            }
        }
        PslStatsMode::PerAlignment => {
            if !opts.tsv {
                write!(writer, "#")?;
            }
            writeln!(
                writer,
                "qName\tqSize\ttName\ttStart\ttEnd\tident\tqCover\trepMatch\ttCover"
            )?;

            let has_queries = !query_stats_tbl.is_empty();
            for psl in iter_psl(reader) {
                let psl = psl?;
                if has_queries {
                    if let Some(entry) = query_stats_tbl.get_mut(&psl.q_name) {
                        writeln!(
                            writer,
                            "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                            psl.q_name,
                            psl.q_size,
                            psl.t_name,
                            psl.t_start,
                            psl.t_end,
                            psl.calc_ident(),
                            psl.calc_q_cover(),
                            psl.calc_rep_match(),
                            psl.calc_t_cover()
                        )?;
                        entry.aln_cnt += 1;
                    }
                } else {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                        psl.q_name,
                        psl.q_size,
                        psl.t_name,
                        psl.t_start,
                        psl.t_end,
                        psl.calc_ident(),
                        psl.calc_q_cover(),
                        psl.calc_rep_match(),
                        psl.calc_t_cover()
                    )?;
                }
            }

            if has_queries {
                let mut keys: Vec<String> = query_stats_tbl.keys().cloned().collect();
                keys.sort();

                for k in keys {
                    let s = &query_stats_tbl[&k];
                    if s.aln_cnt == 0 {
                        writeln!(
                            writer,
                            "{}\t{}\t\t0\t0\t0.0000\t0.0000\t0.0000\t0.0000",
                            s.q_name, s.min_q_size
                        )?;
                    }
                }
            }
        }
    }

    Ok(())
}

/// Iterate PSL records from a buffered reader, skipping header lines and
/// comments. Recognized headers: lines beginning with `#`, `psLayout`,
/// `match`, or `------` (UCSC pslLayout convention); empty lines are also
/// skipped. Unparseable lines are reported as errors.
pub fn iter_psl<R: io::BufRead>(reader: R) -> impl Iterator<Item = anyhow::Result<Psl>> {
    use std::str::FromStr;
    reader.lines().filter_map(|line| match line {
        Ok(line) => {
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("psLayout")
                || line.starts_with("match")
                || line.starts_with("------")
            {
                None
            } else {
                match Psl::from_str(&line) {
                    Ok(psl) => Some(Ok(psl)),
                    Err(err) => Some(Err(anyhow::anyhow!("invalid PSL line: {}", err))),
                }
            }
        }
        Err(err) => Some(Err(anyhow::anyhow!("read error: {}", err))),
    })
}

/// Parse a `chr:start-end` subrange name (1-based, inclusive) into
/// `(chr, start, end)` with `start`/`end` as u32. Returns `None` if `name`
/// is not a valid range.
pub fn parse_subrange(name: &str) -> Option<(String, u32, u32)> {
    let rg = intspan::Range::from_str(name);
    if rg.is_valid() {
        return Some((rg.chr().to_string(), *rg.start() as u32, *rg.end() as u32));
    }
    None
}

/// Compute (min, max) of `func` over a slice of Psl records.
pub fn calc_spread<F>(psls: &[Psl], func: F) -> (f32, f32)
where
    F: Fn(&Psl) -> f32,
{
    let mut min_val = f32::MAX;
    let mut max_val = f32::MIN;

    for psl in psls {
        let val = func(psl);
        if val < min_val {
            min_val = val;
        }
        if val > max_val {
            max_val = val;
        }
    }

    // Handle case where psls is empty (shouldn't happen here)
    if min_val == f32::MAX {
        (0.0, 0.0)
    } else {
        (min_val, max_val)
    }
}

impl Psl {
    /// Lift query coordinates from a fragment subrange to genomic coordinates.
    ///
    /// `sizes` maps chromosome name → real sequence size. Returns `true` if
    /// the query was lifted, `false` if skipped (no subrange, missing size,
    /// or subrange exceeds real size).
    pub fn lift_query(&mut self, sizes: &BTreeMap<String, i32>) -> bool {
        let (name_part, start, end) = match parse_subrange(&self.q_name) {
            Some(v) => v,
            None => return false,
        };
        let start_0 = start.saturating_sub(1);
        let end_0 = end;

        let real_size_i32 = match sizes.get(&name_part).copied() {
            Some(v) => v,
            None => {
                log::warn!("No sizes provided for {name_part}. Skipping query lift.");
                return false;
            }
        };
        let real_size = real_size_i32 as u32;

        if end_0 > real_size {
            log::warn!(
                "Subrange end {end_0} > sequence size {real_size} for {}. Skipping query lift.",
                self.q_name
            );
            return false;
        }

        let is_neg = self.strand.starts_with('-');
        self.q_name = name_part;
        self.q_size = real_size;
        let offset = if is_neg { real_size - end_0 } else { start_0 };
        let offset_i32 = i32::try_from(offset).unwrap_or(i32::MAX);
        self.q_start = self.q_start.saturating_add(offset_i32);
        self.q_end = self.q_end.saturating_add(offset_i32);
        for q_start in &mut self.q_starts {
            *q_start += offset;
        }
        true
    }

    /// Lift target coordinates from a fragment subrange to genomic coordinates.
    ///
    /// `sizes` maps chromosome name → real sequence size. Returns `true` if
    /// the target was lifted, `false` if skipped.
    pub fn lift_target(&mut self, sizes: &BTreeMap<String, i32>) -> bool {
        let (name_part, start, end) = match parse_subrange(&self.t_name) {
            Some(v) => v,
            None => return false,
        };
        let start_0 = start.saturating_sub(1);
        let end_0 = end;

        let real_size_i32 = match sizes.get(&name_part).copied() {
            Some(v) => v,
            None => {
                log::warn!("No sizes provided for {name_part}. Skipping target lift.");
                return false;
            }
        };
        let real_size = real_size_i32 as u32;

        if end_0 > real_size {
            log::warn!(
                "Subrange end {end_0} > sequence size {real_size} for {}. Skipping target lift.",
                self.t_name
            );
            return false;
        }

        let is_neg = if self.strand.len() >= 2 {
            self.strand.chars().nth(1).unwrap() == '-'
        } else {
            false
        };
        self.t_name = name_part;
        self.t_size = real_size;
        let offset = if is_neg { real_size - end_0 } else { start_0 };
        let offset_i32 = i32::try_from(offset).unwrap_or(i32::MAX);
        self.t_start = self.t_start.saturating_add(offset_i32);
        self.t_end = self.t_end.saturating_add(offset_i32);
        for t_start in &mut self.t_starts {
            *t_start += offset;
        }
        true
    }
}

/// Extract block ranges (1-based inclusive) from a PSL record as "name:start-end" strings.
pub fn psl_block_ranges(psl: &Psl, target: bool) -> Vec<String> {
    let (name, size, starts, is_neg) = if target {
        let is_neg = if psl.strand.len() >= 2 {
            psl.strand.chars().nth(1) == Some('-')
        } else {
            false
        };
        (psl.t_name.as_str(), psl.t_size, &psl.t_starts, is_neg)
    } else {
        let is_neg = psl.strand.starts_with('-');
        (psl.q_name.as_str(), psl.q_size, &psl.q_starts, is_neg)
    };

    let mut ranges = Vec::new();
    for (i, &start) in starts.iter().enumerate() {
        let len = psl.block_sizes[i];
        let end = start + len;
        let (final_start, final_end) = if is_neg {
            crate::libs::alignment::reverse_range_1based_pair(
                (start + 1) as usize,
                end as usize,
                size as usize,
            )
        } else {
            ((start + 1) as usize, end as usize)
        };
        ranges.push(format!("{}:{}-{}", name, final_start, final_end));
    }
    ranges
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    #[test]
    fn test_psl_display() {
        let psl = Psl {
            match_count: 59,
            mismatch_count: 13,
            rep_match: 0,
            n_count: 0,
            q_num_insert: 2,
            q_base_insert: 3,
            t_num_insert: 1,
            t_base_insert: 1,
            strand: "+".to_string(),
            q_name: "query".to_string(),
            q_size: 100,
            q_start: 10,
            q_end: 90,
            t_name: "target".to_string(),
            t_size: 200,
            t_start: 50,
            t_end: 130,
            block_count: 2,
            block_sizes: vec![40, 40],
            q_starts: vec![10, 50],
            t_starts: vec![50, 90],
        };

        let output = format!("{}", psl);
        // Note: write_to adds a newline, but Display implementation trims it.
        // Arrays are comma-separated with a trailing comma.
        let expected = "59\t13\t0\t0\t2\t3\t1\t1\t+\tquery\t100\t10\t90\ttarget\t200\t50\t130\t2\t40,40,\t10,50,\t50,90,";
        assert_eq!(output, expected);
    }

    #[test]
    fn test_parse_valid() {
        let line = "59\t13\t0\t0\t2\t3\t1\t1\t+\tquery\t100\t10\t90\ttarget\t200\t50\t130\t2\t40,40,\t10,50,\t50,90,";
        let psl: Psl = line.parse().unwrap();
        assert_eq!(psl.match_count, 59);
        assert_eq!(psl.block_count, 2);
        assert_eq!(psl.block_sizes, vec![40, 40]);
        assert_eq!(psl.q_starts, vec![10, 50]);
        assert_eq!(psl.t_starts, vec![50, 90]);
    }

    #[test]
    fn test_parse_invalid() {
        let line = "59\t13"; // Too short
        let res: Result<Psl, _> = line.parse();
        assert!(res.is_err());

        let line = "invalid\t13\t0\t0\t2\t3\t1\t1\t+\tquery\t100\t10\t90\ttarget\t200\t50\t130\t2\t40,40,\t10,50,\t50,90,";
        let res: Result<Psl, _> = line.parse();
        assert!(res.is_err());
    }

    #[test]
    fn test_score_dna() {
        // match=10, mismatch=2, rep=0, ins=0 -> 10 - 2 = 8
        let mut psl = Psl::default();
        psl.match_count = 10;
        psl.mismatch_count = 2;
        psl.q_size = 100;
        psl.t_size = 100;
        // make sure it's not protein
        psl.block_count = 1;
        psl.block_sizes = vec![10];
        psl.t_starts = vec![0];
        psl.t_start = 0;
        psl.t_end = 10;
        psl.strand = "+".to_string();

        assert_eq!(psl.score(), 8);
    }

    #[test]
    fn test_calc_ident() {
        let mut psl = Psl::default();
        psl.match_count = 90;
        psl.mismatch_count = 10;
        // aligned = 100. ident = 90/100 = 0.9
        assert_eq!(psl.calc_ident(), 0.9);
        assert_eq!(psl.calc_q_cover(), 0.0); // q_size is 0

        psl.q_size = 100;
        assert_eq!(psl.calc_q_cover(), 1.0);
    }

    #[test]
    fn test_swap() {
        let mut psl = Psl::default();
        psl.q_name = "q".to_string();
        psl.t_name = "t".to_string();
        psl.q_size = 100;
        psl.t_size = 200;
        psl.strand = "+".to_string();
        psl.block_count = 1;
        psl.block_sizes = vec![10];
        psl.q_starts = vec![0];
        psl.t_starts = vec![0];
        psl.q_start = 0;
        psl.q_end = 10;
        psl.t_start = 0;
        psl.t_end = 10;

        psl.swap(false);
        assert_eq!(psl.q_name, "t");
        assert_eq!(psl.t_name, "q");
        assert_eq!(psl.q_size, 200);
        assert_eq!(psl.t_size, 100);
    }

    #[test]
    fn test_rc() {
        // Simple case: 1 block, length 10.
        // T: 0-10 (size 100) -> RC: 100-10 = 90, 100-0 = 100. New start 90.
        let mut psl = Psl::default();
        psl.t_size = 100;
        psl.q_size = 50;
        psl.block_count = 1;
        psl.block_sizes = vec![10];
        psl.t_starts = vec![0];
        psl.q_starts = vec![0];
        psl.strand = "++".to_string();

        psl.rc();

        // strand should flip chars: ++ -> --
        assert_eq!(psl.strand, "--");
        // t_start: size(100) - (0 + 10) = 90
        assert_eq!(psl.t_starts[0], 90);
        // q_start: size(50) - (0 + 10) = 40
        assert_eq!(psl.q_starts[0], 40);
    }

    #[test]
    fn test_from_align() {
        // q: AC-G
        // t: ACTG
        let q_seq = "AC-G";
        let t_seq = "ACTG";
        let psl = Psl::from_align("q", 3, 0, 3, q_seq, "t", 4, 0, 4, t_seq, "+").unwrap();

        assert_eq!(psl.block_count, 2);
        assert_eq!(psl.block_sizes, vec![2, 1]); // AC, G
        assert_eq!(psl.q_starts, vec![0, 2]);
        assert_eq!(psl.t_starts, vec![0, 3]); // T is at index 2 in target, G is at 3
        assert_eq!(psl.match_count, 3); // A, C, G
        assert_eq!(psl.t_num_insert, 1);
        assert_eq!(psl.t_base_insert, 1);
    }

    #[test]
    fn test_is_protein() {
        let mut psl = Psl::default();
        psl.block_count = 1;
        psl.block_sizes = vec![10];
        psl.t_starts = vec![0];
        psl.t_start = 0;
        psl.t_end = 30; // 3 * 10
        psl.strand = "+".to_string();
        psl.t_size = 100;

        assert!(psl.is_protein());

        psl.t_end = 10;
        assert!(!psl.is_protein());
    }
}
