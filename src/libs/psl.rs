use std::fmt;
use std::io;

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
            Self::reverse_range(&mut qs, &mut qe, psl.q_size);
        }

        let mut ts = psl.t_start;
        let mut te = psl.t_end;
        let t_strand_rev = if strand.len() >= 2 {
            strand.chars().nth(1) == Some('-')
        } else {
            false
        };

        if t_strand_rev {
            Self::reverse_range(&mut ts, &mut te, psl.t_size);
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

    fn reverse_range(start: &mut i32, end: &mut i32, size: u32) {
        let s = *start;
        let e = *end;
        *start = size as i32 - e;
        *end = size as i32 - s;
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
            self.q_starts[i] = self.q_size - (self.q_starts[i] + self.block_sizes[i]);
            self.t_starts[i] = self.t_size - (self.t_starts[i] + self.block_sizes[i]);
        }
    }

    pub fn is_protein(&self) -> bool {
        if self.block_count == 0 {
            return false;
        }
        let last = (self.block_count as usize) - 1;
        let t_strand = self.strand.chars().nth(1).unwrap_or('+');

        let t_end = self.t_end as u32;
        let t_start = self.t_start as u32;
        let t_size = self.t_size;
        let t_start_last = self.t_starts[last];
        let block_size_last = self.block_sizes[last];

        if t_strand == '+' {
            t_end == t_start_last + 3 * block_size_last
        } else if t_strand == '-' {
            t_start == t_size - (t_start_last + 3 * block_size_last)
        } else {
            false
        }
    }

    /// Reverse-complement a PSL alignment. This makes the target strand explicit.
    pub fn rc(&mut self) {
        let is_prot = self.is_protein();
        let mult = if is_prot { 3 } else { 1 };

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
            self.t_starts[i] = t_size - (self.t_starts[i] + mult * self.block_sizes[i]);
            self.q_starts[i] = q_size - (self.q_starts[i] + self.block_sizes[i]);
        }

        self.t_starts.reverse();
        self.q_starts.reverse();
        self.block_sizes.reverse();
    }

    pub fn score(&self) -> i32 {
        let is_prot = self.is_protein();
        let size_mul = if is_prot { 3 } else { 1 };
        (size_mul * (self.match_count + (self.rep_match >> 1))) as i32
            - (size_mul * self.mismatch_count) as i32
            - self.q_num_insert as i32
            - self.t_num_insert as i32
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
}

impl fmt::Display for Psl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut buf = Vec::new();
        self.write_to(&mut buf).map_err(|_| fmt::Error)?;
        let s = String::from_utf8_lossy(&buf);
        write!(f, "{}", s.trim_end())
    }
}

#[cfg(test)]
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
