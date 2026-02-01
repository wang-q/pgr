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
}
