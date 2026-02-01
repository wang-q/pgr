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

    pub fn write_to<W: io::Write>(&self, w: &mut W) -> io::Result<()> {
        write!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            self.match_count, self.mismatch_count, self.rep_match, self.n_count,
            self.q_num_insert, self.q_base_insert, self.t_num_insert, self.t_base_insert,
            self.strand, self.q_name, self.q_size, self.q_start, self.q_end,
            self.t_name, self.t_size, self.t_start, self.t_end, self.block_count
        )?;
        
        for s in &self.block_sizes { write!(w, "{},", s)?; }
        write!(w, "\t")?;
        for s in &self.q_starts { write!(w, "{},", s)?; }
        write!(w, "\t")?;
        for s in &self.t_starts { write!(w, "{},", s)?; }
        
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
