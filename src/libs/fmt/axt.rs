use std::collections::BTreeMap;
use std::io::BufRead;

use super::maf::{MafAli, MafComp, MafWriter};

#[derive(Debug, Clone, Default)]
pub struct Axt {
    pub id: u64, // The first number in the header line
    pub t_name: String,
    pub t_start: usize, // 0-based
    pub t_end: usize,   // 0-based, half-open
    pub t_strand: char,
    pub q_name: String,
    pub q_start: usize, // 0-based
    pub q_end: usize,   // 0-based, half-open
    pub q_strand: char,
    pub score: Option<i32>,
    pub t_sym: String,
    pub q_sym: String,
}

impl Axt {
    pub fn new() -> Self {
        Self::default()
    }

    /// Convert this AXT record into a MAF alignment block (`MafAli`).
    ///
    /// `t_sizes` / `q_sizes` provide the total sequence lengths (src_size in
    /// MAF). `t_prefix` / `q_prefix` are prepended to the source names.
    pub fn to_maf_ali(
        &self,
        t_sizes: &BTreeMap<String, usize>,
        q_sizes: &BTreeMap<String, usize>,
        t_prefix: &str,
        q_prefix: &str,
    ) -> anyhow::Result<MafAli> {
        let t_size = *t_sizes
            .get(&self.t_name)
            .ok_or_else(|| anyhow::anyhow!("Target size not found for {}", self.t_name))?;
        let q_size = *q_sizes
            .get(&self.q_name)
            .ok_or_else(|| anyhow::anyhow!("Query size not found for {}", self.q_name))?;

        let t_comp = MafComp {
            src: format!("{}{}", t_prefix, self.t_name),
            start: self.t_start,
            size: self.t_end - self.t_start,
            strand: self.t_strand,
            src_size: t_size,
            text: self.t_sym.clone(),
        };
        let q_comp = MafComp {
            src: format!("{}{}", q_prefix, self.q_name),
            start: self.q_start,
            size: self.q_end - self.q_start,
            strand: self.q_strand,
            src_size: q_size,
            text: self.q_sym.clone(),
        };

        Ok(MafAli {
            score: self.score.map(|s| s as f64),
            components: vec![t_comp, q_comp],
        })
    }
}

/// Convert an AXT file to MAF format.
///
/// * `input` / `output` — input AXT path and output MAF path (or `stdout`).
/// * `t_sizes` / `q_sizes` — sequence sizes for target/query.
/// * `t_prefix` / `q_prefix` — optional name prefixes.
/// * `t_split` — when true, write one MAF file per target sequence into the
///   `output` directory (expects input grouped by target name).
pub fn axt_to_maf(
    input: &str,
    output: &str,
    t_sizes: &BTreeMap<String, usize>,
    q_sizes: &BTreeMap<String, usize>,
    t_prefix: &str,
    q_prefix: &str,
    t_split: bool,
) -> anyhow::Result<()> {
    use std::collections::HashMap;
    use std::path::Path;

    let reader = crate::libs::io::reader(input)?;
    let axt_reader = AxtReader::new(reader);

    let mut current_t_name = String::new();
    let mut single_writer: Option<MafWriter<Box<dyn std::io::Write>>> = None;

    if t_split {
        if !Path::new(output).exists() {
            std::fs::create_dir_all(output)?;
        }
    } else {
        let writer = crate::libs::io::writer(output)?;
        let mut writer = MafWriter::new(writer);
        writer.write_header("blastz")?;
        single_writer = Some(writer);
    }

    let mut split_writers: HashMap<String, MafWriter<Box<dyn std::io::Write>>> = HashMap::new();

    for result in axt_reader {
        let axt = result?;

        let writer = if t_split {
            if axt.t_name != current_t_name {
                // C axtToMaf keeps only one file open and overwrites on tName change;
                // input is assumed grouped (sorted) by target name.
                if !split_writers.contains_key(&axt.t_name) {
                    split_writers.clear();
                    let path = Path::new(output).join(format!("{}.maf", axt.t_name));
                    let path_str = path.to_str().ok_or_else(|| {
                        anyhow::anyhow!("path is not valid UTF-8: {}", path.display())
                    })?;
                    let mut w = MafWriter::new(crate::libs::io::writer(path_str)?);
                    w.write_header("blastz")?;
                    split_writers.insert(axt.t_name.clone(), w);
                }
                current_t_name = axt.t_name.clone();
            }
            split_writers.get_mut(&axt.t_name).unwrap()
        } else {
            single_writer.as_mut().unwrap()
        };

        let ali = axt.to_maf_ali(t_sizes, q_sizes, t_prefix, q_prefix)?;
        writer.write_ali(&ali)?;
    }

    Ok(())
}

pub struct AxtReader<R> {
    reader: std::io::BufReader<R>,
    line_buf: String,
    pub headers: Vec<String>,
}

impl<R: std::io::Read> AxtReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            reader: std::io::BufReader::new(inner),
            line_buf: String::new(),
            headers: Vec::new(),
        }
    }

    fn read_line(&mut self) -> std::io::Result<usize> {
        self.line_buf.clear();
        self.reader.read_line(&mut self.line_buf)
    }
}

impl<R: std::io::Read> Iterator for AxtReader<R> {
    type Item = anyhow::Result<Axt>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Read first line (header)
            match self.read_line() {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    let line = self.line_buf.trim();
                    if line.is_empty() {
                        continue;
                    }
                    if line.starts_with('#') {
                        self.headers.push(line.to_string());
                        continue;
                    }

                    // Parse header
                    // Format: id tName tStart tEnd tStrand? qName qStart qEnd qStrand score?
                    // Example: 0 chr19 3001012 3001075 chr11 70568380 70568443 - 3500
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() < 8 {
                        return Some(Err(anyhow::anyhow!("Invalid AXT header: {}", line)));
                    }

                    let id: u64 = match parts[0].parse() {
                        Ok(v) => v,
                        Err(_) => {
                            return Some(Err(anyhow::anyhow!("Invalid AXT id: {}", parts[0])))
                        }
                    };

                    let t_name = parts[1].to_string();
                    let t_start: usize = match parts[2].parse::<usize>() {
                        Ok(v) => {
                            if v > 0 {
                                v - 1
                            } else {
                                0
                            }
                        } // 1-based to 0-based
                        Err(_) => {
                            return Some(Err(anyhow::anyhow!("Invalid tStart: {}", parts[2])))
                        }
                    };
                    let t_end: usize = match parts[3].parse::<usize>() {
                        Ok(v) => v,
                        Err(_) => return Some(Err(anyhow::anyhow!("Invalid tEnd: {}", parts[3]))),
                    };

                    let q_name = parts[4].to_string();
                    let q_start: usize = match parts[5].parse::<usize>() {
                        Ok(v) => {
                            if v > 0 {
                                v - 1
                            } else {
                                0
                            }
                        } // 1-based to 0-based
                        Err(_) => {
                            return Some(Err(anyhow::anyhow!("Invalid qStart: {}", parts[5])))
                        }
                    };
                    let q_end: usize = match parts[6].parse::<usize>() {
                        Ok(v) => v,
                        Err(_) => return Some(Err(anyhow::anyhow!("Invalid qEnd: {}", parts[6]))),
                    };

                    let q_strand = parts[7].chars().next().unwrap_or('?');

                    let score = if parts.len() > 8 {
                        parts[8].parse().ok()
                    } else {
                        None
                    };

                    let mut axt = Axt {
                        id,
                        t_name,
                        t_start,
                        t_end,
                        t_strand: '+',
                        q_name,
                        q_start,
                        q_end,
                        q_strand,
                        score,
                        t_sym: String::new(),
                        q_sym: String::new(),
                    };

                    // Read tSym
                    match self.read_line() {
                        Ok(0) => return Some(Err(anyhow::anyhow!("Unexpected EOF after header"))),
                        Ok(_) => axt.t_sym = self.line_buf.trim().to_string(),
                        Err(e) => return Some(Err(anyhow::Error::new(e))),
                    }

                    // Read qSym
                    match self.read_line() {
                        Ok(0) => return Some(Err(anyhow::anyhow!("Unexpected EOF after tSym"))),
                        Ok(_) => axt.q_sym = self.line_buf.trim().to_string(),
                        Err(e) => return Some(Err(anyhow::Error::new(e))),
                    }

                    if axt.t_sym.len() != axt.q_sym.len() {
                        return Some(Err(anyhow::anyhow!(
                            "Alignment lengths differ: {} vs {}",
                            axt.t_sym.len(),
                            axt.q_sym.len()
                        )));
                    }

                    return Some(Ok(axt));
                }
                Err(e) => return Some(Err(anyhow::Error::new(e))),
            }
        }
    }
}

pub fn write_axt<W: std::io::Write>(writer: &mut W, axt: &Axt) -> std::io::Result<()> {
    // Write header
    // Convert 0-based to 1-based
    let t_start = axt.t_start + 1;
    let t_end = axt.t_end;
    let q_start = axt.q_start + 1;
    let q_end = axt.q_end;

    let score_str = match axt.score {
        Some(s) => format!("{}", s),
        None => String::new(),
    };

    writeln!(
        writer,
        "{} {} {} {} {} {} {} {} {}",
        axt.id,
        axt.t_name,
        t_start,
        t_end,
        axt.q_name,
        q_start,
        q_end,
        axt.q_strand,
        score_str // Placeholder
    )?;

    writeln!(writer, "{}", axt.t_sym)?;
    writeln!(writer, "{}", axt.q_sym)?;
    writeln!(writer)?; // Blank line

    Ok(())
}

/// Convert AXT query coordinates (0-based) to forward-strand 1-based coordinates.
pub fn axt_query_to_forward_coords(
    q_start: usize,
    q_end: usize,
    q_strand: char,
    q_len: i32,
) -> (i32, i32) {
    if q_strand == '-' {
        let q_s_1 = (q_start + 1) as i32;
        let q_e_1 = q_end as i32;
        (q_len - q_e_1 + 1, q_len - q_s_1 + 1)
    } else {
        ((q_start + 1) as i32, q_end as i32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_axt() {
        let input = "\
0 chr19 3001012 3001075 chr11 70568380 70568443 - 3500
TCAGCTCTGGACGCAATACTCGCTCCCAGTCCAGATTCCTTCCTGATACTCGTCATGTGAGGA
TCTGTTCGTTGCACAT---TCGCTCCCAGTCCAGATTCCTTCCTGATACTCGTCATGTGAGGA

";
        let reader = AxtReader::new(input.as_bytes());
        let axts: Vec<Axt> = reader.collect::<Result<Vec<_>, _>>().unwrap();

        assert_eq!(axts.len(), 1);
        let a = &axts[0];
        assert_eq!(a.id, 0);
        assert_eq!(a.t_name, "chr19");
        assert_eq!(a.t_start, 3001011); // 3001012 - 1
        assert_eq!(a.t_end, 3001075);
        assert_eq!(a.q_name, "chr11");
        assert_eq!(a.q_start, 70568379); // 70568380 - 1
        assert_eq!(a.q_end, 70568443);
        assert_eq!(a.q_strand, '-');
        assert_eq!(a.score, Some(3500));
        assert_eq!(a.t_sym.len(), 63);
        assert_eq!(a.q_sym.len(), 63);
    }
}
