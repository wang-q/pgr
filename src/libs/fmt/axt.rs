use super::maf::{MafAli, MafComp, MafWriter};
use super::psl::Psl;
use std::collections::BTreeMap;
use std::io::{BufRead, Write};

#[derive(Debug, Clone, Default)]
pub struct Axt {
    pub id: u64, // The first number in the header line
    pub t_name: String,
    pub t_start: usize, // 0-based
    pub t_end: usize,   // 0-based, half-open
    pub t_strand: char, // AXT target strand is always '+'
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
            size: self.t_end.checked_sub(self.t_start).ok_or_else(|| {
                anyhow::anyhow!("t_end {} < t_start {}", self.t_end, self.t_start)
            })?,
            strand: self.t_strand,
            src_size: t_size,
            text: self.t_sym.clone(),
        };
        let q_comp = MafComp {
            src: format!("{}{}", q_prefix, self.q_name),
            start: self.q_start,
            size: self.q_end.checked_sub(self.q_start).ok_or_else(|| {
                anyhow::anyhow!("q_end {} < q_start {}", self.q_end, self.q_start)
            })?,
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
        let writer: Box<dyn std::io::Write> = Box::new(crate::libs::io::writer(output)?);
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
                    let w: Box<dyn std::io::Write> = Box::new(crate::libs::io::writer(path_str)?);
                    let mut w = MafWriter::new(w);
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
                    // Format: id tName tStart tEnd qName qStart qEnd qStrand score?
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
                        Ok(v) if v > 0 => v - 1, // 1-based to 0-based
                        Ok(_) => {
                            return Some(Err(anyhow::anyhow!(
                                "Invalid tStart: {} (AXT coordinates are 1-based)",
                                parts[2]
                            )))
                        }
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
                        Ok(v) if v > 0 => v - 1, // 1-based to 0-based
                        Ok(_) => {
                            return Some(Err(anyhow::anyhow!(
                                "Invalid qStart: {} (AXT coordinates are 1-based)",
                                parts[5]
                            )))
                        }
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
                        match parts[8].parse::<i32>() {
                            Ok(v) => Some(v),
                            Err(_) => {
                                return Some(Err(anyhow::anyhow!(
                                    "Invalid AXT score: {}",
                                    parts[8]
                                )))
                            }
                        }
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

    if let Some(score) = axt.score {
        writeln!(
            writer,
            "{} {} {} {} {} {} {} {} {}",
            axt.id, axt.t_name, t_start, t_end, axt.q_name, q_start, q_end, axt.q_strand, score
        )?;
    } else {
        writeln!(
            writer,
            "{} {} {} {} {} {} {} {}",
            axt.id, axt.t_name, t_start, t_end, axt.q_name, q_start, q_end, axt.q_strand
        )?;
    }

    writeln!(writer, "{}", axt.t_sym)?;
    writeln!(writer, "{}", axt.q_sym)?;
    writeln!(writer)?; // Blank line

    Ok(())
}

/// Sort key for `sort_axts`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxtSortBy {
    Target,
    Query,
    Score,
}

/// Sort axts in place by the given key. If `renumber`, reassign ids starting from 0.
pub fn sort_axts(axts: &mut [Axt], by: AxtSortBy, renumber: bool) {
    match by {
        AxtSortBy::Score => {
            // Sort by score descending (higher is better).
            axts.sort_by_key(|b| std::cmp::Reverse(b.score.unwrap_or(0)));
        }
        AxtSortBy::Query => {
            axts.sort_by(|a, b| a.q_name.cmp(&b.q_name).then(a.q_start.cmp(&b.q_start)));
        }
        AxtSortBy::Target => {
            axts.sort_by(|a, b| a.t_name.cmp(&b.t_name).then(a.t_start.cmp(&b.t_start)));
        }
    }

    if renumber {
        for (i, axt) in axts.iter_mut().enumerate() {
            axt.id = i as u64;
        }
    }
}

/// Convert AXT query coordinates (0-based) to forward-strand 1-based coordinates.
pub fn axt_query_to_forward_coords(
    q_start: usize,
    q_end: usize,
    q_strand: char,
    q_len: i32,
) -> anyhow::Result<(i32, i32)> {
    if q_strand == '-' {
        let q_s_1 = (q_start + 1) as i32;
        let q_e_1 = q_end as i32;
        if q_e_1 > q_len || q_s_1 > q_len {
            anyhow::bail!("AXT query coordinate exceeds q_len {}", q_len);
        }
        Ok((q_len - q_e_1 + 1, q_len - q_s_1 + 1))
    } else {
        Ok(((q_start + 1) as i32, q_end as i32))
    }
}

/// Convert an AXT stream to PSL, looking up sequence sizes from the given maps.
///
/// For each AXT record, query coordinates on the `-` strand are reversed to
/// forward-strand coordinates (per `axtToPsl.c` convention) before invoking
/// `Psl::from_align`. Records with invalid coordinates are skipped with a
/// `log::warn!`.
pub fn axt_to_psl<R: std::io::Read, W: Write>(
    reader: R,
    writer: &mut W,
    t_sizes: &BTreeMap<String, usize>,
    q_sizes: &BTreeMap<String, usize>,
) -> anyhow::Result<()> {
    let reader = AxtReader::new(reader);

    for result in reader {
        let axt = result?;

        let q_size = *q_sizes
            .get(&axt.q_name)
            .ok_or_else(|| anyhow::anyhow!("Query size not found for {}", axt.q_name))?;
        let t_size = *t_sizes
            .get(&axt.t_name)
            .ok_or_else(|| anyhow::anyhow!("Target size not found for {}", axt.t_name))?;

        // libs/axt.rs returns 0-based half-open coordinates
        let mut q_start = i32::try_from(axt.q_start)
            .map_err(|_| anyhow::anyhow!("q_start {} exceeds i32 range", axt.q_start))?;
        let mut q_end = i32::try_from(axt.q_end)
            .map_err(|_| anyhow::anyhow!("q_end {} exceeds i32 range", axt.q_end))?;
        let q_size_i32 = i32::try_from(q_size)
            .map_err(|_| anyhow::anyhow!("q_size {} exceeds i32 range", q_size))?;
        let t_start_i32 = i32::try_from(axt.t_start)
            .map_err(|_| anyhow::anyhow!("t_start {} exceeds i32 range", axt.t_start))?;
        let t_end_i32 = i32::try_from(axt.t_end)
            .map_err(|_| anyhow::anyhow!("t_end {} exceeds i32 range", axt.t_end))?;
        let q_size_u32 = u32::try_from(q_size)
            .map_err(|_| anyhow::anyhow!("q_size {} exceeds u32 range", q_size))?;
        let t_size_u32 = u32::try_from(t_size)
            .map_err(|_| anyhow::anyhow!("t_size {} exceeds u32 range", t_size))?;

        // axtToPsl.c logic: "if (axt->qStrand == '-') reverseIntRange(&qStart, &qEnd, qSize);"
        // This converts strand-relative coordinates (as in AXT) to positive strand coordinates
        // which pslFromAlign expects (so it can reverse them back internally).
        if axt.q_strand == '-' {
            crate::reverse_range(&mut q_start, &mut q_end, q_size_i32);
        }

        // Construct strand string for PSL (e.g. "-")
        // Note: PSL usually tracks target strand implicitly as +, so strand field is just query strand?
        // axtToPsl.c: strand[0] = axt->qStrand; strand[1] = '\0';
        // So it's just "+" or "-"
        let strand = axt.q_strand.to_string();

        if let Some(psl) = Psl::from_align(
            &axt.q_name,
            q_size_u32,
            q_start,
            q_end,
            &axt.q_sym,
            &axt.t_name,
            t_size_u32,
            t_start_i32,
            t_end_i32,
            &axt.t_sym,
            &strand,
        ) {
            psl.write_to(writer)?;
        } else {
            log::warn!(
                "skipping alignment (invalid coordinates): {} vs {}",
                axt.q_name,
                axt.t_name
            );
        }
    }

    Ok(())
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
