use std::io::{self, Write};

use crate::libs::io::LinesRef;
use crate::libs::paf::cigar::{
    block_identity, cigar_from_alignment, cigar_stats, format_cigar, gap_compressed_identity,
};
use crate::libs::paf::record::PafRecord;

// MAF
// https://genome.ucsc.edu/FAQ/FAQformat.html#format5
// https://github.com/joelarmstrong/maf_stream/blob/master/multiple_alignment_format/src/parser.rs

/// One component (an `s` line) of a MAF alignment block.
#[derive(Debug, PartialEq, Eq, Clone, Default)]
pub struct MafComp {
    /// The sequence name.
    pub src: String,
    /// Start of the aligned region within this sequence (0-based).
    pub start: usize,
    /// Length of the aligned region (not including gaps).
    pub size: usize,
    /// Which strand the aligned sequence is on.
    pub strand: char,
    /// The total length of this sequence (including regions outside this alignment).
    pub src_size: usize,
    /// Actual sequence of bases/amino acids, including gaps.
    pub text: String,
}

impl MafComp {
    /// Build a range string `src(strand):start-end` (1-based inclusive).
    pub fn to_range(&self) -> String {
        // adjust coordinates to be one-based inclusive
        let mut start = self.start + 1;
        let mut end = start + self.size - 1;

        // If the strand field is "-" then this is the start relative to the reverse-complemented source sequence
        if self.strand == '-' {
            crate::libs::alignment::coords::reverse_range_1based(
                &mut start,
                &mut end,
                self.src_size,
            );
        }

        format!("{}({}):{}-{}", self.src, self.strand, start, end)
    }
}

/// A MAF alignment block (an `a` line plus its `s`/`i`/`e`/`q` lines).
#[derive(Debug, PartialEq, Default)]
pub struct MafAli {
    /// Score extracted from the `a` line (`score=` field).
    pub score: Option<f64>,
    pub components: Vec<MafComp>,
}

/// Parse a strand token (`+` or `-`) into a `char`.
fn parse_strand(strand: &str) -> Result<char, io::Error> {
    match strand {
        "+" => Ok('+'),
        "-" => Ok('-'),
        _ => Err(io::Error::other("Strand not valid")),
    }
}

/// Read the next `MafAli` block from `input`. Returns `Err` on EOF or malformed input.
pub fn next_maf_block<T: io::BufRead + ?Sized>(mut input: &mut T) -> Result<MafAli, io::Error> {
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
                // MAF comment
                continue;
            } else if line.starts_with('a') {
                // Start of a block
                header = Some(line);
                break;
            } else {
                // Shouldn't see this.
                return Err(io::Error::other("Unexpected line"));
            }
        }
    }
    let block = parse_maf_block(
        header.ok_or(io::Error::other("EOF"))?,
        LinesRef { buf: &mut input },
    )?;
    Ok(block)
}

fn parse_s_line(fields: &mut Vec<&str>, components: &mut Vec<MafComp>) -> Result<(), io::Error> {
    let text = fields
        .pop()
        .ok_or(io::Error::other("s line incomplete"))?
        .to_string();
    let src_size = fields
        .pop()
        .ok_or(io::Error::other("s line incomplete"))
        .and_then(|s| {
            s.parse::<usize>()
                .map_err(|_| io::Error::other("invalid sequence size"))
        })?;
    let strand = fields
        .pop()
        .ok_or(io::Error::other("s line incomplete"))
        .and_then(parse_strand)?;
    let size = fields
        .pop()
        .ok_or(io::Error::other("s line incomplete"))
        .and_then(|s| {
            s.parse::<usize>()
                .map_err(|_| io::Error::other("invalid aligned length"))
        })?;
    let start = fields
        .pop()
        .ok_or(io::Error::other("s line incomplete"))
        .and_then(|s| {
            s.parse::<usize>()
                .map_err(|_| io::Error::other("invalid start"))
        })?;
    let src = fields
        .pop()
        .ok_or(io::Error::other("s line incomplete"))?
        .to_string();
    components.push(MafComp {
        src,
        start,
        size,
        strand,
        src_size,
        text,
    });
    Ok(())
}

/// Parse a MAF block starting from the given `a` header line, drawing subsequent lines from `iter`.
pub fn parse_maf_block(
    header: String,
    iter: impl Iterator<Item = Result<String, io::Error>>,
) -> Result<MafAli, io::Error> {
    let mut block_lines = vec![];
    block_lines.push(header);

    for line_res in iter {
        let line: String = line_res?;
        if line.is_empty() {
            // Blank lines terminate the "paragraph".
            break;
        }
        block_lines.push(line);
    }
    let mut components: Vec<MafComp> = vec![];
    let mut score: Option<f64> = None;

    for line in block_lines {
        let mut fields: Vec<_> = line.split_whitespace().collect();
        match fields[0] {
            "a" => {
                for f in &fields[1..] {
                    if let Some(val) = f.strip_prefix("score=") {
                        score = val.parse::<f64>().ok();
                    }
                }
            }
            "s" => parse_s_line(&mut fields, &mut components)?,
            "i" => (),
            "e" => (),
            "q" => (),
            "track" => (),
            _ => return Err(io::Error::other("BadLineType")),
        };
    }

    Ok(MafAli { score, components })
}

pub struct MafWriter<W: Write> {
    writer: W,
}

impl<W: Write> MafWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    pub fn write_header(&mut self, program: &str) -> io::Result<()> {
        writeln!(self.writer, "##maf version=1 scoring={}", program)
    }

    pub fn write_ali(&mut self, ali: &MafAli) -> io::Result<()> {
        writeln!(self.writer, "a score={:.1}", ali.score.unwrap_or(0.0))?;
        for comp in &ali.components {
            writeln!(
                self.writer,
                "s {:<20} {:10} {:10} {} {:10} {}",
                comp.src, comp.start, comp.size, comp.strand, comp.src_size, comp.text
            )?;
        }
        writeln!(self.writer)?;
        Ok(())
    }
}

/// Convert a two-sequence MAF block into a PAF record.
///
/// Returns `Ok(None)` if the block does not have exactly two components.
/// Multi-sequence blocks are skipped (caller should log a warning if desired).
pub fn maf_block_to_paf(block: &MafAli) -> anyhow::Result<Option<PafRecord>> {
    if block.components.len() < 2 {
        return Ok(None);
    }
    if block.components.len() > 2 {
        return Ok(None); // caller logs warning
    }

    let ref_entry = &block.components[0];
    let qry_entry = &block.components[1];

    let cigar_ops = cigar_from_alignment(ref_entry.text.as_bytes(), qry_entry.text.as_bytes())?;
    let stats = cigar_stats(&cigar_ops);
    let gi = gap_compressed_identity(&cigar_ops);
    let bi = block_identity(&cigar_ops);
    let cigar_str = format_cigar(&cigar_ops);

    let mut tags = vec![
        format!("gi:f:{gi:.6}"),
        format!("bi:f:{bi:.6}"),
        format!("cg:Z:{cigar_str}"),
    ];
    if let Some(s) = block.score {
        tags.push(format!("ms:i:{}", s as u64));
    }

    let rec = PafRecord {
        query_name: qry_entry.src.clone(),
        query_length: qry_entry.src_size as u32,
        query_start: qry_entry.start as u32,
        query_end: (qry_entry.start + qry_entry.size) as u32,
        strand: qry_entry.strand,
        target_name: ref_entry.src.clone(),
        target_length: ref_entry.src_size as u32,
        target_start: ref_entry.start as u32,
        target_end: (ref_entry.start + ref_entry.size) as u32,
        matches: stats.matches,
        block_length: crate::libs::paf::cigar::block_length(&stats),
        mapq: 255,
        tags,
    };

    Ok(Some(rec))
}

#[cfg(test)]
mod maf_tests {
    use super::*;
    use std::io::{BufRead, BufReader};

    #[test]
    fn parse_comment() {
        let str = "##maf version=1";
        let mut reader = BufReader::new(str.as_bytes());
        let res = next_maf_block(&mut reader);
        eprintln!("got error {:?}", res.as_ref().err());
        assert!(matches!(res.unwrap_err().kind(), io::ErrorKind::Other));
    }

    #[test]
    fn parse_blank_comment() {
        let str = "#";
        let mut reader = BufReader::new(str.as_bytes());
        let res = next_maf_block(&mut reader);
        assert!(matches!(res.unwrap_err().kind(), io::ErrorKind::Other));
    }

    #[test]
    fn parse_err_unexpected() {
        let str = "#\nUnexpected";
        let mut reader = BufReader::new(str.as_bytes());
        let res = next_maf_block(&mut reader);
        eprintln!("got error {:?}", res.as_ref().err());
        assert!(matches!(res.unwrap_err().kind(), io::ErrorKind::Other));
    }

    #[test]
    fn parse_err_s() {
        let str = "#\na\ns 123";
        let mut reader = BufReader::new(str.as_bytes());
        let res = next_maf_block(&mut reader);
        eprintln!("got error {:?}", res.as_ref().err());
        assert!(matches!(res.unwrap_err().kind(), io::ErrorKind::Other));
    }

    #[test]
    fn parse_block_a() {
        let str = "#\na score=23262.0 pass=2";
        let mut reader = BufReader::new(str.as_bytes());
        match next_maf_block(&mut reader) {
            Err(e) => panic!("Got error {:?}", e),
            Ok(val) => assert_eq!(
                val,
                MafAli {
                    components: vec![],
                    score: Some(23262.0)
                }
            ),
        }
    }

    #[test]
    fn parse_block_a_empty() {
        let str = "#\na";
        let mut reader = BufReader::new(str.as_bytes());
        match next_maf_block(&mut reader) {
            Err(e) => panic!("Got error {:?}", e),
            Ok(val) => assert_eq!(
                val,
                MafAli {
                    components: vec![],
                    score: None
                }
            ),
        }
    }

    #[test]
    fn parse_block_s_lines() {
        let str = "a meta1=val1 meta2=val2
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s baboon         249182 12 -   4622798 gcagctgaa-aca
i baboon       I 234 n 19
s mm4.chr6     53310102 12 + 151104725 ACAGCTGA-AATA

this line is a canary to ensure it stops after a 'paragraph'";
        let mut lines = BufReader::new(str.as_bytes()).lines();
        let header = lines.next().unwrap().unwrap();
        match parse_maf_block(header, lines) {
            Err(e) => panic!("got error {:?}", e),
            Ok(val) => assert_eq!(
                val,
                MafAli {
                    components: vec![
                        MafComp {
                            src: "hg16.chr7".to_owned(),
                            start: 27707221,
                            size: 13,
                            strand: '+',
                            src_size: 158545518,
                            text: "gcagctgaaaaca".to_owned(),
                        },
                        MafComp {
                            src: "baboon".to_owned(),
                            start: 249182,
                            size: 12,
                            strand: '-',
                            src_size: 4622798,
                            text: "gcagctgaa-aca".to_owned(),
                        },
                        MafComp {
                            src: "mm4.chr6".to_owned(),
                            start: 53310102,
                            size: 12,
                            strand: '+',
                            src_size: 151104725,
                            text: "ACAGCTGA-AATA".to_owned(),
                        },
                    ],
                    score: None,
                }
            ),
        }
    }

    #[test]
    fn parse_block_s_range() {
        let str = "##maf version=1 scoring=multiz
a score=514600.0
s S288c.VIII          13376 34 + 562643 TTACTCGTCTTGCGGCCAAAACTCGAAGAAAAAC
s RM11_1a.scaffold_12  3529 34 + 536628 TTACTCGTCTTGCGGCCAAAACTCGAAGAAAAAC
s EC1118.FN393072_1    8746 34 + 161280 TTACTCGTCTTGCGGCCAAAACTCGAAGAAAAAC
s Spar.gi_29362578      637 33 -  73522 TTACCCGTCTTGCGTCCAAAACTCGAA-AAAAAC

a score=36468.0
s S288c.VIII          193447  99 + 562643 CG--GCATAATTTTTTCCAGGCACTTTCCGCTGCAG---TTGTTGTGCTGACAATAGTCCCATCTAGGTCAAAAAGACAAAGATCTACTGAAAATTGTGGCAtt
s RM11_1a.scaffold_12 189216 101 + 536628 CGTAACACAACTTGGTCCATGC---TTTCTCTGCGGCCACTGTTGTACTCACTATGGTACCATCTAGGTCAAAAAGACATAGATCAGCTGAAAATTCTGCCATT
s EC1118.FN393073_1    25682  99 +  44323 CG--GCATAATTTTTTCCAGGCACTTTCCGCTGCAG---TTGTTGTGCTGACAATAGTCCCATCTAGGTCAAAAAGACAAAGATCTACTGAAAATTGTGGCAtt
s Spar.gi_29362604    100946  97 - 143114 CG--ACATAGTTTTTTCCAGGCACTTTCAGCTGCGG---TTGTTGTGCTAACAATGGTCCCATCTAGGTCAAAAAGGCAGAGATCTACTGAAAATTGTGGCA--
";
        let mut reader = BufReader::new(str.as_bytes());
        let block = next_maf_block(&mut reader).unwrap();
        assert_eq!(
            block.components.first().unwrap().to_range(),
            "S288c.VIII(+):13377-13410".to_string()
        );
        assert_eq!(
            block.components.get(3).unwrap().to_range(),
            "Spar.gi_29362578(-):72853-72885".to_string()
        );

        let block = next_maf_block(&mut reader).unwrap();
        assert_eq!(
            block.components.get(1).unwrap().to_range(),
            "RM11_1a.scaffold_12(+):189217-189317".to_string()
        );
        assert_eq!(
            block.components.get(3).unwrap().to_range(),
            "Spar.gi_29362604(-):42072-42168".to_string()
        );
    }
}
