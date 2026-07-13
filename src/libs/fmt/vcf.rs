//! VCF (Variant Call Format) header and record writers shared by commands
//! that emit SNP-only VCF output.

use std::collections::BTreeMap;
use std::io::Write;

/// Write the VCF header: `##fileformat` first, optional `##contig` lines
/// (when `contigs` is Some), `##FORMAT=GT`, and the `#CHROM` line with `samples`.
pub fn write_vcf_header<W: Write>(
    writer: &mut W,
    contigs: Option<&BTreeMap<String, i32>>,
    samples: &[String],
) -> anyhow::Result<()> {
    writer.write_all(b"##fileformat=VCFv4.2\n")?;
    if let Some(contigs) = contigs {
        for (chr, len) in contigs.iter() {
            writer.write_all(format!("##contig=<ID={},length={}>\n", chr, len).as_ref())?;
        }
    }
    writer.write_all(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")?;
    let mut header = String::from("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for name in samples {
        header.push('\t');
        header.push_str(name);
    }
    header.push('\n');
    writer.write_all(header.as_ref())?;
    Ok(())
}

/// Write a single SNP VCF row.
///
/// * `chrom` — chromosome name.
/// * `pos` — 1-based position.
/// * `ref_base` — reference allele (uppercase).
/// * `alt_bases` — deduplicated list of alternate alleles (uppercase, no ref).
/// * `sample_bases` — one base per sample (uppercase); `0` = ref, `1..N` =
///   index into `alt_bases`+1, `.` = non-ACGT.
pub fn write_snp_row<W: Write>(
    writer: &mut W,
    chrom: &str,
    pos: i32,
    ref_base: char,
    alt_bases: &[char],
    sample_bases: &[u8],
) -> anyhow::Result<()> {
    use itertools::Itertools;

    let alt_str = if alt_bases.is_empty() {
        ".".to_string()
    } else {
        alt_bases.iter().map(|c| c.to_string()).join(",")
    };

    let mut row = String::new();
    row.push_str(chrom);
    row.push('\t');
    row.push_str(&pos.to_string());
    row.push_str("\t.\t");
    row.push(ref_base);
    row.push('\t');
    row.push_str(&alt_str);
    row.push_str("\t.\t.\t.\tGT");

    for &b in sample_bases {
        row.push('\t');
        let c = char::from(b).to_ascii_uppercase();
        let gt = if !matches!(c, 'A' | 'C' | 'G' | 'T') {
            ".".to_string()
        } else if c == ref_base {
            "0".to_string()
        } else {
            match alt_bases.iter().position(|&x| x == c) {
                Some(idx) => (idx + 1).to_string(),
                None => ".".to_string(),
            }
        };
        row.push_str(&gt);
    }

    row.push('\n');
    writer.write_all(row.as_ref())?;
    Ok(())
}
