use clap::*;
use std::io::Write;

use pgr::libs::paf::index::PafIndex;
use pgr::libs::paf::index::QueryResult;

use super::query;
use super::to_maf::{build_msa_entries, load_fasta_tsv, FastaStore};

pub fn make_subcommand() -> Command {
    query::add_query_args(
        Command::new("to-vcf")
            .arg(
                Arg::new("fasta_tsv")
                    .long("fasta-tsv")
                    .short('f')
                    .required(true)
                    .num_args(1)
                    .help("TSV file: genome_name <tab> bgzf_fasta_path"),
            )
            .arg(
                Arg::new("match_score")
                    .long("match")
                    .num_args(1)
                    .default_value("5")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA match score (default: 5)"),
            )
            .arg(
                Arg::new("mismatch_score")
                    .long("mismatch")
                    .num_args(1)
                    .default_value("-4")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA mismatch score (default: -4)"),
            )
            .arg(
                Arg::new("gap_open")
                    .long("gap-open")
                    .num_args(1)
                    .default_value("-8")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA gap open penalty (default: -8)"),
            )
            .arg(
                Arg::new("gap_extend")
                    .long("gap-extend")
                    .num_args(1)
                    .default_value("-6")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA gap extend penalty (default: -6)"),
            )
            .arg(
                Arg::new("outfile")
                    .long("outfile")
                    .short('o')
                    .num_args(1)
                    .default_value("stdout")
                    .help("Output filename. [stdout] for screen"),
            ),
    )
    .about("Query PAF index and output multi-way VCF via POA MSA")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs a VCF file with substitutions and indels called from a POA
multiple sequence alignment.

For each region, all homologous fragments (target first, then each
query, '-' strand reverse-complemented) are fed into the POA engine to
produce a multi-way MSA. Three variant classes are emitted:

* SNP: single target non-gap column where >=1 query differs. REF is
  the target base, ALT are the distinct non-REF bases.
* INS: consecutive target gap columns. REF is the 1bp anchor (target
  base just before the gap), ALT is anchor + inserted bases per sample.
* DEL: consecutive target non-gap columns where >=1 query has a gap.
  REF is the target segment, ALT is the per-sample non-gap concatenation.

GT fields encode each sample's allele (0=REF, 1..=N=ALT index, '.'=gap
or non-ACGT). DEL is not left-aligned -- use `bcftools norm -f ref.fa`
to normalize if needed.

Recommended with --transitive to gather all homologous fragments of
each region.

-f/--fasta-tsv (required):
  TSV with two columns: genome_name <tab> bgzf_fasta_path
  Each genome_name must match a query/target name in the PAF index.
  All genome names in the PAF index must be present in the TSV.

Notes:
* Input PAF files should contain cg:Z: tags (used for query projection)
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'

Examples:
1. Single region to VCF:
   pgr paf to-vcf alignments.paf chr1:1000-5000 -f genomes.tsv

2. Multi-way VCF with transitive BFS:
   pgr paf to-vcf alignments.paf chr1:1000-5000 -t -f genomes.tsv

3. Batch query from BED regions:
   pgr paf to-vcf alignments.paf.idx -b regions.bed -f genomes.tsv

"###,
    )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let tsv_path = args.get_one::<String>("fasta_tsv").unwrap();
    let seq_to_file = load_fasta_tsv(tsv_path)?;

    let (idx, all_results) = query::run_query(args)?;

    // Validate: every name in the PAF index must be present in the TSV.
    let mut missing: Vec<&str> = idx
        .names
        .keys()
        .filter(|n| !seq_to_file.contains_key(*n))
        .map(|n| n.as_str())
        .collect();
    missing.sort_unstable();
    if !missing.is_empty() {
        anyhow::bail!(
            "FASTA TSV is missing {} genome(s) present in PAF index: {}",
            missing.len(),
            missing.join(", ")
        );
    }

    let mut fasta_store = FastaStore::new(&seq_to_file)?;

    let match_score = *args.get_one::<i32>("match_score").unwrap();
    let mismatch_score = *args.get_one::<i32>("mismatch_score").unwrap();
    let gap_open = *args.get_one::<i32>("gap_open").unwrap();
    let gap_extend = *args.get_one::<i32>("gap_extend").unwrap();

    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());

    output_vcf(
        &mut writer,
        &idx,
        &all_results,
        &mut fasta_store,
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    )?;

    writer.flush()?;
    Ok(())
}

// Emit one VCF row. ref_allele is REF; alt_alleles are distinct non-REF
// alleles (joined by ','); sample_alleles[i] is sample i's allele string
// (empty or non-ACGT -> GT='.'). GT: 0=REF, 1..=N=ALT index, '.'=other.
fn emit_vcf_row<W: Write>(
    writer: &mut W,
    chrom: &str,
    pos: i32,
    ref_allele: &str,
    alt_alleles: &[String],
    sample_alleles: &[String],
) -> anyhow::Result<()> {
    let alt = alt_alleles.join(",");
    let mut row = String::new();
    row.push_str(chrom);
    row.push('\t');
    row.push_str(&pos.to_string());
    row.push_str("\t.\t");
    row.push_str(ref_allele);
    row.push('\t');
    row.push_str(&alt);
    row.push_str("\t.\t.\t.\tGT");
    for allele in sample_alleles {
        row.push('\t');
        let gt = if allele.is_empty() || allele == "-" {
            ".".to_string()
        } else if allele == ref_allele {
            "0".to_string()
        } else {
            match alt_alleles.iter().position(|a| a == allele) {
                Some(i) => (i + 1).to_string(),
                None => ".".to_string(),
            }
        };
        row.push_str(&gt);
    }
    row.push('\n');
    writer.write_all(row.as_bytes())?;
    Ok(())
}

// Output VCF records from POA MSA of each region. Three variant classes
// are emitted: substitutions (single target non-gap column with ≥1 differing
// query), INS (consecutive target gap columns; REF = 1bp anchor at the
// preceding non-gap column, ALT = anchor + inserted bases), and DEL
// (consecutive target non-gap columns where ≥1 query has gap; REF = anchor
// + target segment, ALT = anchor + per-query non-gap bases; a fully-deleted
// sample gets ALT = anchor). Neither INS nor DEL is left-aligned — use
// `bcftools norm -f ref.fa` to normalize if needed.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn output_vcf<W: Write>(
    writer: &mut W,
    idx: &PafIndex,
    all_results: &[((String, i32, i32), Vec<QueryResult>)],
    fasta_store: &mut FastaStore,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> anyhow::Result<()> {
    let params = pgr::libs::poa::AlignmentParams {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    };

    let mut header_written = false;

    for ((tname_region, _, _), results) in all_results {
        if results.is_empty() {
            continue;
        }

        let entries = build_msa_entries(idx, tname_region, results, fasta_store)?;

        // Run POA MSA.
        let mut poa =
            pgr::libs::poa::Poa::new(params.clone(), pgr::libs::poa::AlignmentType::Global);
        for e in &entries {
            poa.add_sequence(&e.seq);
        }
        let msa = poa.msa();

        if msa.is_empty() {
            continue;
        }

        // Write VCF header once, using the entry names as sample columns.
        if !header_written {
            writer.write_all(b"##fileformat=VCFv4.2\n")?;
            writer
                .write_all(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")?;
            let mut header = String::from("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
            for e in &entries {
                header.push('\t');
                header.push_str(&e.name);
            }
            header.push('\n');
            writer.write_all(header.as_bytes())?;
            header_written = true;
        }

        let n_seq = entries.len();
        let aln_len = msa[0].len();
        let target = &entries[0];

        // Walk MSA columns with a while loop so we can advance past indel
        // regions. t_aln_pos counts target non-gap columns processed (used
        // to derive VCF POS from target.start). Three cases:
        //   1. INS: consecutive target gap columns (REF = 1bp anchor at the
        //      preceding non-gap column, ALT = anchor + inserted bases).
        //   2. DEL: consecutive target non-gap columns where ≥1 query has
        //      gap (REF = target segment, ALT = per-query non-gap concat).
        //   3. SNP: single target non-gap column with no gaps.
        let mut col: usize = 0;
        let mut t_aln_pos: i32 = 0;
        while col < aln_len {
            let t_base = msa[0].as_bytes()[col];
            if t_base == b'-' {
                // INS region: collect consecutive target gap columns.
                let col_start = col;
                while col < aln_len && msa[0].as_bytes()[col] == b'-' {
                    col += 1;
                }
                let col_end = col;
                // Anchor = previous non-gap target column. Skip if none.
                if col_start == 0 {
                    continue;
                }
                let anchor_byte = msa[0].as_bytes()[col_start - 1];
                if anchor_byte == b'-' {
                    continue;
                }
                let anchor = String::from(anchor_byte.to_ascii_uppercase() as char);
                let ref_allele = anchor.clone();
                // Per-sample allele: anchor + inserted bases (gaps dropped).
                let sample_alleles: Vec<String> = msa
                    .iter()
                    .take(n_seq)
                    .map(|seq| {
                        let mut s = anchor.clone();
                        for c in col_start..col_end {
                            let b = seq.as_bytes()[c].to_ascii_uppercase();
                            if matches!(b, b'A' | b'C' | b'G' | b'T') {
                                s.push(b as char);
                            }
                        }
                        s
                    })
                    .collect();
                let mut alt_alleles: Vec<String> = Vec::new();
                for a in &sample_alleles {
                    if a != &ref_allele && !alt_alleles.contains(a) {
                        alt_alleles.push(a.clone());
                    }
                }
                if alt_alleles.is_empty() {
                    continue;
                }
                // POS: anchor column's 1-based target coordinate. t_aln_pos
                // has already counted the anchor column (col_start - 1 is
                // target non-gap), so anchor's 0-based offset = t_aln_pos - 1,
                // and 1-based POS = target.start + t_aln_pos.
                let pos = target.start + t_aln_pos;
                emit_vcf_row(
                    writer,
                    &target.name,
                    pos,
                    &ref_allele,
                    &alt_alleles,
                    &sample_alleles,
                )?;
            } else {
                // target non-gap: check if any query has a gap here.
                let col_has_gap = msa.iter().take(n_seq).any(|s| s.as_bytes()[col] == b'-');
                if col_has_gap {
                    // DEL region: collect consecutive target non-gap columns
                    // where ≥1 query has a gap.
                    let col_start = col;
                    while col < aln_len {
                        let tb = msa[0].as_bytes()[col];
                        if tb == b'-' {
                            break;
                        }
                        let cg = msa.iter().take(n_seq).any(|s| s.as_bytes()[col] == b'-');
                        if cg {
                            col += 1;
                        } else {
                            break;
                        }
                    }
                    let col_end = col;
                    // Anchor = previous non-gap target column. Skip if none
                    // (can't represent a deletion without a 1bp anchor in VCF).
                    if col_start == 0 {
                        t_aln_pos += (col_end - col_start) as i32;
                        continue;
                    }
                    let anchor_byte = msa[0].as_bytes()[col_start - 1];
                    if anchor_byte == b'-' {
                        t_aln_pos += (col_end - col_start) as i32;
                        continue;
                    }
                    let anchor = String::from(anchor_byte.to_ascii_uppercase() as char);
                    // REF = anchor + target segment.
                    let mut ref_allele = anchor.clone();
                    for c in col_start..col_end {
                        let b = msa[0].as_bytes()[c].to_ascii_uppercase();
                        if matches!(b, b'A' | b'C' | b'G' | b'T') {
                            ref_allele.push(b as char);
                        }
                    }
                    // Per-sample allele: anchor + non-gap bases in region.
                    // A sample with all gaps -> allele = anchor (the deletion ALT).
                    let sample_alleles: Vec<String> = msa
                        .iter()
                        .take(n_seq)
                        .map(|seq| {
                            let mut s = anchor.clone();
                            for c in col_start..col_end {
                                let b = seq.as_bytes()[c].to_ascii_uppercase();
                                if matches!(b, b'A' | b'C' | b'G' | b'T') {
                                    s.push(b as char);
                                }
                            }
                            s
                        })
                        .collect();
                    let mut alt_alleles: Vec<String> = Vec::new();
                    for a in &sample_alleles {
                        if a != &ref_allele && !alt_alleles.contains(a) {
                            alt_alleles.push(a.clone());
                        }
                    }
                    // POS: anchor column's 1-based target coordinate.
                    let pos = target.start + t_aln_pos;
                    t_aln_pos += (col_end - col_start) as i32;
                    if alt_alleles.is_empty() {
                        continue;
                    }
                    emit_vcf_row(
                        writer,
                        &target.name,
                        pos,
                        &ref_allele,
                        &alt_alleles,
                        &sample_alleles,
                    )?;
                } else {
                    // SNP: single target non-gap column, no gaps.
                    let ref_base = t_base.to_ascii_uppercase();
                    let ref_allele = String::from(ref_base as char);
                    let sample_alleles: Vec<String> = msa
                        .iter()
                        .take(n_seq)
                        .map(|seq| {
                            let b = seq.as_bytes()[col].to_ascii_uppercase();
                            if matches!(b, b'A' | b'C' | b'G' | b'T') {
                                String::from(b as char)
                            } else {
                                String::new()
                            }
                        })
                        .collect();
                    let mut alt_alleles: Vec<String> = Vec::new();
                    for a in &sample_alleles {
                        if !a.is_empty() && a != &ref_allele && !alt_alleles.contains(a) {
                            alt_alleles.push(a.clone());
                        }
                    }
                    let pos = target.start + t_aln_pos + 1;
                    t_aln_pos += 1;
                    col += 1;
                    if alt_alleles.is_empty() {
                        continue;
                    }
                    emit_vcf_row(
                        writer,
                        &target.name,
                        pos,
                        &ref_allele,
                        &alt_alleles,
                        &sample_alleles,
                    )?;
                }
            }
        }
    }

    Ok(())
}
