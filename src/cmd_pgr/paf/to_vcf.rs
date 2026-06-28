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
outputs a VCF file with substitutions called from a POA multiple
sequence alignment.

For each region, all homologous fragments (target first, then each
query, '-' strand reverse-complemented) are fed into the POA engine to
produce a multi-way MSA. Substitution columns (where target has a
non-gap base and at least one query differs) become VCF records:
REF is the target base, ALT are the distinct non-REF bases, and GT
fields encode each sample's base (0=REF, 1..=ALT index, .=gap).

Recommended with --transitive to gather all homologous fragments of
each region.

-f/--fasta-tsv (required):
  TSV with two columns: genome_name <tab> bgzf_fasta_path
  Each genome_name must match a query/target name in the PAF index.
  All genome names in the PAF index must be present in the TSV.

Notes:
* Substitutions only; indels (gap columns) are skipped
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

// Output VCF records from POA MSA of each region. Substitutions only:
// columns where target (first sequence) has a non-gap ACGT base and at
// least one query differs become VCF rows. Indel columns (target gap or
// query gap) are skipped.
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

        // Walk MSA columns. Skip columns where target has a gap (would be an
        // indel — out of scope for substitution-only VCF). For non-gap
        // target columns, call substitutions: REF = target base, ALT =
        // distinct non-REF non-gap bases among queries.
        let mut t_aln_pos: i32 = 0; // target offset within its aligned seq
        for col in 0..aln_len {
            let t_base = msa[0].as_bytes()[col];
            if t_base == b'-' {
                // Target gap — indel column, skip without advancing t_aln_pos.
                continue;
            }

            let ref_base = t_base.to_ascii_uppercase();

            // Collect distinct ALT bases (non-REF, non-gap ACGT) across all
            // sequences (including target — but target == REF by definition).
            let mut alt_bases: Vec<u8> = Vec::new();
            for seq in msa.iter().take(n_seq) {
                let b = seq.as_bytes()[col].to_ascii_uppercase();
                if matches!(b, b'A' | b'C' | b'G' | b'T')
                    && b != ref_base
                    && !alt_bases.contains(&b)
                {
                    alt_bases.push(b);
                }
            }

            // POS in target coordinates: target.start + t_aln_pos (0-based
            // MAF start is the forward-strand position of the first
            // non-gap base; VCF POS is 1-based).
            let pos = target.start + t_aln_pos + 1;
            t_aln_pos += 1;

            // Only emit a row if there is at least one ALT.
            if alt_bases.is_empty() {
                continue;
            }

            let alt_str: Vec<String> = alt_bases
                .iter()
                .map(|b| String::from_utf8_lossy(&[*b]).to_string())
                .collect();
            let alt = alt_str.join(",");

            let chrom = &target.name;
            let mut row = String::new();
            row.push_str(chrom);
            row.push('\t');
            row.push_str(&pos.to_string());
            row.push_str("\t.\t");
            row.push_str(&String::from_utf8_lossy(&[ref_base]));
            row.push('\t');
            row.push_str(&alt);
            row.push_str("\t.\t.\t.\tGT");

            for seq in msa.iter().take(n_seq) {
                row.push('\t');
                let b = seq.as_bytes()[col].to_ascii_uppercase();
                let gt: String = if !matches!(b, b'A' | b'C' | b'G' | b'T') {
                    ".".to_string()
                } else if b == ref_base {
                    "0".to_string()
                } else {
                    match alt_bases.iter().position(|&x| x == b) {
                        Some(i) => (i + 1).to_string(),
                        None => ".".to_string(),
                    }
                };
                row.push_str(&gt);
            }

            row.push('\n');
            writer.write_all(row.as_bytes())?;
        }
    }

    Ok(())
}
