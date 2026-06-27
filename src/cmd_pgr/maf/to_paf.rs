use clap::*;
use pgr::libs::paf::cigar::{
    block_identity, cigar_from_alignment, format_cigar, gap_compressed_identity,
};
use pgr::libs::paf::record::PafRecord;
use pgr::libs::paf::writer::write_paf_record;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("to-paf")
        .about("Convert two-sequence MAF files to PAF format")
        .after_help(
            r###"
Converts MAF (Multiple Alignment Format) files containing pairwise alignments
into PAF (Pairwise mApping Format).

Only blocks with exactly two `s` lines are converted.  Multi-sequence blocks
are skipped with a warning.

Custom PAF tags:
* `cg:Z:` – CIGAR string derived from the MAF alignment strings
* `gi:f:` – gap-compressed identity
* `bi:f:` – block identity
* `ms:i:` – MAF score (from the `a` line `score=` field)

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Convert a MAF file to PAF:
   pgr maf to-paf ref_vs_query.maf -o ref_vs_query.paf

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input MAF file(s) to process"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile);

        while let Ok(block) = pgr::libs::fmt::fas::next_maf_block(&mut reader) {
            if block.entries.len() < 2 {
                continue;
            }
            if block.entries.len() > 2 {
                // Multi-sequence block — skip with warning
                eprintln!(
                    "Warning: skipping block with {} sequences (only two-sequence blocks are supported)",
                    block.entries.len()
                );
                continue;
            }

            let ref_entry = &block.entries[0];
            let qry_entry = &block.entries[1];

            // Generate CIGAR from alignment vectors
            let cigar_ops = cigar_from_alignment(&ref_entry.alignment, &qry_entry.alignment);
            let cigar_str = format_cigar(&cigar_ops);

            // Match counts
            let (matches, _mismatches, _ins, _del) = count_cigar_stats(&cigar_ops);
            let block_len = qry_entry.alignment.len() as u32;

            // Build PAF record
            let rec = PafRecord {
                query_name: qry_entry.src.clone(),
                query_length: qry_entry.src_size as u32,
                query_start: qry_entry.start as u32,
                query_end: (qry_entry.start + qry_entry.size) as u32,
                strand: qry_entry.strand.chars().next().unwrap_or('+'),
                target_name: ref_entry.src.clone(),
                target_length: ref_entry.src_size as u32,
                target_start: ref_entry.start as u32,
                target_end: (ref_entry.start + ref_entry.size) as u32,
                matches,
                block_length: block_len,
                mapq: 255,
            };

            let gi = gap_compressed_identity(&cigar_ops);
            let bi = block_identity(&cigar_ops);
            let score = block.score.map(|s| s as u64);

            write_paf_record(&mut writer, &rec, gi, bi, &cigar_str, score)?;
        }
    }

    Ok(())
}

/// Count matches, mismatches, insertions, and deletions from CIGAR ops.
///
/// Note: M is counted as matches; X as mismatches; I/D as per their type.
fn count_cigar_stats(ops: &[pgr::libs::paf::cigar::CigarOp]) -> (u32, u32, u32, u32) {
    let mut m = 0u32;
    let mut x = 0u32;
    let mut i = 0u32;
    let mut d = 0u32;
    for op in ops {
        let len = op.len();
        match op.op() {
            'M' | '=' => m += len,
            'X' => x += len,
            'I' => i += len,
            'D' => d += len,
            _ => {}
        }
    }
    (m, x, i, d)
}
