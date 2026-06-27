use clap::*;
use pgr::libs::paf::cigar::{
    block_identity, cigar_from_alignment, cigar_stats, format_cigar, gap_compressed_identity,
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

fn build_tags(gi: f64, bi: f64, cigar: &str, score: Option<f64>) -> Vec<String> {
    let mut tags = vec![
        format!("gi:f:{gi:.6}"),
        format!("bi:f:{bi:.6}"),
        format!("cg:Z:{cigar}"),
    ];
    if let Some(s) = score {
        tags.push(format!("ms:i:{}", s as u64));
    }
    tags
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
                eprintln!(
                    "Warning: skipping block with {} sequences (only two-sequence blocks are supported)",
                    block.entries.len()
                );
                continue;
            }

            let ref_entry = &block.entries[0];
            let qry_entry = &block.entries[1];

            let cigar_ops = cigar_from_alignment(&ref_entry.alignment, &qry_entry.alignment);
            let stats = cigar_stats(&cigar_ops);
            let gi = gap_compressed_identity(&cigar_ops);
            let bi = block_identity(&cigar_ops);
            let cigar_str = format_cigar(&cigar_ops);
            let tags = build_tags(gi, bi, &cigar_str, block.score);

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
                matches: stats.matches,
                block_length: pgr::libs::paf::cigar::block_length(&stats),
                mapq: 255,
                tags,
            };

            write_paf_record(&mut writer, &rec)?;
        }
    }

    Ok(())
}
