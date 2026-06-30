use clap::*;

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
        .arg(crate::cmd_pgr::args::infiles_arg("MAF"))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        while let Ok(block) = pgr::libs::fmt::maf::next_maf_block(&mut reader) {
            if block.components.len() > 2 {
                log::warn!(
                    "skipping block with {} sequences (only two-sequence blocks are supported)",
                    block.components.len()
                );
                continue;
            }

            if let Some(rec) = pgr::libs::paf::maf_import::maf_block_to_paf(&block)? {
                write_paf_record(&mut writer, &rec)?;
            }
        }
    }

    Ok(())
}
