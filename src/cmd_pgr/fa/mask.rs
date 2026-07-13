use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for mask.
pub fn make_subcommand() -> Command {
    Command::new("mask")
        .about("Masks regions in FASTA file(s)")
        .after_help(
            r###"
This command masks regions in FASTA files based on a region file (BED/GFF/etc.).

Masking modes:
* Soft-masking (default): Convert to lowercase
* Hard-masking (--hard): Replace with N's

Input format (runlist.json):
{
    "seq1": "1-100,200-300",    # Mask positions 1-100 and 200-300
    "seq2": "50-150",           # Mask positions 50-150
    "seq3": "1-50,90-100,..."   # Multiple regions allowed
}

Notes:
* 1-based coordinates
* Inclusive ranges
* Sequences not in runlist remain unchanged
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Out-of-range spans cause an error (not silently ignored)

Examples:
1. Soft-mask regions:
   pgr fa mask input.fa --runlist regions.json -o output.fa

2. Hard-mask regions:
   pgr fa mask input.fa --runlist regions.json --hard -o output.fa

3. Process gzipped files:
   pgr fa mask input.fa.gz --runlist regions.json -o output.fa.gz

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTA file to process",
        ))
        .arg(crate::cmd_pgr::args::runlist_arg())
        .arg(
            Arg::new("hard")
                .long("hard")
                .action(ArgAction::SetTrue)
                .help("Hard-mask regions (replace with N's)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the mask command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let mut fa_in = pgr::libs::fmt::fa::reader(infile)
        .with_context(|| format!("Failed to open reader for {}", infile))?;

    let runlists = pgr::libs::io::read_runlist(args.get_one::<String>("runlist").unwrap())?;

    let is_hard = args.get_flag("hard");

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut fa_out = pgr::libs::fmt::fa::writer(outfile)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();

        if let Some(ints) = runlists.get(&name) {
            let seq_str = String::from_utf8(seq[..].into())?;
            let seq_out = pgr::libs::fmt::fa::mask_sequence(&seq_str, ints, is_hard)?;
            let record_out =
                pgr::libs::fmt::fa::new_record_preserving_desc(&name, &record, seq_out.as_bytes());
            fa_out.write_record(&record_out)?;
        } else {
            fa_out.write_record(&record)?;
        }
    }

    fa_out.get_mut().flush()?;

    Ok(())
}
