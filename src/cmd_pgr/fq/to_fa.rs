use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for to-fa.
pub fn make_subcommand() -> Command {
    Command::new("to-fa")
        .about("Converts FASTQ to FASTA format")
        .after_help(
            r###"
This command converts FASTQ format sequences to FASTA format.

Features:
* Automatic format detection
* Preserves sequence names
* Supports compressed input
* Processes multiple input files

Examples:
1. Convert a FASTQ file to FASTA:
   pgr fq to-fa input.fq -o output.fa

2. Convert multiple FASTQ files to a single FASTA:
   pgr fq to-fa input1.fq input2.fq -o output.fa

3. Convert and write to stdout:
   pgr fq to-fa input.fq
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("FASTQ"))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the to-fa command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        args.get_many::<String>("infiles")
            .unwrap()
            .map(|s| s.as_str()),
    )?;
    let mut fa_out = pgr::libs::fmt::fa::writer(outfile)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut seq_in = pgr::libs::fmt::seq::SeqReader::new(infile)
            .with_context(|| format!("Failed to open reader for {}", infile))?;
        let mut rec = pgr::libs::fmt::seq::SeqRecord::new();
        while seq_in.read_record(&mut rec)? {
            // Output FASTA format
            let name = std::str::from_utf8(rec.name())?;
            let record_out = pgr::libs::fmt::fa::new_record(name, rec.sequence());
            fa_out.write_record(&record_out)?;
        }
    }

    fa_out.get_mut().flush()?;

    Ok(())
}
