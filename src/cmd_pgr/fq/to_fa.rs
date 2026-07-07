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
* Supports compressed input/output
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
    let mut fa_out = pgr::libs::fmt::fa::writer(outfile)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
        let mut seq_in = noodles_fastq::io::Reader::new(reader);

        for result in seq_in.records() {
            // obtain record or fail with error
            let record = result?;

            // Output FASTA format
            let name = std::str::from_utf8(record.name())?;
            let record_out = pgr::libs::fmt::fa::new_record(name, record.sequence());
            fa_out.write_record(&record_out)?;
        }
    }

    fa_out.get_mut().flush()?;

    Ok(())
}
