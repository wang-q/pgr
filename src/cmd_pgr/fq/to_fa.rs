use clap::{ArgMatches, Command};

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
    let mut fa_out = pgr::libs::fmt::fa::writer(crate::cmd_pgr::args::get_outfile(args))?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile)?;
        let mut seq_in = noodles_fastq::io::Reader::new(reader);

        for result in seq_in.records() {
            // obtain record or fail with error
            let record = result?;

            // Output FASTA format
            let name = String::from_utf8(record.name().to_vec())?;
            let record_out = pgr::libs::fmt::fa::new_record(&name, record.sequence());
            fa_out.write_record(&record_out)?;
        }
    }

    Ok(())
}
