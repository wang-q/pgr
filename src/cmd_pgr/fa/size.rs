use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("size")
        .about("Counts total bases in FASTA file(s)")
        .after_help(
            r###"
This command counts the total number of bases in one or more FASTA files. It outputs the sequence name
and its length in a tab-separated format.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Count bases in a single FASTA file:
   pgr fa size input.fa

2. Count bases in multiple FASTA files:
   pgr fa size input1.fa input2.fa

3. Save the output to a file:
   pgr fa size input.fa -o output.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("FASTA"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("no_ns")
                .long("no-ns")
                .action(ArgAction::SetTrue)
                .help("Output size without Ns"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let no_ns = args.get_flag("no_ns");

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut fa_in = pgr::libs::fmt::fa::reader(infile)?;

        for result in fa_in.records() {
            let record = result?;
            let name = String::from_utf8(record.name().into())?;
            let seq = record.sequence();

            let len = if no_ns {
                seq.get(..)
                    .unwrap()
                    .iter()
                    .filter(|&&b| !pgr::libs::nt::is_n(b))
                    .count()
            } else {
                seq.len()
            };

            writer.write_fmt(format_args!("{}\t{}\n", name, len))?;
        }
    }

    Ok(())
}
