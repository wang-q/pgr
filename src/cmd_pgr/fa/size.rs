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
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input FASTA file(s) to process"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
        .arg(
            Arg::new("no_ns")
                .long("no-ns")
                .action(ArgAction::SetTrue)
                .help("Output size without Ns"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());
    let no_ns = args.get_flag("no_ns");

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile);
        let mut fa_in = noodles_fasta::io::Reader::new(reader);

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
