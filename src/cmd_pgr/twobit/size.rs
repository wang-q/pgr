use clap::*;
use pgr::libs::twobit::TwoBitFile;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("size")
        .about("Get sequence sizes from 2bit file")
        .after_help(
            r###"
This command retrieves the sequence sizes from a 2bit file.

Examples:
1. Get sizes from a 2bit file:
   pgr 2bit size input.2bit

2. Save the output to a file:
   pgr 2bit size input.2bit -o output.tsv

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input 2bit file(s) to process"),
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
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let no_ns = args.get_flag("no_ns");

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut tb = TwoBitFile::open(infile)?;

        // Get all sequence names
        let mut names = tb.get_sequence_names();
        // Sort names to be deterministic (optional but good)
        names.sort();

        for name in names {
            if no_ns {
                let len = tb.get_sequence_len_no_ns(&name)?;
                writer.write_fmt(format_args!("{}\t{}\n", name, len))?;
            } else {
                let len = tb.get_sequence_len(&name)?;
                writer.write_fmt(format_args!("{}\t{}\n", name, len))?;
            }
        }
    }

    Ok(())
}
