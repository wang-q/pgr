use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use pgr::libs::fmt::twobit::TwoBitFile;

pub fn make_subcommand() -> Command {
    Command::new("to-fa")
        .about("Converts 2bit to FASTA")
        .after_help(
            r###"
This command converts 2bit files to FASTA format.

Notes:
* 2bit files are binary and require random access (seeking)
* Does not support stdin or gzipped inputs

Examples:
  # Convert entire 2bit file to FASTA
  pgr 2bit to-fa input.2bit -o output.fa

  # No masking (all uppercase)
  pgr 2bit to-fa input.2bit --no-mask -o out.fa

  # Set line width (default 60)
  pgr 2bit to-fa input.2bit -l 80
  pgr 2bit to-fa input.2bit -l 0  # no wrapping
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input 2bit file to process",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(crate::cmd_pgr::args::line_arg(Some("60")))
        .arg(
            Arg::new("no_mask")
                .long("no-mask")
                .action(ArgAction::SetTrue)
                .help("Convert sequence to all upper case"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("infile").unwrap();
    let output_path = crate::cmd_pgr::args::get_outfile(args);
    let no_mask = args.get_flag("no_mask");
    let line_width = args.get_one::<usize>("line").copied().unwrap_or(60);

    let mut tb = TwoBitFile::open(input_path)?;
    let mut writer = if line_width == 0 {
        pgr::libs::fmt::fa::writer(output_path)?
    } else {
        pgr::libs::fmt::fa::writer_with_wrap(output_path, line_width)?
    };

    let names = tb.get_sequence_names();
    for name in names {
        let seq = tb.read_sequence(&name, None, None, no_mask)?;
        let record = pgr::libs::fmt::fa::new_record(&name, seq.as_bytes());
        writer.write_record(&record)?;
    }

    Ok(())
}
