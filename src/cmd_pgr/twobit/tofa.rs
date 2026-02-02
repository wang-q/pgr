use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::twobit::TwoBitFile;
use std::io::Write;

pub fn make_subcommand() -> Command {
    Command::new("tofa")
        .about("Convert 2bit to FASTA")
        .after_help(
            r###"
Examples:
  # Convert entire 2bit file to FASTA
  pgr twobit tofa input.2bit -o output.fa

  # No masking (all uppercase)
  pgr twobit tofa input.2bit --no-mask -o out.fa
"###,
        )
        .arg(
            Arg::new("input")
                .help("Input 2bit file")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output FASTA file")
                .default_value("stdout"),
        )
        .arg(
            Arg::new("no_mask")
                .long("no-mask")
                .action(ArgAction::SetTrue)
                .help("Convert sequence to all upper case"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("input").unwrap();
    let output_path = args.get_one::<String>("output").unwrap();
    let no_mask = args.get_flag("no_mask");

    let mut tb = TwoBitFile::open(input_path)?;
    let mut writer = intspan::writer(output_path);

    let names = tb.get_sequence_names();
    for name in names {
        let seq = tb.read_sequence(&name, None, None, no_mask)?;
        
        writeln!(writer, ">{}", name)?;
        
        let mut idx = 0;
        let len = seq.len();
        while idx < len {
            let next_idx = (idx + 60).min(len);
            writeln!(writer, "{}", &seq[idx..next_idx])?;
            idx = next_idx;
        }
    }

    Ok(())
}
