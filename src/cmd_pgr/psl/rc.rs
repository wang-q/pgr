use clap::{Arg, Command};
use intspan::*;
use std::io::BufRead;

use pgr::libs::psl::Psl;

pub fn make_subcommand() -> Command {
    Command::new("rc")
        .about("Reverse-complement psl")
        .after_help(
            r###"
Reverse-complement psl.

Examples:
  pgr psl rc in.psl -o out.psl
"###,
        )
        .arg(
            Arg::new("input")
                .help("Input PSL file")
                .default_value("stdin")
                .index(1),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Output PSL file")
                .default_value("stdout"),
        )
}

pub fn execute(args: &clap::ArgMatches) -> anyhow::Result<()> {
    let input = args.get_one::<String>("input").unwrap();
    let output = args.get_one::<String>("output").unwrap();

    let reader = reader(input);
    let mut writer = writer(output);

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let mut psl: Psl = line.parse()?;
        psl.rc();
        psl.write_to(&mut writer)?;
    }

    Ok(())
}
