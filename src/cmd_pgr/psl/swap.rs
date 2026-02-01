use clap::{Arg, Command};
use intspan::*;
use std::io::BufRead;

use pgr::libs::psl::Psl;

pub fn make_subcommand() -> Command {
    Command::new("swap")
        .about("Reverse target and query in psls")
        .after_help(
            r###"
Reverse target and query in psls.

Examples:
  pgr psl swap in.psl -o out.psl
  pgr psl swap --no-rc in.psl -o out.psl
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
        .arg(
            Arg::new("no_rc")
                .short('n')
                .long("no-rc")
                .action(clap::ArgAction::SetTrue)
                .help("Don't reverse-complement PSL if needed, instead make target strand explict"),
        )
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    let input = matches.get_one::<String>("input").unwrap();
    let output = matches.get_one::<String>("output").unwrap();
    let no_rc = matches.get_flag("no_rc");

    let reader = reader(input);
    let mut writer = writer(output);

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let mut psl: Psl = line.parse()?;
        psl.swap(no_rc);
        psl.write_to(&mut writer)?;
    }

    Ok(())
}
