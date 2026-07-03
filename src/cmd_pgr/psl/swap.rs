use clap::{Arg, Command};
use intspan::*;
use std::io::BufRead;

use pgr::libs::fmt::psl::Psl;

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
        .arg(crate::cmd_pgr::args::infile_arg().help("Input PSL file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("no_rc")
                .short('n')
                .long("no-rc")
                .action(clap::ArgAction::SetTrue)
                .help("Don't reverse-complement PSL if needed, instead make target strand explict"),
        )
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(matches);
    let output = crate::cmd_pgr::args::get_outfile(matches);
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
