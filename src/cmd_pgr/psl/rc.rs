use clap::Command;
use intspan::{reader, writer};
use std::io::BufRead;

use pgr::libs::fmt::psl::Psl;

pub fn make_subcommand() -> Command {
    Command::new("rc")
        .about("Reverse-complements psl")
        .after_help(
            r###"
Reverse-complement psl.

Examples:
  pgr psl rc in.psl -o out.psl
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input PSL file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(args: &clap::ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);

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
