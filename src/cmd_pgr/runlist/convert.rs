//! `pgr runlist convert` — convert runlist JSON to `.rg` range lines.

use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for convert.
pub fn make_subcommand() -> Command {
    Command::new("convert")
        .about("Converts runlist JSON files to .rg range lines")
        .after_help(
            r###"
Writes `chr:start-end` lines (one per span) for every chromosome of each
input runlist. With `--longest` only the longest span per chromosome is kept.

Examples:
1. Convert to ranges:
   pgr runlist convert in.json -o out.rg
2. Keep only the longest span per chromosome:
   pgr runlist convert in.json --longest -o out.rg
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("longest")
                .long("longest")
                .action(ArgAction::SetTrue)
                .help("Only keep the longest range"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the convert command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let is_longest = args.get_flag("longest");
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        args.get_many::<String>("infiles")
            .unwrap()
            .map(String::as_str),
    )?;
    let mut writer = pgr::writer(outfile)?;
    for infile in args.get_many::<String>("infiles").unwrap() {
        let json = pgr::libs::runlist::read_json(infile)?;
        let set_of = pgr::libs::runlist::json_to_sets(&json)?;
        for set in set_of.values() {
            for line in pgr::libs::runlist::convert_set(set, is_longest) {
                writeln!(writer, "{}", line)?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}
