//! `pgr rg count` — count overlaps between `.rg` ranges.

use clap::{Arg, ArgMatches, Command};
use std::io::{BufRead, Write};

/// Build the clap subcommand for count.
pub fn make_subcommand() -> Command {
    Command::new("count")
        .about("Counts overlaps between ranges in a target file and other range files")
        .after_help(
            r###"
Counts, for each range in the target `.rg` file, how many ranges in the other
`.rg` files overlap it, appending the count as an extra tab-separated field.
Lines without a valid range are skipped.

Examples:
1. Count overlaps between two .rg files:
   pgr rg count target.rg intervals.rg
2. Count overlaps with intervals from stdin:
   pgr rg count target.rg stdin
"###,
        )
        .arg(
            Arg::new("target")
                .required(true)
                .index(1)
                .num_args(1)
                .help("Target .rg file to count ranges for"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .index(2)
                .num_args(1..)
                .help("Input .rg files to count overlaps with"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the count command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let files: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let index = pgr::libs::runlist::RgIndex::from_files(&files)?;

    let mut writer = pgr::writer(outfile)?;
    let reader = pgr::reader(args.get_one::<String>("target").unwrap())?;
    for line in reader.lines() {
        let line = line?;
        let range = pgr::libs::ds::Range::from_str(&line);
        if !range.is_valid() || range.start() > range.end() {
            continue;
        }
        let n = index.count(range.chr(), *range.start(), *range.end());
        writeln!(writer, "{}\t{}", line, n)?;
    }
    writer.flush()?;
    Ok(())
}
