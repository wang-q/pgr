//! `pgr rg runlist` — filter `.rg` lines by comparing with a runlist.

use clap::{builder, Arg, ArgAction, ArgMatches, Command};
use std::io::{BufRead, Write};

/// Build the clap subcommand for runlist.
pub fn make_subcommand() -> Command {
    Command::new("runlist")
        .about("Filters .rg lines by comparing with a runlist file")
        .after_help(
            r###"
Keeps `.rg` lines whose range overlaps, does not overlap, or is fully contained
by the runlist, according to `--op`. `--op superset` keeps only lines whose
range is entirely inside the runlist (the runlist is a superset of the range).
Lines without a valid range are skipped.

Examples:
1. Keep lines overlapping the runlist:
   pgr rg runlist intergenic.json a.rg
2. Keep lines outside the runlist:
   pgr rg runlist intergenic.json a.rg --op non-overlap
3. Keep lines fully contained in the runlist:
   pgr rg runlist intergenic.json a.rg --op superset
"###,
        )
        .arg(
            Arg::new("runlist")
                .required(true)
                .index(1)
                .num_args(1)
                .help("Runlist JSON file to compare against"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .index(2)
                .num_args(1..)
                .help("Input .rg files to process"),
        )
        .arg(
            Arg::new("op")
                .long("op")
                .num_args(1)
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("overlap"),
                    builder::PossibleValue::new("non-overlap"),
                    builder::PossibleValue::new("superset"),
                ])
                .default_value("overlap")
                .help("Filter operation: overlap, non-overlap, or superset (fully contained)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the runlist command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let op = args.get_one::<String>("op").unwrap().as_str();

    let json = pgr::libs::runlist::read_json(args.get_one::<String>("runlist").unwrap())?;
    let set = pgr::libs::runlist::json_to_set(&json)?;
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        args.get_many::<String>("infiles")
            .unwrap()
            .map(String::as_str),
    )?;

    // Validate that every input can be opened before creating the output,
    // so a missing or unreadable input fails without truncating an existing
    // output file. Files are re-opened for streaming below so the number of
    // open descriptors stays bounded regardless of the input count.
    for infile in args.get_many::<String>("infiles").unwrap() {
        pgr::reader(infile)?;
    }
    let mut writer = pgr::writer(outfile)?;
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile)?;
        for line in reader.lines() {
            let line = line?;
            if line.trim_start().starts_with('#') {
                continue;
            }
            let range = pgr::libs::ds::Range::from_str(&line);
            if !pgr::libs::runlist::usable_range(&range) {
                continue;
            }
            // `IntSpan::covered` is O(log n + k) per line; all three ops
            // derive from the covered size.
            let start = *range.start();
            let end = *range.end();
            let size = set.get(range.chr()).map_or(0, |s| s.covered(start, end));
            let length = end - start + 1;
            let keep = match op {
                "overlap" => size > 0,
                "non-overlap" => size == 0,
                // Fully contained: every base of the range is covered.
                "superset" => size == length,
                _ => unreachable!("invalid runlist op"),
            };
            if keep {
                writeln!(writer, "{}", line)?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}
