//! `pgr rg runlist` — filter `.rg` lines by comparing with a runlist.

use clap::{builder, Arg, ArgAction, ArgMatches, Command};
use std::io::{BufRead, Write};

/// Build the clap subcommand for runlist.
pub fn make_subcommand() -> Command {
    Command::new("runlist")
        .about("Filters .rg lines by comparing with a runlist file")
        .after_help(
            r###"
Keeps `.rg` lines whose range overlaps, does not overlap, or is contained by
the runlist, according to `--op`. Lines without a valid range are skipped.

Examples:
1. Keep lines overlapping the runlist:
   pgr rg runlist intergenic.json a.rg
2. Keep lines outside the runlist:
   pgr rg runlist intergenic.json a.rg --op non-overlap
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
                .help("Filter operation: overlap, non-overlap or superset"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the runlist command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let op = args.get_one::<String>("op").unwrap().as_str();

    let json = pgr::libs::runlist::read_json(args.get_one::<String>("runlist").unwrap())?;
    let set = pgr::libs::runlist::json_to_set(&json)?;
    let index = pgr::libs::runlist::SpanIndex::from_set(&set);

    let mut writer = pgr::writer(outfile)?;
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile)?;
        for line in reader.lines() {
            let line = line?;
            let range = pgr::libs::ds::Range::from_str(&line);
            if !pgr::libs::runlist::usable_range(&range) {
                continue;
            }
            // `overlap` is O(log n + k) per line via binary search over the
            // sorted spans; all three ops derive from the covered size.
            let (size, length) = index.overlap(range.chr(), *range.start(), *range.end());
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
