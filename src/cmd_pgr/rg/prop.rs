//! `pgr rg prop` — proportion of `.rg` ranges intersecting a runlist.

use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::{BufRead, Write};

/// Build the clap subcommand for prop.
pub fn make_subcommand() -> Command {
    Command::new("prop")
        .about("Proportion of ranges intersecting a runlist file")
        .after_help(
            r###"
For each range in the `.rg` files, appends the proportion of the range covered
by the runlist (intersection size / range length, 4 decimals). With `--full`,
the range length and the intersection size are appended as well. Lines without
a valid range are skipped.

Examples:
1. Intersection proportion against a runlist:
   pgr rg prop intergenic.json a.rg
2. Also append length and intersection size:
   pgr rg prop intergenic.json a.rg --full
"###,
        )
        .arg(
            Arg::new("runlist")
                .required(true)
                .index(1)
                .num_args(1)
                .help("Runlist JSON file to intersect against"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .index(2)
                .num_args(1..)
                .help("Input .rg files to process"),
        )
        .arg(
            Arg::new("full")
                .long("full")
                .action(ArgAction::SetTrue)
                .help("Also append `length` and `size` fields"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the prop command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let is_full = args.get_flag("full");

    let json = pgr::libs::runlist::read_json(args.get_one::<String>("runlist").unwrap())?;
    let set = pgr::libs::runlist::json_to_set(&json)?;

    let mut writer = pgr::writer(outfile)?;
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile)?;
        for line in reader.lines() {
            let line = line?;
            let range = pgr::libs::ds::Range::from_str(&line);
            if !pgr::libs::runlist::usable_range(&range) {
                continue;
            }
            let (prop, length, size) =
                pgr::libs::runlist::range_prop(&set, range.chr(), *range.start(), *range.end());
            if is_full {
                writeln!(writer, "{}\t{:.4}\t{}\t{}", line, prop, length, size)?;
            } else {
                writeln!(writer, "{}\t{:.4}", line, prop)?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}
