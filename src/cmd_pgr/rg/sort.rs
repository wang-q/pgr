//! `pgr rg sort` — sort `.rg` lines by chromosome, start and strand.

use clap::{Arg, ArgMatches, Command};
use std::io::{BufRead, Write};

/// Build the clap subcommand for sort.
pub fn make_subcommand() -> Command {
    Command::new("sort")
        .about("Sorts .rg lines by chromosome, start and strand")
        .after_help(
            r###"
Sorts `.rg` lines by the parsed (chromosome, start, strand) key. Lines without
a valid range are written to the end of the output in their original order.
Lines starting with `#` are treated as comments and skipped. Lines with
identical keys keep their input order (stable sort).

Examples:
1. Sort a .rg file:
   pgr rg sort a.rg
2. Sort multiple files:
   pgr rg sort a.rg b.rg -o sorted.rg
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input .rg file(s) to process"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the sort command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut rows: Vec<(String, pgr::libs::ds::Range)> = Vec::new();
    let mut invalids: Vec<String> = Vec::new();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile)?;
        for line in reader.lines() {
            let line = line?;
            if line.trim_start().starts_with('#') {
                continue;
            }
            let range = pgr::libs::ds::Range::from_str(&line);
            if pgr::libs::runlist::usable_range(&range) {
                rows.push((line, range));
            } else {
                invalids.push(line);
            }
        }
    }

    rows.sort_by(|(_, a), (_, b)| {
        (a.chr(), a.start(), a.strand()).cmp(&(b.chr(), b.start(), b.strand()))
    });

    let mut writer = pgr::writer(outfile)?;
    for (line, _) in &rows {
        writeln!(writer, "{}", line)?;
    }
    for line in &invalids {
        writeln!(writer, "{}", line)?;
    }
    writer.flush()?;
    Ok(())
}
