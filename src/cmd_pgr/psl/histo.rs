use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for histo.
pub fn make_subcommand() -> Command {
    Command::new("histo")
        .about("Collects counts on PSL alignments for making histograms")
        .after_help(
            r###"
These can then be analyzed with R, textHistogram, etc.

The 'field' argument determines what data to collect, the following are currently supported:

* alignsPerQuery - number of alignments per query. Output is one line per query with the number of alignments.

* coverSpread - difference between the highest and lowest coverage for alignments of a query. Output line per query, with the difference. Use --multi-only to omit queries with a single alignment.

* idSpread - difference between the highest and lowest fraction identity for alignments of a query. Output line per query, with the difference.

Examples:
1. Collect alignment counts per query:
   pgr psl histo --field alignsPerQuery in.psl -o out.histo
"###,
        )
        .arg(
            Arg::new("field")
                .long("field")
                .required(true)
                .value_parser(["alignsPerQuery", "coverSpread", "idSpread"])
                .help("What data to collect"),
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input PSL file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("multi_only")
                .long("multi-only")
                .action(ArgAction::SetTrue)
                .help("Omit queries with only one alignment"),
        )
        .arg(
            Arg::new("non_zero")
                .long("non-zero")
                .short('z')
                .action(ArgAction::SetTrue)
                .help("Omit queries with zero values"),
        )
}
/// Execute the histo command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let what = args.get_one::<String>("field").unwrap();
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let multi_only = args.get_flag("multi_only");
    let non_zero = args.get_flag("non_zero");

    let reader =
        pgr::reader(input).with_context(|| format!("Failed to open reader for {}", input))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    pgr::libs::fmt::psl::histogram(reader, &mut writer, what, multi_only, non_zero)?;

    writer.flush()?;
    Ok(())
}
