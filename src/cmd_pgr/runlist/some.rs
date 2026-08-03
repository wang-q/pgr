//! `pgr runlist some` — extract selected keys from a runlist JSON.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for some.
pub fn make_subcommand() -> Command {
    Command::new("some")
        .about("Extracts some records from a runlist JSON file")
        .after_help(
            r###"
Keeps only the top-level keys listed in the names file (one per line).

Examples:
1. Extract chromosomes of interest:
   pgr runlist some in.json names.txt -o out.json
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Sets the input file to use"),
        )
        .arg(
            Arg::new("list")
                .required(true)
                .num_args(1)
                .index(2)
                .help("File of names to keep, one per line"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the some command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let json = pgr::libs::runlist::read_json(args.get_one::<String>("infile").unwrap())?;
    let names: std::collections::BTreeSet<String> =
        pgr::libs::io::read_names(args.get_one::<String>("list").unwrap())?;
    let out = pgr::libs::runlist::some_json(&json, &names);
    pgr::libs::ds::intspan::write_json(outfile, &out)?;
    Ok(())
}
