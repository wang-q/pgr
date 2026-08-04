//! `pgr runlist combine` — combine multiple sets of a multi runlist JSON.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for combine.
pub fn make_subcommand() -> Command {
    Command::new("combine")
        .about("Combines multiple sets of runlists in a JSON file")
        .after_help(
            r###"
Combines all sets of a multi runlist JSON into one, applying the operation
between the first set and each subsequent set.

Examples:
1. Union of all sets:
   pgr runlist combine in.json -o out.json
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Sets the input file to use"),
        )
        .arg(
            Arg::new("op")
                .long("op")
                .num_args(1)
                .default_value("union")
                .value_parser(["intersect", "union", "diff", "xor"])
                .help("Operations: intersect, union, diff or xor"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the combine command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let op = match args.get_one::<String>("op").unwrap().as_str() {
        "intersect" => pgr::libs::runlist::CompareOp::Intersect,
        "union" => pgr::libs::runlist::CompareOp::Union,
        "diff" => pgr::libs::runlist::CompareOp::Diff,
        "xor" => pgr::libs::runlist::CompareOp::Xor,
        _ => unreachable!("invalid combine op"),
    };
    // The output is a runlist JSON, not the input runlist; refuse to
    // overwrite it.
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        std::iter::once(args.get_one::<String>("infile").unwrap().as_str()),
    )?;
    let json = pgr::libs::runlist::read_json(args.get_one::<String>("infile").unwrap())?;
    let set_of = pgr::libs::runlist::json_to_sets(&json)?;
    let res = pgr::libs::runlist::combine_sets(&set_of, op);
    let json = pgr::libs::ds::intspan::set2json(&res);
    pgr::libs::ds::intspan::write_json(outfile, &json)?;
    Ok(())
}
