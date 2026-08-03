//! `pgr runlist compare` — set operations between runlist JSON files.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for compare.
pub fn make_subcommand() -> Command {
    Command::new("compare")
        .about("Compares runlist JSON files")
        .after_help(
            r###"
Applies an operation between the first file (which may hold multiple runlist
sets) and each of the other files. Missing chromosomes are treated as empty.

Examples:
1. Intersection of several runlists:
   pgr runlist compare a.json b.json c.json --op intersect -o out.json
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Sets the input file to use"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .index(2)
                .num_args(1..)
                .help("Sets the input file(s) to use"),
        )
        .arg(
            Arg::new("op")
                .long("op")
                .num_args(1)
                .default_value("intersect")
                .value_parser(["intersect", "union", "diff", "xor"])
                .help("Operations: intersect, union, diff or xor"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the compare command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let op = match args.get_one::<String>("op").unwrap().as_str() {
        "intersect" => pgr::libs::runlist::CompareOp::Intersect,
        "union" => pgr::libs::runlist::CompareOp::Union,
        "diff" => pgr::libs::runlist::CompareOp::Diff,
        "xor" => pgr::libs::runlist::CompareOp::Xor,
        _ => unreachable!("invalid compare op"),
    };

    let first = pgr::libs::runlist::read_json(args.get_one::<String>("infile").unwrap())?;
    let first = pgr::libs::runlist::json_to_sets(&first)?;
    let mut others = Vec::new();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let json = pgr::libs::runlist::read_json(infile)?;
        others.push(pgr::libs::runlist::json_to_set(&json)?);
    }
    let res = pgr::libs::runlist::compare_sets(&first, &others, op);
    pgr::libs::runlist::write_sets(outfile, &res)
}
