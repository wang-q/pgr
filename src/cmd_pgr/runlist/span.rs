//! `pgr runlist span` — per-span operations on runlist JSON.

use clap::{builder, value_parser, Arg, ArgAction, ArgMatches, Command};

/// Build the clap subcommand for span.
pub fn make_subcommand() -> Command {
    Command::new("span")
        .about("Operates on spans in runlist JSON")
        .after_help(
            r###"
Applies an operation to every chromosome of a runlist JSON (single or multi):

* cover:  a single span from min to max
* holes:  all the holes in the runlist
* trim:   remove N integers from each end of each span
* pad:    add N integers to each end of each span
* excise: remove all spans smaller than N
* fill:   fill in all holes smaller than or equal to N

Examples:
1. Fill holes up to 10 bp:
   pgr runlist span in.json --op fill -n 10 -o out.json
2. Drop fragments shorter than 100 bp:
   pgr runlist span in.json --op excise -n 100 -o out.json
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Runlist JSON file, [stdin] for standard input"),
        )
        .arg(
            Arg::new("op")
                .long("op")
                .num_args(1)
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("cover"),
                    builder::PossibleValue::new("holes"),
                    builder::PossibleValue::new("trim"),
                    builder::PossibleValue::new("pad"),
                    builder::PossibleValue::new("excise"),
                    builder::PossibleValue::new("fill"),
                ])
                .default_value("cover")
                .help("Operations"),
        )
        .arg(
            Arg::new("number")
                .long("number")
                .short('n')
                .num_args(1)
                .value_parser(value_parser!(i32))
                .default_value("0"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the span command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let op = match args.get_one::<String>("op").unwrap().as_str() {
        "cover" => pgr::libs::runlist::SpanOp::Cover,
        "holes" => pgr::libs::runlist::SpanOp::Holes,
        "trim" => pgr::libs::runlist::SpanOp::Trim,
        "pad" => pgr::libs::runlist::SpanOp::Pad,
        "excise" => pgr::libs::runlist::SpanOp::Excise,
        "fill" => pgr::libs::runlist::SpanOp::Fill,
        _ => unreachable!("invalid span op"),
    };
    let number = *args.get_one::<i32>("number").unwrap();

    let json = pgr::libs::runlist::read_json(infile)?;
    let set_of = pgr::libs::runlist::json_to_sets(&json);
    let res: std::collections::BTreeMap<_, _> = set_of
        .iter()
        .map(|(name, set)| (name.clone(), pgr::libs::runlist::span_op(set, op, number)))
        .collect();
    pgr::libs::runlist::write_sets(outfile, &res)
}
