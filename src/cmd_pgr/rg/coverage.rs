//! `pgr rg coverage` — depth of coverage over `.rg` ranges.

use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use pgr::libs::ds::IntSpan;
use std::collections::BTreeMap;

/// Build the clap subcommand for coverage.
pub fn make_subcommand() -> Command {
    Command::new("coverage")
        .about("Computes depth of coverage over .rg ranges")
        .after_help(
            r###"
Reads `chr:start-end` lines (1-based inclusive) from one or more `.rg` files,
computes per-position coverage depth with a sweep line, and writes regions
whose depth reaches `--minimum`. With `--detailed`, regions are grouped by
their exact depth instead.

Examples:
1. Regions covered at least 4 times:
   pgr rg coverage a.rg -m 4 -o cov.json
2. Per-depth regions:
   pgr rg coverage a.rg -m 2 -d -o cov.json
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input .rg file(s) to process"),
        )
        .arg(
            Arg::new("minimum")
                .long("minimum")
                .short('m')
                .value_parser(value_parser!(u32))
                .num_args(1)
                .default_value("1")
                .help("Set the minimum depth of coverage"),
        )
        .arg(
            Arg::new("detailed")
                .long("detailed")
                .short('d')
                .action(ArgAction::SetTrue)
                .help("Output detailed depth"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the coverage command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let minimum = *args.get_one::<u32>("minimum").unwrap();
    let is_detailed = args.get_flag("detailed");

    let mut iv_of: BTreeMap<String, Vec<(u32, u32)>> = BTreeMap::new();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile)?;
        for (chr, ivs) in pgr::libs::runlist::rg_to_intervals(reader)? {
            iv_of.entry(chr).or_default().extend(ivs);
        }
    }

    if is_detailed {
        let mut set_of: BTreeMap<String, BTreeMap<String, IntSpan>> = BTreeMap::new();
        for (chr, ivs) in &iv_of {
            let by_level = pgr::libs::runlist::depth_by_level(ivs, minimum);
            for (depth, is) in by_level {
                set_of.entry(depth).or_default().insert(chr.clone(), is);
            }
        }
        let json = pgr::libs::ds::intspan::set2json_m(&set_of);
        pgr::libs::ds::intspan::write_json(outfile, &json)?;
    } else {
        let mut set: BTreeMap<String, IntSpan> = BTreeMap::new();
        for (chr, ivs) in &iv_of {
            set.insert(
                chr.clone(),
                pgr::libs::runlist::depth_at_least(ivs, minimum),
            );
        }
        let json = pgr::libs::ds::intspan::set2json(&set);
        pgr::libs::ds::intspan::write_json(outfile, &json)?;
    }
    Ok(())
}
