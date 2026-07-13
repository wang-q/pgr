use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use std::collections::BTreeMap;
/// Build the clap subcommand for cover.
pub fn make_subcommand() -> Command {
    Command::new("cover")
        .about("Outputs covered regions on chromosomes")
        .after_help(
            r###"
Outputs covered regions on chromosomes from block FA files.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* The output is in JSON format, showing the coverage of sequences on chromosomes
* Optionally, you can specify a species name to limit the output to that species
* For lastz results, use --trim 10

Examples:
1. Calculate coverage for all species:
   pgr fas cover tests/fas/example.fas

2. Calculate coverage for a specific species:
   pgr fas cover tests/fas/example.fas --name S288c

3. Trim alignment borders to avoid overlaps:
   pgr fas cover tests/fas/example.fas --trim 10

4. Output results to a file:
   pgr fas cover tests/fas/example.fas -o output.json

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::fas_name_arg(
            "Only output regions for this species",
        ))
        .arg(
            Arg::new("trim")
                .long("trim")
                .num_args(1)
                .value_parser(value_parser!(i32))
                .default_value("0")
                .help("Trim align borders to avoid overlaps"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the cover command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opt_trim = *args.get_one::<i32>("trim").unwrap();
    let opt_name = args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("");

    let mut res_of: BTreeMap<String, BTreeMap<String, intspan::IntSpan>> = BTreeMap::new();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
        pgr::libs::fmt::fas::aggregate_coverage_into(&mut reader, &mut res_of, opt_name, opt_trim)?;
    }

    let out_json = if !opt_name.is_empty() {
        // Output coverage for a single species
        intspan::set2json(
            res_of
                .values()
                .next()
                .ok_or_else(|| anyhow::anyhow!("no coverage found for --name {}", opt_name))?,
        )
    } else {
        // Output coverage for all species
        intspan::set2json_m(&res_of)
    };
    // Write the JSON output to the specified file or stdout
    intspan::write_json(crate::cmd_pgr::args::get_outfile(args), &out_json)?;

    Ok(())
}
