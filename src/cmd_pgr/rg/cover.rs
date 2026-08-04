//! `pgr rg cover` — merge `.rg` lines into a runlist JSON.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for cover.
pub fn make_subcommand() -> Command {
    Command::new("cover")
        .about("Merges .rg range lines into runlist JSON")
        .after_help(
            r###"
Reads `chr:start-end` lines (1-based inclusive; species/strand prefixes are
dropped) from one or more `.rg` files and writes the per-chromosome union as
a runlist JSON ready for `pgr fa mask`.

Examples:
1. Merge ranges:
   pgr rg cover a.rg b.rg -o out.json
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

/// Execute the cover command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let files: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    // The output is runlist JSON, not `.rg`; refuse to overwrite an input.
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, files.iter().map(String::as_str))?;
    let set = pgr::libs::runlist::rg_files_to_set(&files)?;
    let json = pgr::libs::ds::intspan::set2json(&set);
    pgr::libs::ds::intspan::write_json(outfile, &json)?;
    Ok(())
}
