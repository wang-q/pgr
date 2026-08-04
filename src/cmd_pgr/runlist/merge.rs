//! `pgr runlist merge` — merge runlist JSON files into a multi runlist.

use clap::{Arg, ArgAction, ArgMatches, Command};

/// Build the clap subcommand for merge.
pub fn make_subcommand() -> Command {
    Command::new("merge")
        .about("Merges runlist JSON files")
        .after_help(
            r###"
Reads several runlist JSON files and writes a multi runlist keyed by file
stem. Without `--all` only the first dot-separated segment of the stem is
used as the key.

Examples:
1. Merge with short keys:
   pgr runlist merge a.json b.json -o out.json
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("all")
                .long("all")
                .action(ArgAction::SetTrue)
                .help("Use the full file stem as the key (without --all only the first dot-separated part is used)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the merge command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let files: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    // The output is a multi runlist JSON, not a single runlist; refuse to
    // overwrite an input file.
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, files.iter().map(String::as_str))?;
    let out = pgr::libs::runlist::merge_files(&files, args.get_flag("all"))?;
    pgr::libs::ds::intspan::write_json(outfile, &out)?;
    Ok(())
}
