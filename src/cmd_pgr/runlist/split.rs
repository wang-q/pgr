//! `pgr runlist split` — split a multi runlist JSON into per-key files.

use clap::{Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for split.
pub fn make_subcommand() -> Command {
    Command::new("split")
        .about("Splits a multi runlist JSON file")
        .after_help(
            r###"
Splits a multi runlist JSON into one JSON per top-level key, written to
`<outdir>/<key><suffix>` or printed line by line with `-o stdout`.

Examples:
1. Split into files:
   pgr runlist split in.json -o out_dir
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Sets the input file to use"),
        )
        .arg(
            Arg::new("suffix")
                .long("suffix")
                .short('s')
                .num_args(1)
                .default_value(".json")
                .help("Extensions of output files"),
        )
        .arg(
            Arg::new("outdir")
                .short('o')
                .long("outdir")
                .num_args(1)
                .default_value("stdout")
                .help("Output location. [stdout] for screen"),
        )
}

/// Execute the split command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let json = pgr::libs::runlist::read_json(args.get_one::<String>("infile").unwrap())?;
    let outdir = args.get_one::<String>("outdir").unwrap();
    if outdir != "stdout" {
        std::fs::create_dir_all(outdir)?;
    }
    let suffix = args.get_one::<String>("suffix").unwrap();
    let parts = pgr::libs::runlist::split_json(&json)?;
    if outdir == "stdout" {
        let mut w = pgr::writer("stdout")?;
        for (_, s) in parts {
            writeln!(w, "{}", s)?;
        }
        w.flush()?;
    } else {
        for (key, s) in parts {
            std::fs::write(std::path::Path::new(outdir).join(key + suffix), s + "\n")?;
        }
    }
    Ok(())
}
