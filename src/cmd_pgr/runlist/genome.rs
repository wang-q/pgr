//! `pgr runlist genome` — convert a chromosome sizes file to a runlist JSON.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for genome.
pub fn make_subcommand() -> Command {
    Command::new("genome")
        .about("Converts a chr.sizes file to runlists")
        .after_help(
            r###"
Builds a runlist JSON where every chromosome spans its full length (1..size).

Examples:
1. Full-genome runlist:
   pgr runlist genome chr.sizes -o genome.json
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Sets the input file to use"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the genome command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let sizes = pgr::read_sizes::<i32>(args.get_one::<String>("infile").unwrap())?;
    let set = pgr::libs::runlist::genome_set(&sizes)?;
    let json = pgr::libs::ds::intspan::set2json(&set);
    pgr::libs::ds::intspan::write_json(outfile, &json)?;
    Ok(())
}
