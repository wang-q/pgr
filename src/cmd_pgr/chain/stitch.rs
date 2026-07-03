use anyhow::Result;
use clap::Command;

pub fn make_subcommand() -> Command {
    Command::new("stitch")
        .about("Joins chain fragments with the same chain ID into a single chain per ID")
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input chain file",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg_required())
}

pub fn execute(args: &clap::ArgMatches) -> Result<()> {
    let input_path = args.get_one::<String>("infile").unwrap();
    let output_path = crate::cmd_pgr::args::get_outfile(args);
    let reader = pgr::reader(input_path)?;
    let writer = pgr::writer(output_path)?;
    pgr::libs::chain::stitch_chains(reader, writer)
}
