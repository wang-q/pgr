use clap::{ArgMatches, Command};
/// Build the clap subcommand for stitch.
pub fn make_subcommand() -> Command {
    Command::new("stitch")
        .about("Joins chain fragments with the same chain ID into a single chain per ID")
        .after_help(
            r###"
Joins chain fragments sharing the same chain ID into a single chain per ID,
mirroring the UCSC chainStitchId workflow. Chains with the same ID are
concatenated end-to-end in the order they appear in the input.

Examples:
1. Stitch chain fragments by ID:
   pgr chain stitch in.chain -o stitched.chain

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input chain file",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg_required())
}
/// Execute the stitch command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("infile").unwrap();
    let output_path = crate::cmd_pgr::args::get_outfile(args);
    let reader = pgr::reader(input_path)?;
    let writer = pgr::writer(output_path)?;
    pgr::libs::chain::stitch_chains(reader, writer)
}
