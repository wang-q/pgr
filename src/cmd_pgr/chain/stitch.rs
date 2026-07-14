use anyhow::Context;
use clap::{ArgMatches, Command};
/// Build the clap subcommand for stitch.
pub fn make_subcommand() -> Command {
    Command::new("stitch")
        .about("Joins chain fragments with the same chain ID into a single chain per ID")
        .after_help(
            r###"
Joins chain fragments sharing the same chain ID into a single chain per ID,
mirroring the UCSC chainStitchId workflow.

Processing:
  1. Group input chains by their ID.
  2. For each group, check that target name, query name, and query strand are
     consistent; inconsistent fragments are skipped with a warning.
  3. Convert all fragments to blocks, sort the blocks by (target start, query start),
     and rebuild a single chain.
  4. Sum the scores of all fragments and assign the result to the stitched chain.

Notes:
* Fragments are concatenated in input order before sorting by (target start, query start).
* No overlap or abutment validation is performed between fragments.

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
    let reader = pgr::reader(input_path)
        .with_context(|| format!("Failed to open reader for {}", input_path))?;
    let writer = pgr::writer(output_path)
        .with_context(|| format!("Failed to open writer for {}", output_path))?;
    pgr::libs::chain::stitch_chains(reader, writer)
}
