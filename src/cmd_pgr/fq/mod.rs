pub mod interleave;
pub mod to_fa;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for fq.
pub fn make_subcommand() -> Command {
    Command::new("fq")
        .about("Manipulates FASTQ files")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(to_fa::make_subcommand())
        .subcommand(interleave::make_subcommand())
}
/// Execute the fq command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("to-fa", sub_matches)) => to_fa::execute(sub_matches),
        Some(("interleave", sub_matches)) | Some(("il", sub_matches)) => {
            interleave::execute(sub_matches)
        }
        _ => Ok(()),
    }
}
