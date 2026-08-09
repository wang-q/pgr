pub mod interleave;
pub mod range;
pub mod to_fa;
pub mod trim_q;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for fq.
pub fn make_subcommand() -> Command {
    Command::new("fq")
        .about("Manipulates FASTQ files")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(to_fa::make_subcommand())
        .subcommand(interleave::make_subcommand())
        .subcommand(range::make_subcommand())
        .subcommand(trim_q::make_subcommand())
}
/// Execute the fq command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("to-fa", sub_matches)) => to_fa::execute(sub_matches),
        Some(("interleave", sub_matches)) | Some(("il", sub_matches)) => {
            interleave::execute(sub_matches)
        }
        Some(("range", sub_matches)) => range::execute(sub_matches),
        Some(("trim-q", sub_matches)) => trim_q::execute(sub_matches),
        _ => Ok(()),
    }
}
