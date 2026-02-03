pub mod interleave;
pub mod to_fa;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("fq")
        .about("Manipulate FASTQ files")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(to_fa::make_subcommand())
        .subcommand(interleave::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("to-fa", sub_matches)) => to_fa::execute(sub_matches),
        Some(("interleave", sub_matches)) | Some(("il", sub_matches)) => {
            interleave::execute(sub_matches)
        }
        _ => Ok(()),
    }
}
