pub mod tofa;
pub mod interleave;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("fq")
        .about("Fastq tools")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(tofa::make_subcommand())
        .subcommand(interleave::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("tofa", sub_matches)) => tofa::execute(sub_matches),
        Some(("interleave", sub_matches)) | Some(("il", sub_matches)) => {
            interleave::execute(sub_matches)
        }
        _ => Ok(()),
    }
}
