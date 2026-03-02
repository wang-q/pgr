use clap::Command;

pub mod hh;
pub mod venn;

pub fn make_subcommand() -> Command {
    Command::new("plot")
        .about("Plotting tools")
        .subcommand_required(true)
        .subcommand(hh::make_subcommand())
        .subcommand(venn::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("hh", sub_matches)) => hh::execute(sub_matches),
        Some(("venn", sub_matches)) => venn::execute(sub_matches),
        _ => unreachable!(),
    }
}
