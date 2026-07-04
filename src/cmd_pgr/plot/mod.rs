use clap::{ArgMatches, Command};

pub mod hh;
pub mod nrps;
pub mod venn;
/// Build the clap subcommand for plot.
pub fn make_subcommand() -> Command {
    Command::new("plot")
        .about("Plots figures")
        .subcommand_required(true)
        .subcommand(hh::make_subcommand())
        .subcommand(nrps::make_subcommand())
        .subcommand(venn::make_subcommand())
}
/// Execute the plot command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("hh", sub_matches)) => hh::execute(sub_matches),
        Some(("nrps", sub_matches)) => nrps::execute(sub_matches),
        Some(("venn", sub_matches)) => venn::execute(sub_matches),
        _ => Ok(()),
    }
}
