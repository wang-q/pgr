use clap::*;

pub mod sort;
pub mod tomaf;
pub mod topsl;

pub fn make_subcommand() -> Command {
    Command::new("axt")
        .about("Axt tools")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(sort::make_subcommand())
        .subcommand(tomaf::make_subcommand())
        .subcommand(topsl::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("sort", sub_matches)) => sort::execute(sub_matches),
        Some(("tomaf", sub_matches)) => tomaf::execute(sub_matches),
        Some(("topsl", sub_matches)) => topsl::execute(sub_matches),
        _ => unreachable!(),
    }
}
