use clap::Command;

pub mod compare;
pub mod format;
pub mod nj;
pub mod subset;
pub mod to_pair;
pub mod to_phylip;
pub mod upgma;

pub fn make_subcommand() -> Command {
    Command::new("mat")
        .about("Matrix operations")
        .subcommand_required(true)
        .subcommand(compare::make_subcommand())
        .subcommand(format::make_subcommand())
        .subcommand(to_pair::make_subcommand())
        .subcommand(to_phylip::make_subcommand())
        .subcommand(subset::make_subcommand())
        .subcommand(upgma::make_subcommand())
        .subcommand(nj::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("compare", sub_matches)) => compare::execute(sub_matches),
        Some(("format", sub_matches)) => format::execute(sub_matches),
        Some(("to-pair", sub_matches)) => to_pair::execute(sub_matches),
        Some(("to-phylip", sub_matches)) => to_phylip::execute(sub_matches),
        Some(("subset", sub_matches)) => subset::execute(sub_matches),
        Some(("upgma", sub_matches)) => upgma::execute(sub_matches),
        Some(("nj", sub_matches)) => nj::execute(sub_matches),
        _ => unreachable!(),
    }
}
