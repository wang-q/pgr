pub mod hv;
pub mod seq;
pub mod vector;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("dist")
        .about("Distance/Similarity metrics")
        .after_help(
            r###"Subcommand groups:

* distance: hv / seq / vector

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(hv::make_subcommand())
        .subcommand(seq::make_subcommand())
        .subcommand(vector::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("hv", sub_matches)) => hv::execute(sub_matches),
        Some(("seq", sub_matches)) => seq::execute(sub_matches),
        Some(("vector", sub_matches)) => vector::execute(sub_matches),
        _ => Ok(()),
    }
}
