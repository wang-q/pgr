pub mod chain;
pub mod histo;
pub mod lift;
pub mod rc;
pub mod stats;
pub mod swap;
pub mod to_chain;
pub mod to_range;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("psl")
        .about("Manipulate PSL alignment files")
        .after_help(
            r###"Note:
These utilities are primarily provided for cross-validation with the original UCSC Kent tools,
ensuring the fidelity of the ported libraries.
"###,
        )
        .subcommand(chain::make_subcommand())
        .subcommand(histo::make_subcommand())
        .subcommand(lift::make_subcommand())
        .subcommand(rc::make_subcommand())
        .subcommand(stats::make_subcommand())
        .subcommand(swap::make_subcommand())
        .subcommand(to_chain::make_subcommand())
        .subcommand(to_range::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("chain", sub_matches)) => chain::execute(sub_matches),
        Some(("histo", sub_matches)) => histo::execute(sub_matches),
        Some(("lift", sub_matches)) => lift::execute(sub_matches),
        Some(("rc", sub_matches)) => rc::execute(sub_matches),
        Some(("stats", sub_matches)) => stats::execute(sub_matches),
        Some(("swap", sub_matches)) => swap::execute(sub_matches),
        Some(("to-chain", sub_matches)) => to_chain::execute(sub_matches),
        Some(("to-range", sub_matches)) => to_range::execute(sub_matches),
        _ => Ok(()),
    }
}
