pub mod count;
pub mod filter;
pub mod gz;
pub mod masked;
pub mod n50;
pub mod one;
pub mod order;
pub mod range;
pub mod rc;
pub mod replace;
pub mod size;
pub mod some;
pub mod split;
pub mod to2bit;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("fa")
        .about("Fasta tools")
        .after_help(
            r###"Notes:
* Supports both plain text and gzipped (.gz) files

"###,
        )
        .subcommand(count::make_subcommand())
        .subcommand(filter::make_subcommand())
        .subcommand(gz::make_subcommand())
        .subcommand(masked::make_subcommand())
        .subcommand(n50::make_subcommand())
        .subcommand(one::make_subcommand())
        .subcommand(order::make_subcommand())
        .subcommand(range::make_subcommand())
        .subcommand(rc::make_subcommand())
        .subcommand(replace::make_subcommand())
        .subcommand(size::make_subcommand())
        .subcommand(some::make_subcommand())
        .subcommand(split::make_subcommand())
        .subcommand(to2bit::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("count", sub_matches)) => count::execute(sub_matches),
        Some(("filter", sub_matches)) => filter::execute(sub_matches),
        Some(("gz", sub_matches)) => gz::execute(sub_matches),
        Some(("masked", sub_matches)) => masked::execute(sub_matches),
        Some(("n50", sub_matches)) => n50::execute(sub_matches),
        Some(("one", sub_matches)) => one::execute(sub_matches),
        Some(("order", sub_matches)) => order::execute(sub_matches),
        Some(("range", sub_matches)) => range::execute(sub_matches),
        Some(("rc", sub_matches)) => rc::execute(sub_matches),
        Some(("replace", sub_matches)) => replace::execute(sub_matches),
        Some(("size", sub_matches)) => size::execute(sub_matches),
        Some(("some", sub_matches)) => some::execute(sub_matches),
        Some(("split", sub_matches)) => split::execute(sub_matches),
        Some(("to2bit", sub_matches)) => to2bit::execute(sub_matches),
        _ => Ok(()),
    }
}
