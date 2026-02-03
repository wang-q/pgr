pub mod anti_repeat;
pub mod net;
pub mod pre_net;
pub mod sort;
pub mod split;
pub mod stitch;

use clap::{ArgMatches, Command};

pub fn make_subcommand() -> Command {
    Command::new("chain")
        .about("Chain tools")
        .subcommand_required(true)
        .subcommand(anti_repeat::make_subcommand())
        .subcommand(split::make_subcommand())
        .subcommand(stitch::make_subcommand())
        .subcommand(sort::make_subcommand())
        .subcommand(pre_net::make_subcommand())
        .subcommand(net::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("anti-repeat", sub_matches)) => anti_repeat::execute(sub_matches),
        Some(("split", sub_matches)) => split::execute(sub_matches),
        Some(("stitch", sub_matches)) => stitch::execute(sub_matches),
        Some(("sort", sub_matches)) => sort::execute(sub_matches),
        Some(("pre-net", sub_matches)) => pre_net::execute(sub_matches),
        Some(("net", sub_matches)) => net::execute(sub_matches),
        _ => unreachable!(),
    }
}
