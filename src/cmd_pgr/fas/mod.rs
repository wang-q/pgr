pub mod concat;
pub mod consensus;
pub mod cover;
pub mod join;
pub mod link;
pub mod name;
pub mod refine;
pub mod replace;
pub mod separate;
pub mod slice;
pub mod split;
pub mod stat;
pub mod subset;
pub mod variation;
pub mod toxlsx;
pub mod tovcf;
pub mod filter;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("fas")
        .about("Block FA tools")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(concat::make_subcommand())
        .subcommand(consensus::make_subcommand())
        .subcommand(cover::make_subcommand())
        .subcommand(filter::make_subcommand())
        .subcommand(join::make_subcommand())
        .subcommand(link::make_subcommand())
        .subcommand(name::make_subcommand())
        .subcommand(refine::make_subcommand())
        .subcommand(replace::make_subcommand())
        .subcommand(separate::make_subcommand())
        .subcommand(slice::make_subcommand())
        .subcommand(split::make_subcommand())
        .subcommand(stat::make_subcommand())
        .subcommand(variation::make_subcommand())
        .subcommand(toxlsx::make_subcommand())
        .subcommand(tovcf::make_subcommand())
        .subcommand(subset::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("concat", sub_matches)) => concat::execute(sub_matches),
        Some(("consensus", sub_matches)) => consensus::execute(sub_matches),
        Some(("cover", sub_matches)) => cover::execute(sub_matches),
        Some(("filter", sub_matches)) => filter::execute(sub_matches),
        Some(("join", sub_matches)) => join::execute(sub_matches),
        Some(("link", sub_matches)) => link::execute(sub_matches),
        Some(("name", sub_matches)) => name::execute(sub_matches),
        Some(("refine", sub_matches)) => refine::execute(sub_matches),
        Some(("replace", sub_matches)) => replace::execute(sub_matches),
        Some(("separate", sub_matches)) => separate::execute(sub_matches),
        Some(("slice", sub_matches)) => slice::execute(sub_matches),
        Some(("split", sub_matches)) => split::execute(sub_matches),
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        Some(("variation", sub_matches)) => variation::execute(sub_matches),
        Some(("toxlsx", sub_matches)) => toxlsx::execute(sub_matches),
        Some(("tovcf", sub_matches)) => tovcf::execute(sub_matches),
        Some(("subset", sub_matches)) => subset::execute(sub_matches),
        _ => unreachable!(),
    }
}
