pub mod ir;
pub mod p2m;
pub mod prefilter;
pub mod rept;
pub mod trf;
pub mod ucsc;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("pl")
        .about("Run integrated pipelines")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(p2m::make_subcommand())
        .subcommand(prefilter::make_subcommand())
        .subcommand(trf::make_subcommand())
        .subcommand(ir::make_subcommand())
        .subcommand(rept::make_subcommand())
        .subcommand(ucsc::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("p2m", sub_matches)) => p2m::execute(sub_matches),
        Some(("prefilter", sub_matches)) => prefilter::execute(sub_matches),
        Some(("trf", sub_matches)) => trf::execute(sub_matches),
        Some(("ir", sub_matches)) => ir::execute(sub_matches),
        Some(("rept", sub_matches)) => rept::execute(sub_matches),
        Some(("ucsc", sub_matches)) => ucsc::execute(sub_matches),
        _ => unreachable!(),
    }
}
