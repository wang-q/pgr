use clap::*;

pub mod sort;
pub mod to_fas;
pub mod to_maf;
pub mod to_psl;

pub fn make_subcommand() -> Command {
    Command::new("axt")
        .about("Manipulate AXT alignment files")
        .after_help(
            r###"
AXT is a format for representing pairwise genomic alignments.

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(sort::make_subcommand())
        .subcommand(to_maf::make_subcommand())
        .subcommand(to_fas::make_subcommand())
        .subcommand(to_psl::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("sort", sub_matches)) => sort::execute(sub_matches),
        Some(("to-maf", sub_matches)) => to_maf::execute(sub_matches),
        Some(("to-fas", sub_matches)) => to_fas::execute(sub_matches),
        Some(("to-psl", sub_matches)) => to_psl::execute(sub_matches),
        _ => unreachable!(),
    }
}
