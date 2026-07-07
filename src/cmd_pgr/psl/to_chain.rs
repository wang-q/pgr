use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use std::io::Write;
/// Build the clap subcommand for to-chain.
pub fn make_subcommand() -> Command {
    Command::new("to-chain")
        .about("Converts PSL to Chain format")
        .after_help(
            r###"
Examples:
1. Convert PSL to Chain:
   pgr psl to-chain in.psl -o out.chain
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input PSL file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("fix_strand")
                .long("fix-strand")
                .action(clap::ArgAction::SetTrue)
                .help("Fix '-' target strand by reverse complementing the record"),
        )
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(clap::ArgAction::SetTrue)
                .help("Fail on parse errors instead of skipping malformed lines"),
        )
}
/// Execute the to-chain command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let fix_strand = args.get_flag("fix_strand");
    let strict = args.get_flag("strict");

    let reader =
        pgr::reader(input).with_context(|| format!("Failed to open reader for {}", input))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    pgr::libs::fmt::psl::to_chain(reader, &mut writer, fix_strand, strict)?;

    writer.flush()?;
    Ok(())
}
