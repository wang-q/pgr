use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
/// Build the clap subcommand for to-psl.
pub fn make_subcommand() -> Command {
    Command::new("to-psl")
        .about("Converts from lav to psl format")
        .after_help(
            r###"
Convert blastz lav to psl format.

Examples:
1. Convert lav to psl:
   pgr lav to-psl in.lav -o out.psl
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input LAV file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("target_strand")
                .long("target-strand")
                .help("Set the target strand (default is no strand)"),
        )
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(clap::ArgAction::SetTrue)
                .help("Fail on unknown LAV stanzas instead of warning and skipping"),
        )
}
/// Execute the to-psl command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let target_strand = args.get_one::<String>("target_strand");
    let strict = args.get_flag("strict");

    let reader =
        pgr::reader(input).with_context(|| format!("Failed to open reader for {}", input))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    pgr::libs::fmt::lav::lav_to_psl(
        reader,
        &mut writer,
        target_strand.map(|s| s.as_str()),
        strict,
    )?;

    Ok(())
}
