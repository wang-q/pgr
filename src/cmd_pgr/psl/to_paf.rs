use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};

/// Build the clap subcommand for to-paf.
pub fn make_subcommand() -> Command {
    Command::new("to-paf")
        .about("Converts PSL alignments to PAF")
        .after_help(
            r###"
Convert UCSC PSL alignments to PAF format (12 mandatory columns).

* The PAF strand is the first character of the PSL strand.
* The PAF block length is the sum of the PSL block sizes.
* Mapping quality is set to 255.

Examples:
    pgr psl to-paf input.psl -o output.paf

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PSL file. [stdin] for standard input",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(ArgAction::SetTrue)
                .help("Fail on malformed records instead of skipping"),
        )
}

/// Execute the to-paf command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let strict = args.get_flag("strict");

    let reader = pgr::reader(infile)?;
    let mut writer = pgr::writer(outfile)?;
    pgr::libs::fmt::psl::to_paf(reader, &mut writer, strict)
        .with_context(|| format!("failed to convert {infile} to PAF"))?;
    Ok(())
}
