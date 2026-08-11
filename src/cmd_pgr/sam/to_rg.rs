use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for to-rg.
pub fn make_subcommand() -> Command {
    Command::new("to-rg")
        .about("Extracts coordinates from SAM as .rg range lines")
        .after_help(
            r###"
Extract alignment coordinates from SAM files and output in .rg format
(chr:start-end, 1-based inclusive). This is useful for depth calculation
with `pgr rg coverage`, e.g. deriving per-base coverage from the mapped SAM
of `pgr asm map`.

Notes:
* Header lines and unmapped records are skipped.
* The range spans the full CIGAR (M/D/N/=/X operations consume reference).
* Malformed lines are skipped unless `--strict`.

Examples:
1. Convert mapped reads to ranges:
   pgr sam to-rg mapped.sam > mapped.rg

2. Derive per-base coverage from an asm map SAM (anchr anchors step):
   pgr asm map UT.fasta R1.fq.gz R2.fq.gz --outm mapped.sam --outu unmapped.sam
   pgr sam to-rg mapped.sam | pgr rg coverage stdin -m 2 -o cov.json
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input SAM file. [stdin] for standard input",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(ArgAction::SetTrue)
                .help("Fail on parse errors instead of skipping malformed lines"),
        )
}

/// Execute the to-rg command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let output = crate::cmd_pgr::args::get_outfile(args);
    let strict = args.get_flag("strict");

    let reader =
        pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    pgr::libs::fmt::sam::to_ranges(reader, &mut writer, strict)?;

    writer.flush()?;
    Ok(())
}
