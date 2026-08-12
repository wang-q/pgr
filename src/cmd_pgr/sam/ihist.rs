use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for ihist.
pub fn make_subcommand() -> Command {
    Command::new("ihist")
        .about("Computes an insert-size histogram from a paired SAM")
        .after_help(
            r###"
Reads a paired SAM (e.g. the `--paired` output of `anchr asm map`) and writes
the insert-size histogram in the BBTools `reformat.sh ihist` text format:
`#Mean`/`#Median`/`#Mode`/`#STDev`/`#PercentOfPairs` lines followed by
`#InsertSize	Count` rows.

Notes:
* Pairs are grouped by read name (first whitespace token, trailing
  `/1`/`/2` stripped)
* Only proper FR pairs (both ends mapped, same reference, opposite
  strands, pointing inward) contribute an insert size
* `#PercentOfPairs` is the fraction of pairs contributing to the histogram
* The median is the lower median; the mode is the most frequent size
  (ties -> the smallest)

Examples:
1. Insert-size histogram from a paired mapping (anchr 2_insert_size step):
   pgr sam ihist mapped.sam -o insert_size.ihist.txt

2. From stdin:
   pgr sam ihist stdin -o insert_size.ihist.txt
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input SAM file. [stdin] for standard input",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the ihist command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let output = crate::cmd_pgr::args::get_outfile(args);

    let reader =
        pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    pgr::libs::fmt::sam::ihist(reader, &mut writer)?;

    writer.flush()?;
    Ok(())
}
