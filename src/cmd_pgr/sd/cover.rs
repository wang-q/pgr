//! `pgr sd cover` — mark core duplicons via greedy set cover.

use anyhow::Context;
use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for cover.
pub fn make_subcommand() -> Command {
    Command::new("cover")
        .about("Marks core duplicons via greedy set cover over elementary SDs")
        .after_help(
            r###"
Reads the refined SD hits (`hits.paf` from `pgr sd align`) and the elementary
SD BED (`elems.bed`, merged `pgr sd decompose` output with genome
coordinates). An elementary SD set (same `set_id`) covers a hit if any copy
overlaps the hit's query or target interval. A greedy set cover selects the
smallest elementary sets that cover all hits; selected rows are marked CORE.

Examples:
1. Mark core duplicons:
   pgr sd cover hits.paf elems.bed -o out.elem.bed
"###,
        )
        .arg(
            Arg::new("paf")
                .index(1)
                .required(true)
                .help("Refined SD hits PAF file"),
        )
        .arg(
            Arg::new("elems")
                .index(2)
                .required(true)
                .help("Elementary SD BED file (merged decompose output)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the cover command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let paf = args.get_one::<String>("paf").unwrap();
    let elems = args.get_one::<String>("elems").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let hits_reader = pgr::reader(paf).with_context(|| format!("failed to open PAF {}", paf))?;
    let elems_reader =
        pgr::reader(elems).with_context(|| format!("failed to open BED {}", elems))?;
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("failed to open writer for {}", outfile))?;
    pgr::libs::sd::cover::run_cover(hits_reader, elems_reader, &mut writer)?;
    Ok(())
}
