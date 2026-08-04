//! `pgr rg merge` — cluster nearly-identical `.rg` ranges.

use clap::{Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for merge.
pub fn make_subcommand() -> Command {
    Command::new("merge")
        .about("Clusters nearly-identical .rg ranges and emits mappings")
        .after_help(
            r###"
Clusters `.rg` ranges whose reciprocal overlap reaches `--coverage` and emits
`range<TAB>merged` mapping lines for ranges in multi-member clusters. The
merged representative is the union cover `chr(+):min-max`. Ranges not joined
with any other are omitted. Migrated from the external `rgr merge` (adapted
to single-column `.rg` input; rgr's multi-part TSV use case is out of scope).

Examples:
1. Cluster with the default 0.95 reciprocal overlap:
   pgr rg merge a.rg
2. Looser threshold:
   pgr rg merge a.rg b.rg --coverage 0.90 -o map.tsv
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input .rg files to process"),
        )
        .arg(
            Arg::new("coverage")
                .long("coverage")
                .short('c')
                .num_args(1)
                .default_value("0.95")
                .value_parser(clap::value_parser!(f32))
                .help("Minimum reciprocal overlap to join two ranges"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the merge command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let coverage = *args.get_one::<f32>("coverage").unwrap();
    let files: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    // The output is a `part<TAB>merged` mapping, not `.rg`; refuse to
    // overwrite an input file.
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, files.iter().map(String::as_str))?;
    let mapping = pgr::libs::runlist::rg_merge_mapping(&files, coverage)?;

    let mut writer = pgr::writer(outfile)?;
    for (part, merged) in &mapping {
        writeln!(writer, "{}\t{}", part, merged)?;
    }
    writer.flush()?;
    Ok(())
}
