//! `pgr sd search` — LASTZ-based putative SD detection.

use clap::{value_parser, Arg, ArgMatches, Command};

/// Build the clap subcommand for search.
pub fn make_subcommand() -> Command {
    Command::new("search")
        .about("Detects putative segmental duplications via lastz self-alignment")
        .after_help(
            r###"
Runs `lastz --self` on a genome, converts the LAV output to PSL, and keeps
hits meeting the T2T-CHM13 SD standard (> 1 kbp, > 90% identity; see
notes/references/biser.md §4.2.1). The output PSL is NOT chained - feed it
through `pgr pl ucsc` (without --syn) for chain/net refinement, then
`pgr maf to-paf` for the downstream PAF.

Notes:
* Requires lastz in PATH.
* Input FASTA may be plain or gzipped (.gz).
* Coordinates in the output PSL are 0-based half-open on the input genome.
* Identity is computed as (matches + rep_matches) / block_length, where
  block_length includes insert bases (unlike gap-compressed identity).

Examples:
1. Detect SDs in a bacterial genome:
   pgr sd search genome.fa -o hits.psl

2. Relax/tighten filters:
   pgr sd search genome.fa -o hits.psl --min-len 500 --min-identity 0.85
"###,
        )
        .arg(
            Arg::new("genome")
                .index(1)
                .required(true)
                .help("Genome FASTA file (plain or .gz)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("preset")
                .long("preset")
                .value_parser(clap::builder::PossibleValuesParser::new(
                    pgr::libs::lastz::preset_names(),
                ))
                .help("lastz parameter set (default: set01)"),
        )
        .arg(
            Arg::new("min_len")
                .long("min-len")
                .default_value("1000")
                .value_parser(value_parser!(u32))
                .help("Minimum alignment block length in bp"),
        )
        .arg(
            Arg::new("min_identity")
                .long("min-identity")
                .default_value("0.90")
                .value_parser(value_parser!(f64))
                .help("Minimum block identity (0.0-1.0)"),
        )
        .arg(
            Arg::new("query_depth")
                .long("query-depth")
                .default_value("50")
                .value_parser(value_parser!(usize))
                .help("lastz query depth threshold"),
        )
        .arg(crate::cmd_pgr::args::parallel_arg_with_default("4"))
}
/// Execute the search command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let genome = args.get_one::<String>("genome").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let opts = pgr::libs::sd::search_lastz::SearchLastzOptions {
        preset: args.get_one::<String>("preset").cloned(),
        min_len: *args.get_one::<u32>("min_len").unwrap(),
        min_identity: *args.get_one::<f64>("min_identity").unwrap(),
        query_depth: *args.get_one::<usize>("query_depth").unwrap(),
        parallel: *args.get_one::<usize>("parallel").unwrap(),
    };

    let workdir = tempfile::tempdir()?;
    let hits =
        pgr::libs::sd::search_lastz::search_lastz(genome, workdir.path().to_str().unwrap(), &opts)?;

    let mut writer = pgr::writer(outfile)?;
    for psl in &hits {
        psl.write_to(&mut writer)?;
    }
    Ok(())
}
