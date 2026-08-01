//! `pgr sd cross` — cross-genome SD mapping (BISER cross_search/cross_align).

use clap::{Arg, ArgMatches, Command};
use pgr::libs::sd::search_lastz::{lastz_to_hits, SearchLastzOptions};
use std::io::Write;

/// Build the clap subcommand for cross.
pub fn make_subcommand() -> Command {
    Command::new("cross")
        .about("Maps SD-like homology from one genome to another")
        .after_help(
            r###"
Cross-genome counterpart of search+align (BISER cross_search / cross_align on
the external-alignment route): runs lastz with the second genome as query
against the first as target, filters to the T2T-CHM13 SD standard, then
refines via chain/net (no --syn) and merges the MAF into one PAF.

Examples:
1. Map genome B's homology onto genome A:
   pgr sd cross A.fa B.fa -o cross.paf
"###,
        )
        .arg(
            Arg::new("target")
                .index(1)
                .required(true)
                .help("Target genome FASTA"),
        )
        .arg(
            Arg::new("query")
                .index(2)
                .required(true)
                .help("Query genome FASTA"),
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
                .value_parser(clap::value_parser!(u32))
                .help("Minimum alignment block length in bp"),
        )
        .arg(
            Arg::new("min_identity")
                .long("min-identity")
                .default_value("0.90")
                .value_parser(clap::value_parser!(f64))
                .help("Minimum block identity (0.0-1.0)"),
        )
}
/// Execute the cross command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let target = args.get_one::<String>("target").unwrap();
    let query = args.get_one::<String>("query").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let opts = SearchLastzOptions {
        preset: args.get_one::<String>("preset").cloned(),
        min_len: *args.get_one::<u32>("min_len").unwrap(),
        min_identity: *args.get_one::<f64>("min_identity").unwrap(),
        ..Default::default()
    };

    let workdir = tempfile::tempdir()?;
    let hits = lastz_to_hits(
        target,
        query,
        false,
        workdir.path().to_str().unwrap(),
        &opts,
    )?;
    let psl_path = workdir.path().join("hits.psl");
    let mut w = std::io::BufWriter::new(std::fs::File::create(&psl_path)?);
    for psl in &hits {
        psl.write_to(&mut w)?;
    }
    w.flush()?;
    super::chainnet_to_paf(target, query, psl_path.to_str().unwrap(), outfile)
}
