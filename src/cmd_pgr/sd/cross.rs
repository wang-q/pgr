//! `pgr sd cross` — cross-genome SD mapping (BISER cross_search/cross_align).

use clap::{Arg, ArgMatches, Command};
use pgr::libs::sd::search_lastz::{lastz_to_hits, SearchLastzOptions};
use pgr::libs::sd::search_pgi::{pgi_to_hits, SearchPgiOptions};
use std::io::Write;

/// Build the clap subcommand for cross.
pub fn make_subcommand() -> Command {
    Command::new("cross")
        .about("Maps SD-like homology from one genome to another")
        .after_help(
            r###"
Cross-genome counterpart of search+align (BISER cross_search / cross_align):
aligns the second genome as query against the first as target, filters to the
T2T-CHM13 SD standard, then refines via chain/net (no --syn) and merges the
MAF into one PAF.

Two engines (same as `pgr sd search`):
* `pgi` (default): native `pgr align pgi`; no external tools.
* `lastz`: external `lastz`; `--preset`/`--query-depth` apply to this engine.

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
                .help("lastz parameter set (default: set01; lastz engine only)"),
        )
        .arg(
            Arg::new("query_depth")
                .long("query-depth")
                .default_value("50")
                .value_parser(clap::value_parser!(usize))
                .help("lastz query depth threshold (lastz engine only)"),
        )
        .arg(
            Arg::new("engine")
                .long("engine")
                .value_parser(["pgi", "lastz"])
                .default_value("pgi")
                .help("Alignment engine: native pgi (default) or external lastz"),
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
        .arg(crate::cmd_pgr::args::parallel_arg_with_default("4"))
}
/// Execute the cross command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let target = args.get_one::<String>("target").unwrap();
    let query = args.get_one::<String>("query").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let engine = args.get_one::<String>("engine").unwrap().as_str();
    let min_len = *args.get_one::<u32>("min_len").unwrap();
    let min_identity = *args.get_one::<f64>("min_identity").unwrap();
    let parallel = *args.get_one::<usize>("parallel").unwrap();

    let workdir = tempfile::tempdir()?;
    let hits = match engine {
        "pgi" => {
            let opts = SearchPgiOptions {
                min_len,
                min_identity,
                parallel,
            };
            pgi_to_hits(
                target,
                query,
                false,
                workdir.path().to_str().unwrap(),
                &opts,
            )?
        }
        _ => {
            let opts = SearchLastzOptions {
                preset: args.get_one::<String>("preset").cloned(),
                query_depth: *args.get_one::<usize>("query_depth").unwrap(),
                min_len,
                min_identity,
                parallel,
            };
            lastz_to_hits(
                target,
                query,
                false,
                workdir.path().to_str().unwrap(),
                &opts,
            )?
        }
    };
    let psl_path = workdir.path().join("hits.psl");
    let mut w = std::io::BufWriter::new(std::fs::File::create(&psl_path)?);
    for psl in &hits {
        psl.write_to(&mut w)?;
    }
    w.flush()?;
    super::chainnet_to_paf(target, query, psl_path.to_str().unwrap(), outfile)
}
