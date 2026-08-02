//! `pgr sd search` — putative SD detection via pgi or lastz self-alignment.

use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use cmd_lib::run_cmd;
use std::io::BufRead;

/// Build the clap subcommand for search.
pub fn make_subcommand() -> Command {
    Command::new("search")
        .about("Detects putative segmental duplications via self-alignment")
        .after_help(
            r###"
Detects putative segmental duplications by self-aligning a genome and keeping
hits meeting the T2T-CHM13 SD standard (> 1 kbp, > 90% identity; see
notes/references/biser.md §4.2.1). The output PSL is NOT chained - feed it
through `pgr sd align` for chain/net refinement and the downstream PAF.

Two engines:
* `pgi` (default): native `pgr align pgi` self-alignment (FastGA-style
  syncmer seeds + tube chaining + banded/wave extension). No external tools.
* `lastz`: external `lastz --self` (LAV -> PSL); requires lastz in PATH.
  `--preset`/`--query-depth` apply to this engine only.

Notes:
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
                .help("lastz parameter set (default: set01; lastz engine only)"),
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
                .help("lastz query depth threshold (lastz engine only)"),
        )
        .arg(
            Arg::new("engine")
                .long("engine")
                .value_parser(["pgi", "lastz"])
                .default_value("pgi")
                .help("Self-alignment engine: native pgi (default) or external lastz"),
        )
        .arg(crate::cmd_pgr::args::parallel_arg_with_default("4"))
}
/// Execute the search command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let genome = args.get_one::<String>("genome").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let engine = args.get_one::<String>("engine").unwrap().as_str();
    let min_len = *args.get_one::<u32>("min_len").unwrap();
    let min_identity = *args.get_one::<f64>("min_identity").unwrap();
    let parallel = *args.get_one::<usize>("parallel").unwrap();

    let workdir = tempfile::tempdir()?;
    let hits = match engine {
        "pgi" => search_pgi(genome, workdir.path(), parallel, min_len, min_identity)?,
        _ => {
            let opts = pgr::libs::sd::search_lastz::SearchLastzOptions {
                preset: args.get_one::<String>("preset").cloned(),
                query_depth: *args.get_one::<usize>("query_depth").unwrap(),
                min_len,
                min_identity,
                parallel,
            };
            pgr::libs::sd::search_lastz::search_lastz(
                genome,
                workdir.path().to_str().unwrap(),
                &opts,
            )?
        }
    };

    let mut writer = pgr::writer(outfile)?;
    for psl in &hits {
        psl.write_to(&mut writer)?;
    }
    Ok(())
}

/// Run the native pgi self-alignment and keep hits passing the SD filters.
fn search_pgi(
    genome: &str,
    workdir: &std::path::Path,
    parallel: usize,
    min_len: u32,
    min_identity: f64,
) -> anyhow::Result<Vec<pgr::libs::fmt::psl::Psl>> {
    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_sd_search_pgi_")?;
    let pgr = ctx.pgr.clone();
    let abs_genome = ctx.abs_path(genome)?;
    let raw = ctx.abs_path(&workdir.join("hits.raw.psl").to_string_lossy())?;
    let _cwd_guard = ctx.enter()?;
    run_cmd!(${pgr} align pgi ${abs_genome} -o ${raw} --parallel ${parallel})?;

    let mut hits = Vec::new();
    let mut reader =
        pgr::reader(&raw).with_context(|| format!("failed to open pgi SD hits {}", raw))?;
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if let Some(psl) = pgr::libs::fmt::psl::parse_or_warn(line.trim_end(), false)? {
            if pgr::libs::sd::search_lastz::passes_sd_filters(&psl, min_len, min_identity) {
                hits.push(psl);
            }
        }
    }
    Ok(hits)
}
