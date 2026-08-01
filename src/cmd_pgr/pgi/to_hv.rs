//! `pgr pgi to-hv` — project an index onto a hypervector.

use clap::{value_parser, Arg, ArgMatches, Command};

/// Build the clap subcommand for to-hv.
pub fn make_subcommand() -> Command {
    Command::new("to-hv")
        .about("Projects a .pgi index onto a hypervector (.hv)")
        .after_help(
            r###"
Projects the index's unique k-mer keys onto a fixed-dimension hypervector,
enabling O(dim) distance comparisons for very large cohorts (see
notes/design/pbit-index-extension.md §6).

Examples:
1. Project to a 1024-dim hypervector:
   pgr pgi to-hv genome.pgi -o genome.hv --dim 1024
"###,
        )
        .arg(
            Arg::new("infile")
                .index(1)
                .required(true)
                .help(".pgi index file"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("dim")
                .long("dim")
                .default_value("1024")
                .value_parser(value_parser!(usize))
                .help("Hypervector dimension (multiple of 32)"),
        )
}
/// Execute the to-hv command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();
    let dim = *args.get_one::<usize>("dim").unwrap();
    anyhow::ensure!(
        dim > 0 && dim.is_multiple_of(32),
        "--dim must be a positive multiple of 32"
    );

    let mut reader = pgr::reader(infile)?;
    let idx = pgr::libs::pgi::PgiIndex::read(&mut reader)?;
    let hv = pgr::libs::pgi::to_hv::index_to_hv(&idx, dim);
    let name = pgr::libs::io::get_basename(infile).unwrap_or_else(|| infile.clone());
    let mut writer = pgr::writer(outfile)?;
    pgr::libs::pgi::to_hv::write_hv(&mut writer, &name, idx.k, dim, &hv)?;
    log::info!("wrote {}-dim hypervector to {}", dim, outfile);
    Ok(())
}
