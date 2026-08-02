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
notes/design/pbit.md). The projection is sparse: each
k-mer updates --sparse random dimensions, so the shared-k-mer signal stays
dominant for large k-mer sets and cosine similarity on the result
approximates the k-mer set overlap.

Examples:
1. Project to a 4096-dim hypervector (default):
   pgr pgi to-hv genome.pgi -o genome.hv
2. Higher dimension for closer approximation of `pgr dist pgi`:
   pgr pgi to-hv genome.pgi -o genome.hv --dim 16384
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
                .default_value("4096")
                .value_parser(value_parser!(usize))
                .help("Hypervector dimension (multiple of 32)"),
        )
        .arg(
            Arg::new("sparse")
                .short('s')
                .long("sparse")
                .default_value("3")
                .value_parser(value_parser!(usize))
                .help("Dimensions updated per k-mer"),
        )
}
/// Execute the to-hv command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();
    let dim = *args.get_one::<usize>("dim").unwrap();
    let sparse = *args.get_one::<usize>("sparse").unwrap();
    anyhow::ensure!(
        dim > 0 && dim.is_multiple_of(32),
        "--dim must be a positive multiple of 32"
    );
    anyhow::ensure!(sparse > 0, "--sparse must be positive");

    let mut reader = pgr::reader(infile)?;
    let idx = pgr::libs::pgi::PgiIndex::read(&mut reader)?;
    let hv = pgr::libs::pgi::to_hv::index_to_hv(&idx, dim, sparse);
    let name = pgr::libs::io::get_basename(infile).unwrap_or_else(|| infile.clone());
    let mut writer = pgr::writer(outfile)?;
    pgr::libs::pgi::to_hv::write_hv(
        &mut writer,
        &name,
        idx.k,
        dim,
        sparse,
        idx.n_unique() as usize,
        &hv,
    )?;
    log::info!("wrote {}-dim sparse hypervector to {}", dim, outfile);
    Ok(())
}
