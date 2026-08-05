//! `pgr pgi stat` — show index parameters and sizes.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for stat.
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Shows .pgi index statistics")
        .after_help(
            r###"
Prints the index parameters (k, syncmer, contigs, unique k-mers, positions)
and the file size for a .pgi file.

Examples:
1. Inspect an index:
   pgr pgi stat genome.pgi
"###,
        )
        .arg(
            Arg::new("infile")
                .index(1)
                .required(true)
                .help(".pgi index file"),
        )
}
/// Execute the stat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let idx = pgr::libs::pgi::PgiMmap::open(std::path::Path::new(infile))?;
    let size = std::fs::metadata(infile)?.len();
    println!("File: {}", infile);
    println!("K-mer size: {}", idx.k());
    println!("Syncmer: {}/{}", idx.smer(), idx.window());
    println!("Contigs: {}", idx.contigs().len());
    println!("Unique k-mers: {}", pgr::libs::pgi::count_unique(&idx));
    println!("Positions: {}", idx.n_records());
    println!("File size: {} bytes", size);
    Ok(())
}
