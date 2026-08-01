//! `pgr dist pgi` — deterministic distance between two .pgi indexes.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for pgi.
pub fn make_subcommand() -> Command {
    Command::new("pgi")
        .about("Computes distances between two .pgi indexes by k-mer merge")
        .after_help(
            r###"
Computes deterministic Jaccard / containment / Mash distances between two
genome indexes by merging their sorted k-mer tables. Both indexes must use
identical sampling parameters (k, smer, window), checked against the headers.

Output (tab-separated):
    <idx1> <idx2> <total1> <total2> <inter> <union> <mash> <jaccard> <containment>

Examples:
1. Distance between two indexes:
   pgr dist pgi a.pgi b.pgi
"###,
        )
        .arg(
            Arg::new("idx1")
                .index(1)
                .required(true)
                .help("First .pgi index"),
        )
        .arg(
            Arg::new("idx2")
                .index(2)
                .required(true)
                .help("Second .pgi index"),
        )
}
/// Execute the pgi command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let idx1 = args.get_one::<String>("idx1").unwrap();
    let idx2 = args.get_one::<String>("idx2").unwrap();

    let mut r1 = pgr::reader(idx1)?;
    let mut r2 = pgr::reader(idx2)?;
    let a = pgr::libs::pgi::PgiIndex::read(&mut r1)?;
    let b = pgr::libs::pgi::PgiIndex::read(&mut r2)?;
    let d = pgr::libs::pgi::dist::dist_between(&a, &b)?;

    let n1 = pgr::libs::io::get_basename(idx1).unwrap_or_else(|| idx1.clone());
    let n2 = pgr::libs::io::get_basename(idx2).unwrap_or_else(|| idx2.clone());
    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}",
        n1, n2, d.total1, d.total2, d.inter, d.union, d.mash, d.jaccard, d.containment
    );
    Ok(())
}
