use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("cc")
        .about("Connected components clustering (ignoring weights)")
        .after_help(
            r###"
Ignores scores and writes all connected components.

Output formats:
    * cluster: Each line contains points of one cluster.
    * pair: Each line contains a (representative point, cluster member) pair.

Note:
    For the 'pair' format, the representative point is the alphabetically first member of the cluster.

"###,
        )
        .arg(
            crate::cmd_pgr::args::infile_arg_required_with_help(
                "Input file containing pairwise relations (weights ignored) in .tsv format",
            ),
        )
        .arg(crate::cmd_pgr::args::format_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // 1. Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let opt_format = args.get_one::<String>("clust_format").unwrap();
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    //----------------------------
    // 2. Load Graph & Clustering
    //----------------------------
    let reader = pgr::reader(infile)?;
    let (names_vec, mut scc) = pgr::libs::clust::connected_components(reader)?;

    //----------------------------
    // 3. Output
    //----------------------------
    let out =
        pgr::libs::clust::format::format_flat_clusters(&mut scc, &names_vec, opt_format, |c| {
            c.first().copied()
        })?;
    writer.write_all(out.as_bytes())?;

    Ok(())
}
