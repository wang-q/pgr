use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for dbscan.
pub fn make_subcommand() -> Command {
    Command::new("dbscan")
        .about("Clusters entries via DBSCAN")
        .after_help(
            r###"
Density-based spatial clustering of applications with noise (DBSCAN).

Note: The input file should contain pairwise distances (lower is better), NOT similarities.

Output formats:
    * cluster: Each line contains points of one cluster.
    * pair: Each line contains a (representative point, cluster member) pair.

Note:
The representative point is selected by --rep and affects both output formats:
    * For 'pair' format: it is the first column (representative ID).
    * For 'cluster' format: it is placed in the first column.
    * medoid (default): point with minimum sum of distances to other cluster members.
    * first: alphabetically first member of the cluster.

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input file containing pairwise distances in .tsv format",
        ))
        .arg(crate::cmd_pgr::args::format_arg())
        .arg(crate::cmd_pgr::args::flat_rep_arg())
        .arg(crate::cmd_pgr::args::same_arg("0.0"))
        .arg(crate::cmd_pgr::args::missing_arg("1.0"))
        .arg(crate::cmd_pgr::args::eps_arg())
        .arg(crate::cmd_pgr::args::min_points_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the dbscan command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    // 1. Args
    let infile = args
        .get_one::<String>("infile")
        .ok_or_else(|| anyhow::anyhow!("missing required argument: infile"))?;

    // Remaining arguments have clap default values, so unwrap is safe.
    let opt_format = args.get_one::<String>("clust_format").unwrap();
    let opt_rep = args.get_one::<String>("flat_rep").unwrap().as_str();
    let opt_same = *args.get_one::<f32>("same").unwrap();
    let opt_missing = *args.get_one::<f32>("missing").unwrap();
    let opt_eps = *args.get_one::<f32>("eps").unwrap();
    let opt_min_points = *args.get_one::<usize>("min_points").unwrap();

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    // 2. Load Matrix

    // Load matrix from pairwise distances
    let (matrix, names) =
        pgr::libs::pairmat::ScoringMatrix::from_pair_scores(infile, opt_same, opt_missing)?;

    // 3. Clustering
    let mut dbscan = pgr::libs::clust::dbscan::Dbscan::new(opt_eps, opt_min_points);
    dbscan.perform_clustering(&matrix);
    let mut clusters = dbscan.results_cluster();

    // 4. Output
    let out = if opt_rep == "first" {
        pgr::libs::clust::format::format_flat_clusters(&mut clusters, &names, opt_format, |c| {
            c.first().copied()
        })?
    } else {
        pgr::libs::clust::format::format_flat_clusters(&mut clusters, &names, opt_format, |c| {
            pgr::libs::clust::medoid::find_medoid(&matrix, c, false)
        })?
    };
    writer.write_all(out.as_bytes())?;

    writer.flush()?;
    Ok(())
}
