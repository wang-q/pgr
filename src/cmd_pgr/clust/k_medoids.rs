use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for k-medoids.
pub fn make_subcommand() -> Command {
    Command::new("k-medoids")
        .about("Clusters entries via K-Medoids")
        .visible_alias("km")
        .after_help(
            r###"
K-Medoids clustering algorithm (PAM/Lloyd-like).

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
        .arg(crate::cmd_pgr::args::k_arg().required(true))
        .arg(crate::cmd_pgr::args::format_arg())
        .arg(crate::cmd_pgr::args::flat_rep_arg())
        .arg(crate::cmd_pgr::args::same_arg("0.0"))
        .arg(crate::cmd_pgr::args::missing_arg("1.0"))
        .arg(crate::cmd_pgr::args::runs_arg())
        .arg(crate::cmd_pgr::args::max_iter_arg())
        .arg(crate::cmd_pgr::args::seed_arg(
            None,
            None,
            "Random seed for reproducible initialization",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the k-medoids command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    // 1. Args
    let infile = args
        .get_one::<String>("infile")
        .ok_or_else(|| anyhow::anyhow!("missing required argument: infile"))?;
    let opt_k = *args
        .get_one::<usize>("k")
        .ok_or_else(|| anyhow::anyhow!("missing required argument: k"))?;
    // Remaining arguments have clap default values, so unwrap is safe.
    let opt_format = args.get_one::<String>("clust_format").unwrap();
    let opt_rep = args.get_one::<String>("flat_rep").unwrap().as_str();
    let opt_same = *args.get_one::<f32>("same").unwrap();
    let opt_missing = *args.get_one::<f32>("missing").unwrap();
    let runs = *args.get_one::<usize>("runs").unwrap();
    let max_iter = *args.get_one::<usize>("max_iter").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    // 2. Load Matrix
    let (sm, names): (pgr::libs::pairmat::ScoringMatrix<f32>, Vec<String>) =
        pgr::libs::pairmat::ScoringMatrix::from_pair_scores(infile, opt_same, opt_missing)?;

    // 3. Clustering
    let mut kmedoids = pgr::libs::clust::k_medoids::KMedoids::new(opt_k, max_iter, runs);
    if let Some(&seed) = args.get_one::<u64>("seed") {
        kmedoids = kmedoids.with_seed(seed);
    }
    let mut clusters = kmedoids.perform_clustering(&sm);

    // 4. Output
    let out = if opt_rep == "first" {
        pgr::libs::clust::format::format_flat_clusters(&mut clusters, &names, opt_format, |c| {
            c.first().copied()
        })?
    } else {
        pgr::libs::clust::format::format_flat_clusters(&mut clusters, &names, opt_format, |c| {
            pgr::libs::clust::medoid::find_medoid(&sm, c, false)
        })?
    };
    writer.write_all(out.as_bytes())?;

    writer.flush()?;
    Ok(())
}
