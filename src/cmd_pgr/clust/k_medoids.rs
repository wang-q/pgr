use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("k-medoids")
        .about("K-Medoids clustering")
        .visible_alias("km")
        .after_help(
            r###"
K-Medoids clustering algorithm (PAM/Lloyd-like).

Note: The input file should contain pairwise distances (lower is better), NOT similarities.

Output formats:
    * cluster: Each line contains points of one cluster.
    * pair: Each line contains a (representative point, cluster member) pair.

Note:
For the 'pair' format, the representative point is the medoid (point with minimum sum of distances to other cluster members).
If there are ties, the alphabetically first member is chosen.
"###,
        )
        .arg(
            crate::cmd_pgr::args::infile_arg_required_with_help(
                "Input file containing pairwise distances in .tsv format",
            ),
        )
        .arg(crate::cmd_pgr::args::k_arg().required(true))
        .arg(crate::cmd_pgr::args::format_arg())
        .arg(crate::cmd_pgr::args::same_arg("0.0"))
        .arg(crate::cmd_pgr::args::missing_arg("1.0"))
        .arg(
            Arg::new("runs")
                .long("runs")
                .num_args(1)
                .default_value("10")
                .value_parser(value_parser!(usize))
                .help("Number of random initializations"),
        )
        .arg(crate::cmd_pgr::args::max_iter_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // 1. Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let opt_k = *args.get_one::<usize>("k").unwrap();
    let opt_format = args.get_one::<String>("clust_format").unwrap();
    let opt_same = *args.get_one::<f32>("same").unwrap();
    let opt_missing = *args.get_one::<f32>("missing").unwrap();
    let runs = *args.get_one::<usize>("runs").unwrap();
    let max_iter = *args.get_one::<usize>("max_iter").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let mut writer = pgr::writer(outfile)?;

    //----------------------------
    // 2. Load Matrix
    //----------------------------
    let (sm, names): (pgr::libs::pairmat::ScoringMatrix<f32>, Vec<String>) =
        pgr::libs::pairmat::ScoringMatrix::from_pair_scores(infile, opt_same, opt_missing)?;

    //----------------------------
    // 3. Clustering
    //----------------------------
    let kmedoids = pgr::libs::clust::k_medoids::KMedoids::new(opt_k, max_iter, runs);
    let mut clusters = kmedoids.perform_clustering(&sm);

    //----------------------------
    // 4. Output
    //----------------------------
    let out =
        pgr::libs::clust::format::format_flat_clusters(&mut clusters, &names, opt_format, |c| {
            pgr::libs::clust::medoid::find_medoid(&sm, c, false)
        })?;
    writer.write_all(out.as_bytes())?;

    Ok(())
}
