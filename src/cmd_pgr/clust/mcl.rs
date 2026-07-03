use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("mcl")
        .about("Markov Clustering Algorithm (MCL)")
        .after_help(
            r###"
MCL is a fast and scalable unsupervised cluster algorithm for graphs (also known as networks) based on simulation of (stochastic) flow in graphs.

It is particularly useful for clustering protein interaction networks or similarity networks.

Note: The input file should contain similarity scores (higher is better), NOT distances.

Output formats:
    * cluster: Each line contains points of one cluster.
    * pair: Each line contains a (representative point, cluster member) pair.

Note:
For the 'pair' format, the representative point is the medoid (point with maximum sum of similarities to other cluster members).
If there are ties, the alphabetically first member is chosen.

Reference:
Stijn van Dongen, Graph Clustering by Flow Simulation. PhD thesis, University of Utrecht, May 2000.
"###,
        )
        .arg(
            crate::cmd_pgr::args::infile_arg_required_with_help(
                "Input file containing pairwise similarities (edge weights) in .tsv format",
            ),
        )
        .arg(crate::cmd_pgr::args::format_arg())
        .arg(crate::cmd_pgr::args::same_arg("1.0"))
        .arg(crate::cmd_pgr::args::missing_arg("0.0"))
        .arg(
            Arg::new("inflation")
                .long("inflation")
                .short('I')
                .num_args(1)
                .default_value("2.0")
                .value_parser(value_parser!(f64))
                .help("Inflation parameter. Controls the granularity of clusters. Higher values = tighter/more clusters."),
        )
        .arg(
            Arg::new("prune")
                .long("prune")
                .num_args(1)
                .default_value("1e-5")
                .value_parser(value_parser!(f64))
                .help("Pruning threshold. Matrix entries smaller than this will be set to zero."),
        )
        .arg(crate::cmd_pgr::args::max_iter_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // 1. Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let opt_format = args.get_one::<String>("format").unwrap();
    let opt_same = *args.get_one::<f32>("same").unwrap();
    let opt_missing = *args.get_one::<f32>("missing").unwrap();
    let inflation = *args.get_one::<f64>("inflation").unwrap();
    let prune = *args.get_one::<f64>("prune").unwrap();
    let max_iter = *args.get_one::<usize>("max_iter").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let mut writer = pgr::writer(outfile)?;

    //----------------------------
    // 2. Load Matrix
    //----------------------------
    // ScoringMatrix::from_pair_scores is only implemented for f32
    let (sm, names) =
        pgr::libs::pairmat::ScoringMatrix::<f32>::from_pair_scores(infile, opt_same, opt_missing)?;

    //----------------------------
    // 3. Clustering
    //----------------------------
    let mut mcl = pgr::libs::clust::mcl::Mcl::new(inflation);
    mcl.set_prune_limit(prune);
    mcl.set_max_iter(max_iter);
    let mut clusters = mcl.perform_clustering(&sm);

    //----------------------------
    // 4. Output
    //----------------------------
    let out =
        pgr::libs::clust::format::format_flat_clusters(&mut clusters, &names, opt_format, |c| {
            pgr::libs::clust::medoid::find_medoid(&sm, c, true)
        })?;
    writer.write_all(out.as_bytes())?;

    Ok(())
}
