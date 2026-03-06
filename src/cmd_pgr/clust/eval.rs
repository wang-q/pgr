use clap::{Arg, ArgMatches, Command};
use pgr::libs::clust::eval::{
    ball_hall_score, c_index_score, calinski_harabasz_score, davies_bouldin_score, dunn_score,
    evaluate, gamma_score, load_batch_partitions, load_partition, pbm_score, silhouette_score,
    tau_score, wemmert_gancarski_score, xie_beni_score, Coordinates, DistanceMatrix, Partition,
    PartitionFormat, TreeDistance,
};
use pgr::libs::pairmat::NamedMatrix;
use pgr::libs::phylo::tree::Tree;
use std::fs::File;
use std::io::{self, Write};

pub fn make_subcommand() -> Command {
    Command::new("eval")
        .about("Evaluates clustering quality")
        .after_help(
            r###"Calculates clustering evaluation metrics.

Modes:
1. External Evaluation (Partition vs Partition):
   Compares two clustering partitions (e.g., ground truth vs result).
   Metrics: ARI, AMI, V-Measure.

2. Internal Evaluation (Partition + Matrix):
   Evaluates a single partition using a distance matrix.
   Metrics: Silhouette Coefficient.

3. Batch Evaluation (Long Format):
   Evaluates multiple partitions (e.g. from parameter scan) against a ground truth or using internal metrics.
   Input file must be in 'long' format (Group, Cluster, Sample).

Examples:
1. Compare result with ground truth:
   $ pgr clust eval result.tsv --other other.tsv -o eval.tsv

2. Evaluate result using distance matrix:
   $ pgr clust eval result.tsv --matrix dist.phy

3. Batch evaluation of scan results:
   $ pgr clust eval scan.tsv --format long --matrix dist.phy
"###,
        )
        .arg(
            Arg::new("p1")
                .required(true)
                .index(1)
                .help("Partition file"),
        )
        .arg(
            Arg::new("other")
                .long("other")
                .alias("truth")
                .num_args(1)
                .help("Other partition file (for external evaluation)"),
        )
        .arg(
            Arg::new("matrix")
                .long("matrix")
                .num_args(1)
                .help("Distance matrix file (for internal evaluation: Silhouette)"),
        )
        .arg(
            Arg::new("tree")
                .long("tree")
                .num_args(1)
                .help("Tree file (for internal evaluation: Silhouette, using patristic distance)"),
        )
        .arg(
            Arg::new("coords")
                .long("coords")
                .num_args(1)
                .help("Coordinate matrix file (for internal evaluation: Davies-Bouldin)"),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .value_parser(["cluster", "pair", "long"])
                .default_value("pair")
                .help("Input format for partition files"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
        .arg(
            Arg::new("no-singletons")
                .long("no-singletons")
                .action(clap::ArgAction::SetTrue)
                .help("Exclude true singletons (from Reference/Ground Truth) from evaluation"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let p1_path = matches.get_one::<String>("p1").unwrap();
    let outfile = matches.get_one::<String>("outfile").unwrap();

    let format_str = matches.get_one::<String>("format").unwrap();
    let format: PartitionFormat = format_str.parse().expect("Invalid format");

    let mut writer: Box<dyn Write> = if outfile == "stdout" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(outfile)?)
    };

    let remove_singletons_flag = matches.get_flag("no-singletons");

    if format == PartitionFormat::Long {
        // Batch Mode
        let batches = load_batch_partitions(p1_path)?;

        // Prepare resources
        let p2 = if let Some(p2_path) = matches.get_one::<String>("other") {
            let mut truth = load_partition(p2_path, PartitionFormat::Pair)?;
            if remove_singletons_flag {
                remove_singletons(&mut truth);
            }
            Some(truth)
        } else {
            None
        };

        let dist_provider: Option<Box<dyn DistanceMatrix>> =
            if let Some(matrix_path) = matches.get_one::<String>("matrix") {
                Some(Box::new(NamedMatrix::from_relaxed_phylip(matrix_path)))
            } else if let Some(tree_path) = matches.get_one::<String>("tree") {
                let trees = Tree::from_file(tree_path)?;
                if trees.len() != 1 {
                    anyhow::bail!("Tree file must contain exactly one tree.");
                }
                Some(Box::new(TreeDistance::new(
                    trees.into_iter().next().unwrap(),
                )))
            } else {
                None
            };

        let coords = if let Some(coords_path) = matches.get_one::<String>("coords") {
            Some(Coordinates::from_path(coords_path)?)
        } else {
            None
        };

        if p2.is_none() && dist_provider.is_none() && coords.is_none() {
            anyhow::bail!(
                "Batch mode requires at least one evaluation target: --other/--truth, --matrix, --tree, or --coords."
            );
        }

        // Write Header
        let mut header = vec!["Group"];
        if p2.is_some() {
            header.extend_from_slice(&[
                "ari",
                "ami",
                "homogeneity",
                "completeness",
                "v_measure",
                "fmi",
                "nmi",
                "mi",
                "ri",
                "jaccard",
                "precision",
                "recall",
            ]);
        }
        if dist_provider.is_some() {
            header.push("silhouette");
            header.push("dunn");
            header.push("c_index");
            header.push("gamma");
            header.push("tau");
        }
        if coords.is_some() {
            header.push("davies_bouldin");
            header.push("calinski_harabasz");
            header.push("pbm");
            header.push("ball_hall");
            header.push("xie_beni");
            header.push("wemmert_gancarski");
        }
        writeln!(writer, "{}", header.join("\t"))?;

        // Process batches
        for (group, p1) in batches {
            let mut row = vec![group];

            if let Some(ref truth) = p2 {
                let metrics = evaluate(&p1, truth);
                row.push(format!("{:.6}", metrics.ari));
                row.push(format!("{:.6}", metrics.ami));
                row.push(format!("{:.6}", metrics.homogeneity));
                row.push(format!("{:.6}", metrics.completeness));
                row.push(format!("{:.6}", metrics.v_measure));
                row.push(format!("{:.6}", metrics.fmi));
                row.push(format!("{:.6}", metrics.nmi));
                row.push(format!("{:.6}", metrics.mi));
                row.push(format!("{:.6}", metrics.ri));
                row.push(format!("{:.6}", metrics.jaccard));
                row.push(format!("{:.6}", metrics.precision));
                row.push(format!("{:.6}", metrics.recall));
            }

            if let Some(ref d) = dist_provider {
                let s_score = silhouette_score(&p1, d.as_ref());
                let d_score = dunn_score(&p1, d.as_ref());
                let c_score = c_index_score(&p1, d.as_ref());
                let g_score = gamma_score(&p1, d.as_ref());
                let t_score = tau_score(&p1, d.as_ref());
                row.push(format!("{:.6}", s_score));
                row.push(format!("{:.6}", d_score));
                row.push(format!("{:.6}", c_score));
                row.push(format!("{:.6}", g_score));
                row.push(format!("{:.6}", t_score));
            }

            if let Some(ref c) = coords {
                let db_score = davies_bouldin_score(&p1, c);
                let ch_score = calinski_harabasz_score(&p1, c);
                let pbm = pbm_score(&p1, c);
                let bh = ball_hall_score(&p1, c);
                let xb = xie_beni_score(&p1, c);
                let wg = wemmert_gancarski_score(&p1, c);
                row.push(format!("{:.6}", db_score));
                row.push(format!("{:.6}", ch_score));
                row.push(format!("{:.6}", pbm));
                row.push(format!("{:.6}", bh));
                row.push(format!("{:.6}", xb));
                row.push(format!("{:.6}", wg));
            }

            writeln!(writer, "{}", row.join("\t"))?;
        }

        return Ok(());
    }

    // Single Mode
    let p1 = load_partition(p1_path, format)?;

    if let Some(p2_path) = matches.get_one::<String>("other") {
        let mut p2 = load_partition(p2_path, format)?;
        if remove_singletons_flag {
            remove_singletons(&mut p2);
        }
        let metrics = evaluate(&p1, &p2);

        writeln!(
            writer,
            "ari\tami\thomogeneity\tcompleteness\tv_measure\tfmi\tnmi\tmi\tri\tjaccard\tprecision\trecall"
        )?;
        writeln!(
            writer,
            "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            metrics.ari,
            metrics.ami,
            metrics.homogeneity,
            metrics.completeness,
            metrics.v_measure,
            metrics.fmi,
            metrics.nmi,
            metrics.mi,
            metrics.ri,
            metrics.jaccard,
            metrics.precision,
            metrics.recall
        )?;
    } else if let Some(matrix_path) = matches.get_one::<String>("matrix") {
        let matrix = NamedMatrix::from_relaxed_phylip(matrix_path);
        let s_score = silhouette_score(&p1, &matrix);
        let d_score = dunn_score(&p1, &matrix);
        let c_score = c_index_score(&p1, &matrix);
        let g_score = gamma_score(&p1, &matrix);
        let t_score = tau_score(&p1, &matrix);

        writeln!(writer, "silhouette\tdunn\tc_index\tgamma\ttau")?;
        writeln!(
            writer,
            "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            s_score, d_score, c_score, g_score, t_score
        )?;
    } else if let Some(tree_path) = matches.get_one::<String>("tree") {
        let trees = Tree::from_file(tree_path)?;
        if trees.len() != 1 {
            anyhow::bail!("Tree file must contain exactly one tree.");
        }
        let dist = TreeDistance::new(trees.into_iter().next().unwrap());
        let s_score = silhouette_score(&p1, &dist);
        let d_score = dunn_score(&p1, &dist);
        let c_score = c_index_score(&p1, &dist);
        let g_score = gamma_score(&p1, &dist);
        let t_score = tau_score(&p1, &dist);

        writeln!(writer, "silhouette\tdunn\tc_index\tgamma\ttau")?;
        writeln!(
            writer,
            "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            s_score, d_score, c_score, g_score, t_score
        )?;
    } else if let Some(coords_path) = matches.get_one::<String>("coords") {
        let coords = Coordinates::from_path(coords_path)?;
        let db_score = davies_bouldin_score(&p1, &coords);
        let ch_score = calinski_harabasz_score(&p1, &coords);
        let pbm = pbm_score(&p1, &coords);
        let bh = ball_hall_score(&p1, &coords);
        let xb = xie_beni_score(&p1, &coords);
        let wg = wemmert_gancarski_score(&p1, &coords);

        writeln!(
            writer,
            "davies_bouldin\tcalinski_harabasz\tpbm\tball_hall\txie_beni\twemmert_gancarski"
        )?;
        writeln!(
            writer,
            "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            db_score, ch_score, pbm, bh, xb, wg
        )?;
    } else {
        anyhow::bail!(
            "Either --other/--truth (for external eval), --matrix, --tree, or --coords (for internal eval) must be provided."
        );
    }

    Ok(())
}

fn remove_singletons(partition: &mut Partition) {
    let mut counts = std::collections::HashMap::new();
    for cid in partition.values() {
        *counts.entry(*cid).or_insert(0) += 1;
    }
    partition.retain(|_, cid| counts[cid] > 1);
}
