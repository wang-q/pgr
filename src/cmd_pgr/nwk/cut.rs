use clap::*;
use pgr::libs::phylo::tree::{cut, Tree};
use std::io::Write;

pub fn make_subcommand() -> Command {
    Command::new("cut")
        .about("Cut a tree into clusters")
        .after_help(
            r###"
Cuts the tree into clusters based on various criteria.

Criteria:
* `--k <K>`: Cut into K clusters (top-down split by height).
* `--height <H>`: Cut at specific height (max distance to leaves).
* `--root-dist <D>`: Cut at specific distance from root.
* `--max-clade <T>`: TreeCluster style (max pairwise distance in clade <= T).

Output formats:
    * cluster: Each line contains points of one cluster. The first point is the representative.
    * pair: Each line contains a (representative point, cluster member) pair.

Note:
The representative point is determined by `--rep` (applies to both 'cluster' and 'pair' formats):
* root (default): Member closest to root (alphabetical tie-break).
* medoid: Member with min sum of distances to others (alphabetical tie-break).
* first: Alphabetically first member.

Examples:
1. Cut into 5 clusters:
   pgr nwk cut tree.nwk --k 5

2. Cut at height 0.5:
   pgr nwk cut tree.nwk --height 0.5

3. Cut where max pairwise distance in cluster <= 0.1:
   pgr nwk cut tree.nwk --max-clade 0.1
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input Newick file"),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("cluster"),
                    builder::PossibleValue::new("pair"),
                ])
                .default_value("cluster")
                .help("Output format for clustering results"),
        )
        .arg(
            Arg::new("k")
                .long("k")
                .short('k')
                .value_parser(value_parser!(usize))
                .help("Number of clusters"),
        )
        .arg(
            Arg::new("height")
                .long("height")
                .value_parser(value_parser!(f64))
                .help("Cut at specific height (max distance to leaves)"),
        )
        .arg(
            Arg::new("root-dist")
                .long("root-dist")
                .value_parser(value_parser!(f64))
                .help("Cut at specific distance from root"),
        )
        .arg(
            Arg::new("max-clade")
                .long("max-clade")
                .value_parser(value_parser!(f64))
                .help("Max pairwise distance in cluster threshold"),
        )
        .arg(
            Arg::new("rep")
                .long("rep")
                .value_parser([
                    builder::PossibleValue::new("root"),
                    builder::PossibleValue::new("first"),
                    builder::PossibleValue::new("medoid"),
                ])
                .default_value("root")
                .help("Representative selection method"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .default_value("stdout")
                .help("Output file path"),
        )
        .arg(
            Arg::new("inconsistent")
                .long("inconsistent")
                .value_parser(value_parser!(f64))
                .help("Cut by inconsistent coefficient threshold"),
        )
        .arg(
            Arg::new("deep")
                .long("deep")
                .value_parser(value_parser!(usize))
                .default_value("2")
                .help("Depth for inconsistent coefficient calculation (default: 2)"),
        )
        .group(
            ArgGroup::new("method")
                .args(["k", "height", "root-dist", "max-clade", "inconsistent"])
                .required(true),
        )
}

fn compute_root_distances(
    tree: &Tree,
) -> std::collections::HashMap<pgr::libs::phylo::node::NodeId, f64> {
    let mut dists = std::collections::HashMap::new();
    if let Some(root) = tree.get_root() {
        let mut stack = vec![(root, 0.0)];
        while let Some((node_id, d)) = stack.pop() {
            dists.insert(node_id, d);
            if let Some(node) = tree.get_node(node_id) {
                for &child in &node.children {
                    let len = tree.get_node(child).and_then(|n| n.length).unwrap_or(0.0);
                    stack.push((child, d + len));
                }
            }
        }
    }
    dists
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let infile = matches.get_one::<String>("infile").unwrap();
    let outfile = matches.get_one::<String>("outfile").unwrap();
    let format = matches.get_one::<String>("format").unwrap();
    let rep_method = matches.get_one::<String>("rep").unwrap().as_str();
    let deep = *matches.get_one::<usize>("deep").unwrap();

    let trees = Tree::from_file(infile)?;

    let mut writer = pgr::writer(outfile);

    let method = if let Some(&k) = matches.get_one::<usize>("k") {
        cut::Method::K(k)
    } else if let Some(&h) = matches.get_one::<f64>("height") {
        cut::Method::Height(h)
    } else if let Some(&d) = matches.get_one::<f64>("root-dist") {
        cut::Method::RootDist(d)
    } else if let Some(&t) = matches.get_one::<f64>("max-clade") {
        cut::Method::MaxClade(t)
    } else if let Some(&t) = matches.get_one::<f64>("inconsistent") {
        cut::Method::Inconsistent(t, deep)
    } else {
        unreachable!("ArgGroup requires one method");
    };

    for tree in trees.iter() {
        let partition = cut::cut(tree, method).map_err(|e| anyhow::anyhow!(e))?;
        let root_dists = compute_root_distances(tree);

        let clusters_map = partition.get_clusters();

        // Convert NodeId to names and group into Vec<Vec<(NodeId, String)>>
        // We need NodeId for representative selection
        let mut clusters: Vec<Vec<(pgr::libs::phylo::node::NodeId, String)>> = Vec::new();

        for members in clusters_map.values() {
            let mut member_info = Vec::new();
            for &mid in members {
                if let Some(node) = tree.get_node(mid) {
                    let name = if let Some(name) = &node.name {
                        name.clone()
                    } else {
                        format!("Leaf_{}", mid)
                    };
                    member_info.push((mid, name));
                }
            }
            // Sort members within each cluster alphabetically
            member_info.sort_by(|a, b| a.1.cmp(&b.1));
            clusters.push(member_info);
        }

        // Sort clusters: first by size (descending), then by first member name
        clusters.sort_by(|a, b| match b.len().cmp(&a.len()) {
            std::cmp::Ordering::Equal => {
                let name_a = a.first().map(|s| s.1.as_str()).unwrap_or("");
                let name_b = b.first().map(|s| s.1.as_str()).unwrap_or("");
                name_a.cmp(name_b)
            }
            other => other,
        });

        // Output
        let find_rep =
            |c: &Vec<(pgr::libs::phylo::node::NodeId, String)>| -> (Option<String>, usize) {
                match rep_method {
                    "first" => {
                        if let Some(first) = c.first() {
                            (Some(first.1.clone()), 0)
                        } else {
                            (None, 0)
                        }
                    }
                    "root" => {
                        if let Some((idx, rep)) =
                            c.iter().enumerate().min_by(|(_, (id1, _)), (_, (id2, _))| {
                                let d1 = root_dists.get(id1).unwrap_or(&f64::MAX);
                                let d2 = root_dists.get(id2).unwrap_or(&f64::MAX);
                                d1.partial_cmp(d2).unwrap_or(std::cmp::Ordering::Equal)
                            })
                        {
                            (Some(rep.1.clone()), idx)
                        } else {
                            (None, 0)
                        }
                    }
                    "medoid" => {
                        if c.len() <= 1 {
                            if let Some(first) = c.first() {
                                (Some(first.1.clone()), 0)
                            } else {
                                (None, 0)
                            }
                        } else {
                            let ids: Vec<_> = c.iter().map(|(id, _)| *id).collect();
                            let mut min_sum_dist = f64::MAX;
                            let mut best_idx = 0;

                            for i in 0..ids.len() {
                                let mut current_sum = 0.0;
                                for j in 0..ids.len() {
                                    if i == j {
                                        continue;
                                    }
                                    let dist = pgr::libs::phylo::tree::query::get_distance(
                                        tree, &ids[i], &ids[j],
                                    )
                                    .map(|(d, _)| d)
                                    .unwrap_or(f64::MAX);
                                    current_sum += dist;
                                }
                                if current_sum < min_sum_dist {
                                    min_sum_dist = current_sum;
                                    best_idx = i;
                                }
                            }
                            (Some(c[best_idx].1.clone()), best_idx)
                        }
                    }
                    _ => unreachable!(),
                }
            };

        match format.as_str() {
            "cluster" => {
                for c in clusters {
                    let (best_rep_name, best_rep_idx) = find_rep(&c);

                    if let Some(_) = best_rep_name {
                        let mut names: Vec<&str> =
                            c.iter().map(|(_, name)| name.as_str()).collect();
                        if best_rep_idx != 0 {
                            names.swap(0, best_rep_idx);
                            // After swap, re-sort the rest to maintain determinism
                            names[1..].sort();
                        }

                        writer.write_fmt(format_args!("{}\n", names.join("\t")))?;
                    }
                }
            }
            "pair" => {
                for c in clusters {
                    let (best_rep_name, _) = find_rep(&c);

                    if let Some(rep_name) = best_rep_name {
                        // For pair format, we might want to ensure the rep is listed first?
                        // But pair format is just rep\tmember. Order of lines usually follows member order.
                        // c is sorted alphabetically.
                        for (_, member_name) in &c {
                            writer.write_fmt(format_args!("{}\t{}\n", rep_name, member_name))?;
                        }
                    }
                }
            }
            _ => unreachable!(),
        }
    }

    Ok(())
}
