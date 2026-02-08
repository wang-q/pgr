use clap::*;
use pgr::libs::phylo::reader;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Prints statistics about trees")
        .after_help(
            r###"
Prints information about the trees in the input.

Input format:
* Newick trees filename or 'stdin'

Output format:
* Key-value pairs (TSV, default):
  Type	cladogram
  nodes	18
  leaves	11
  ...

* Tab-separated values (--style line):
  Type	nodes	leaves	dichotomies	leaf labels	internal labels
  cladogram	18	11	5	11	0

Examples:
1. Default statistics:
   pgr nwk stat data/catarrhini.nw

2. Output to file:
   pgr nwk stat data/catarrhini.nw -o stats.tsv
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input filename. [stdin] for standard input"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
        .arg(
            Arg::new("style")
                .long("style")
                .value_parser(["col", "line"])
                .default_value("col")
                .help("Output style. [col] for key-value pairs, [line] for TSV"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let infile = args.get_one::<String>("infile").unwrap();
    let style = args.get_one::<String>("style").unwrap();

    let trees = reader::from_file(infile);

    if style == "line" {
        writer.write_fmt(format_args!(
            "Type\tnodes\tleaves\tdichotomies\tleaf labels\tinternal labels\n"
        ))?;
    }

    for tree in trees {
        let mut n_edge_w_len = 0;
        let mut n_edge_wo_len = 0;
        let mut n_node = 0;
        let mut n_leaf = 0;
        let mut n_dichotomies = 0;
        let mut n_leaf_label = 0;
        let mut n_internal_label = 0;

        if let Some(root) = tree.get_root() {
            if let Ok(nodes) = tree.preorder(&root) {
                 for id in nodes {
                    let node = tree.get_node(id).unwrap();
                    n_node += 1;
                    if node.is_leaf() {
                        n_leaf += 1;
                    }

                    if node.children.len() == 2 {
                        n_dichotomies += 1;
                    }

                    if node.name.is_some() {
                        if node.is_leaf() {
                            n_leaf_label += 1;
                        } else {
                            n_internal_label += 1;
                        }
                    }
                    
                    if node.length.is_some() {
                        n_edge_w_len += 1;
                    } else {
                        n_edge_wo_len += 1;
                    }
                }
            }
        }

        let tree_type = if n_edge_wo_len == n_node {
            "cladogram"
        } else if n_edge_w_len == n_node || n_edge_w_len == n_node - 1 {
            "phylogram"
        } else {
            "neither"
        };

        if style == "line" {
            writer.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                tree_type, n_node, n_leaf, n_dichotomies, n_leaf_label, n_internal_label
            ))?;
        } else {
            writer.write_fmt(format_args!("Type\t{}\n", tree_type))?;
            writer.write_fmt(format_args!("nodes\t{}\n", n_node))?;
            writer.write_fmt(format_args!("leaves\t{}\n", n_leaf))?;
            writer.write_fmt(format_args!("dichotomies\t{}\n", n_dichotomies))?;
            writer.write_fmt(format_args!("leaf labels\t{}\n", n_leaf_label))?;
            writer.write_fmt(format_args!("internal labels\t{}\n", n_internal_label))?;
        }
    }

    Ok(())
}
