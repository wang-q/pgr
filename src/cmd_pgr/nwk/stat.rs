use clap::*;
use pgr::libs::phylo::reader;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Statistics about the Newick file")
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
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let infile = args.get_one::<String>("infile").unwrap();

    let tree = reader::from_file(infile);

    let mut n_edge_w_len = 0;
    let mut n_edge_wo_len = 0;
    let mut n_node = 0;
    let mut n_leaf = 0;
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

    writer.write_fmt(format_args!("Type\t{}\n", tree_type))?;
    writer.write_fmt(format_args!("nodes\t{}\n", n_node))?;
    writer.write_fmt(format_args!("leaves\t{}\n", n_leaf))?;
    writer.write_fmt(format_args!("leaf labels\t{}\n", n_leaf_label))?;
    writer.write_fmt(format_args!("internal labels\t{}\n", n_internal_label))?;

    Ok(())
}
