use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

/// Build the clap subcommand for comment.
pub fn make_subcommand() -> Command {
    Command::new("comment")
        .about("Adds comments to node(s) in a Newick file")
        .after_help(
            r###"
* Comments are in the NHX-like format
    * :key=value

* For nodes with names, set `--node` to the name
* For nodes without names (e.g., internal nodes), set `--lca` to a combination
  of two node names, separated by commas
    * `--lca A,B`

* Set `--string` to add free-form strings

* The following options are used for visualization
    * `--color`, `--label` and `--comment-text` take 1 argument
    * `--dot`, `--bar` and `--rec` take 1 or 0 argument

* Predefined colors for `--color`, `--dot` and `--bar`
    * {red}{RGB}{188,36,46}
    * {black}{RGB}{26,25,25}
    * {grey}{RGB}{129,130,132}
    * {green}{RGB}{32,128,108}
    * {purple}{RGB}{160,90,150}
* Colors for background rectangles `--rec`
    * {LemonChiffon}{RGB}{251, 248, 204}
    * {ChampagnePink}{RGB}{253, 228, 207}
    * {TeaRose}{RGB}{255, 207, 210}
    * {PinkLavender}{RGB}{241, 192, 232}
    * {Mauve}{RGB}{207, 186, 240}
    * {JordyBlue}{RGB}{163, 196, 243}
    * {NonPhotoBlue}{RGB}{144, 219, 244}
    * {ElectricBlue}{RGB}{142, 236, 245}
    * {Aquamarine}{RGB}{152, 245, 225}
    * {Celadon}{RGB}{185, 251, 192}

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::node_arg())
        .arg(crate::cmd_pgr::args::lca_arg())
        .arg(
            Arg::new("string")
                .long("string")
                .num_args(1)
                .help("Free-form strings"),
        )
        .arg(crate::cmd_pgr::args::color_arg(None, "Color of names"))
        .arg(
            Arg::new("label")
                .long("label")
                .num_args(1)
                .help("Add this label to the south west of the node"),
        )
        .arg(
            Arg::new("comment_text")
                .long("comment-text")
                .num_args(1)
                .help("Comment text after names"),
        )
        .arg(
            Arg::new("dot")
                .long("dot")
                .num_args(0..=1)
                .default_missing_value("black")
                .help("Place a dot in the node; value as color"),
        )
        .arg(
            Arg::new("bar")
                .long("bar")
                .num_args(0..=1)
                .default_missing_value("black")
                .help("Place a bar in the middle of the parent edge; value as color"),
        )
        .arg(
            Arg::new("rec")
                .long("rec")
                .num_args(0..=1)
                .default_missing_value("LemonChiffon")
                .help("Place a rectangle in the background of the subtree; value as color"),
        )
        .arg(
            Arg::new("tri")
                .long("tri")
                .num_args(0..=1)
                .default_missing_value("white")
                .help("Place a triangle at the end of the branch; value as color"),
        )
        .arg(
            Arg::new("remove")
                .long("remove")
                .num_args(1)
                .help("Scan all nodes and remove parts of comments matching the regex"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the comment command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    let opt_string = args.get_one::<String>("string");

    let opt_label = args.get_one::<String>("label");
    let opt_color = args.get_one::<String>("color");
    let opt_comment = args.get_one::<String>("comment_text");

    let opt_dot = args.get_one::<String>("dot");
    let opt_bar = args.get_one::<String>("bar");
    let opt_rec = args.get_one::<String>("rec");
    let opt_tri = args.get_one::<String>("tri");

    let infile = args.get_one::<String>("infile").unwrap();
    let mut trees = Tree::from_file(infile)?;

    for tree in &mut trees {
        // ids with names, name => id
        let id_of = tree.get_name_id();

        // all IDs to be modified
        let mut ids = vec![];

        // ids supplied by --node
        if args.contains_id("node") {
            for name in args.get_many::<String>("node").unwrap() {
                if let Some(id) = id_of.get(name) {
                    ids.push(*id);
                } else {
                    log::warn!("node not found: {}", name);
                }
            }
        }

        // ids supplied by --lca
        if args.contains_id("lca") {
            for lca in args.get_many::<String>("lca").unwrap() {
                let (first, last) = super::common::parse_lca_pair(lca)?;

                match (id_of.get(first), id_of.get(last)) {
                    (Some(id1), Some(id2)) => {
                        let id = tree.get_common_ancestor(id1, id2)?;
                        ids.push(id);
                    }
                    _ => {
                        log::warn!("lca name not found in tree: {} / {}", first, last);
                    }
                }
            }
        }

        for id in &ids {
            if let Some(node) = tree.get_node_mut(*id) {
                if let Some(x) = opt_string {
                    node.add_property_from_str(x);
                }

                if let Some(x) = opt_label {
                    node.add_property("label", x);
                }
                if let Some(x) = opt_color {
                    node.add_property("color", x);
                }
                if let Some(x) = opt_comment {
                    node.add_property("comment", x);
                }

                if let Some(x) = opt_dot {
                    node.add_property("dot", x);
                }
                if let Some(x) = opt_bar {
                    node.add_property("bar", x);
                }
                if let Some(x) = opt_rec {
                    node.add_property("rec", x);
                }
                if let Some(x) = opt_tri {
                    node.add_property("tri", x);
                }
            }
        }

        // Remove parts of comments
        if args.contains_id("remove") {
            let pattern = args.get_one::<String>("remove").unwrap();
            pgr::libs::phylo::tree::ops::remove_properties_matching(tree, pattern)?;
        }

        let out_string = tree.to_newick();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    writer.flush()?;
    Ok(())
}
