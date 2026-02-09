use clap::*;
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("comment")
        .about("Add comments to node(s) in a Newick file")
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
    * `--color`, `--label` and `--comment` take 1 argument
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
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input filename. [stdin] for standard input"),
        )
        .arg(
            Arg::new("node")
                .long("node")
                .short('n')
                .num_args(1)
                .action(ArgAction::Append)
                .help("Node name"),
        )
        .arg(
            Arg::new("lca")
                .long("lca")
                .short('l')
                .num_args(1)
                .action(ArgAction::Append)
                .help("Lowest common ancestor of two nodes"),
        )
        .arg(
            Arg::new("string")
                .long("string")
                .short('s')
                .num_args(1)
                .help("Free-form strings"),
        )
        .arg(
            Arg::new("color")
                .long("color")
                .num_args(1)
                .help("Color of names"),
        )
        .arg(
            Arg::new("label")
                .long("label")
                .num_args(1)
                .help("Add this label to the south west of the node"),
        )
        .arg(
            Arg::new("comment")
                .long("comment")
                .num_args(1)
                .help("comment text after names"),
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
                .short('r')
                .num_args(1)
                .help("Scan all nodes and remove parts of comments matching the regex"),
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
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());

    let opt_string = args.get_one::<String>("string");

    let opt_label = args.get_one::<String>("label");
    let opt_color = args.get_one::<String>("color");
    let opt_comment = args.get_one::<String>("comment");

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
                }
            }
        }

        // ids supplied by --lca
        if args.contains_id("lca") {
            for lca in args.get_many::<String>("lca").unwrap() {
                let parts = lca.split(',').map(|e| e.to_string()).collect::<Vec<_>>();
                if parts.len() != 2 {
                    continue;
                }

                if let (Some(id1), Some(id2)) = (id_of.get(&parts[0]), id_of.get(&parts[1])) {
                    if let Ok(id) = tree.get_common_ancestor(id1, id2) {
                        ids.push(id);
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

        //----------------------------
        // Remove parts of comments
        //----------------------------
        // ids matched with --remove
        if args.contains_id("remove") {
            let regex = args.get_one::<String>("remove").unwrap();
            let re = regex::RegexBuilder::new(regex)
                .case_insensitive(true)
                .unicode(false)
                .build()
                .unwrap();

            // Iterate all nodes
            for i in 0..tree.len() {
                if let Some(node) = tree.get_node_mut(i) {
                    if let Some(props) = &mut node.properties {
                        let mut to_remove = vec![];
                        for (k, v) in props.iter() {
                            let entry = if v.is_empty() {
                                format!("{}", k)
                            } else {
                                format!("{}={}", k, v)
                            };

                            if re.is_match(&entry) {
                                to_remove.push(k.clone());
                            }
                        }

                        for k in to_remove {
                            props.remove(&k);
                        }
                    }
                }
            }
        }

        //----------------------------
        // Output
        //----------------------------
        let out_string = tree.to_newick();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}
