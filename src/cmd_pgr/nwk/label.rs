use clap::*;
use std::collections::BTreeSet;
use super::utils as nwr;
use pgr::libs::phylo::reader;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("label")
        .about("Labels in the Newick file")
        .after_help(
            r###"
Extracts the tree's labels.

By default, prints all labels that occur in the tree, in the same order as
in the Newick, one per line. Empty labels produce no output.

Notes:
* The `-t` option prints labels on a single line, separated by tabs.
* The `-I` and `-L` options filter out internal or leaf nodes.
* Selection options (`-n`, `-f`, `-r`) can be combined.
* With `-D`, descendants of selected internal nodes are also included.
* Monophyly check (`-M`) verifies if the selected nodes form a monophyletic
  group. It checks terminal nodes against the selection.
* Warning: Duplicate node names may affect selection/monophyly checks.
* Extra columns (`-c`) details:
    * `dup` - duplicate the node name
    * `taxid` - `:T=` field in comment
    * `species` - `:S=` field in comment
    * `full` - full comment

Examples:
1. List all labels:
   pgr nwk label tree.nwk

2. Count leaves:
   pgr nwk label tree.nwk -I | wc -l

3. List specific nodes:
   pgr nwk label tree.nwk -n Human -n Chimp

4. List labels matching regex:
   pgr nwk label tree.nwk -r "^Homo"

5. Check monophyly:
   pgr nwk label tree.nwk -n Human -n Chimp -M

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
            Arg::new("Internal")
                .long("Internal")
                .short('I')
                .action(ArgAction::SetTrue)
                .help("Don't print internal labels"),
        )
        .arg(
            Arg::new("Leaf")
                .long("Leaf")
                .short('L')
                .action(ArgAction::SetTrue)
                .help("Don't print leaf labels"),
        )
        .arg(
            Arg::new("node")
                .long("node")
                .short('n')
                .num_args(1)
                .action(ArgAction::Append)
                .help("Select nodes by exact name"),
        )
        .arg(
            Arg::new("file")
                .long("file")
                .short('f')
                .num_args(1)
                .help("Select nodes from a file"),
        )
        .arg(
            Arg::new("regex")
                .long("regex")
                .short('r')
                .num_args(1)
                .action(ArgAction::Append)
                .help("Select nodes by regular expression (case insensitive)"),
        )
        .arg(
            Arg::new("descendants")
                .long("descendants")
                .short('D')
                .action(ArgAction::SetTrue)
                .help("Include all descendants of selected internal nodes"),
        )
        .arg(
            Arg::new("root")
                .long("root")
                .action(ArgAction::SetTrue)
                .help("Only print the root label"),
        )
        .arg(
            Arg::new("tab")
                .long("tab")
                .short('t')
                .action(ArgAction::SetTrue)
                .help("Print labels on a single line, separated by tab stops"),
        )
        .arg(
            Arg::new("monophyly")
                .long("monophyly")
                .short('M')
                .action(ArgAction::SetTrue)
                .help("Only print the labels when they form a monophyletic subtree"),
        )
        .arg(
            Arg::new("column")
                .long("column")
                .short('c')
                .action(ArgAction::Append)
                .value_parser([
                    builder::PossibleValue::new("dup"),
                    builder::PossibleValue::new("taxid"),
                    builder::PossibleValue::new("species"),
                    builder::PossibleValue::new("full"),
                ])
                .help("Add extra columns to output"),
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
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    let infile = args.get_one::<String>("infile").unwrap();
    let trees = reader::from_file(infile)?;
    if trees.is_empty() {
        return Ok(());
    }
    let mut columns = vec![];
    if args.contains_id("column") {
        for column in args.get_many::<String>("column").unwrap() {
            columns.push(column.to_string());
        }
    }

    for tree in &trees {
        // Handle --root option
        if args.get_flag("root") {
            let root_id = tree.get_root().unwrap();
            let root = tree.get_node(root_id).unwrap();
            if let Some(name) = &root.name {
                if !name.is_empty() {
                    writer.write_fmt(format_args!("{}\n", name)).unwrap();
                }
            }
            continue;
        }

        let is_monophyly = args.get_flag("monophyly");

        //----------------------------
        // Operating
        //----------------------------
        // All IDs matching positions
        let ids_pos = nwr::match_positions(&tree, args);

        // All IDs matching names
        let ids_name = nwr::match_names(&tree, args);

        let ids: BTreeSet<usize> = ids_pos.intersection(&ids_name).cloned().collect();

        // Print nothing if check_monophyly() failed
        let ids_vec: Vec<usize> = ids.iter().cloned().collect();
        if is_monophyly && !tree.is_monophyletic(&ids_vec) {
            continue;
        }

        //----------------------------
        // Output
        //----------------------------
        let tab_sep = args.get_flag("tab");
        let mut collected_labels = Vec::new();

        for id in ids.iter() {
            let node = tree.get_node(*id).unwrap();
            if let Some(x) = node.name.clone() {
                let mut out_string: String = x.clone();
                if !columns.is_empty() {
                    for column in columns.iter() {
                        match column.as_str() {
                            "dup" => out_string += format!("\t{}", x).as_str(),
                            "taxid" => {
                                out_string += format!(
                                    "\t{}",
                                    node.get_property("T").map(|s: &String| s.as_str()).unwrap_or("")
                                )
                                .as_str()
                            }
                            "species" => {
                                out_string += format!(
                                    "\t{}",
                                    node.get_property("S").map(|s: &String| s.as_str()).unwrap_or("")
                                )
                                .as_str()
                            }
                            "full" => {
                                let props = node.properties.as_ref().map(|p: &std::collections::BTreeMap<String, String>| {
                                    p.iter().map(|(k,v)| format!(":{}={}", k, v)).collect::<Vec<String>>()
                                });
                                
                                let mut comment = String::new();
                                if let Some(p) = props {
                                    if !p.is_empty() {
                                        comment = format!("[{}]", p.join(" "));
                                    }
                                }

                                out_string += format!(
                                    "\t{}",
                                    comment
                                )
                                .as_str()
                            }
                            _ => unreachable!(),
                        }
                    }
                }

                if tab_sep {
                    collected_labels.push(out_string);
                } else {
                    writer.write_fmt(format_args!("{}\n", out_string)).unwrap();
                }
            }
        }

        if tab_sep && !collected_labels.is_empty() {
            writer.write_fmt(format_args!("{}\n", collected_labels.join("\t"))).unwrap();
        }
    }

    Ok(())
}
