use clap::*;
use pgr::libs::phylo::reader;
use std::collections::BTreeMap;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("replace")
        .about("Replace node names/comments in a Newick file")
        .after_help(
            r###"
Replace node names or append annotations in a Newick file using a TSV file.

Notes:
* <replace.tsv> is a tab-separated file with 2 or more columns:
  <original_name> <replacement> [additional_annotations...]
* The behavior of the 2nd column (<replacement>) depends on `--mode`:
  * label (default): Replaces the node name. Empty string removes the name.
  * taxid:           Appends as NCBI TaxID (`:T=<replacement>`) in NHX.
  * species:         Appends as species name (`:S=<replacement>`) in NHX.
  * asis:            Appends verbatim to comments (e.g. `key=val` or `tag`).
* Columns 3+ are ALWAYS appended to the node's comments/properties.
  Key-value pairs (e.g., `color=red`) are stored as properties.
  Simple tags (e.g., `highlight`) are stored as keys with empty values.

Examples:
1. Basic renaming of nodes:
   pgr nwk replace input.nwk names.tsv > output.nwk

2. Add species and color annotations:
   pgr nwk replace input.nwk annotations.tsv --mode species

3. Remove node names (2nd column is empty):
   pgr nwk replace input.nwk remove.tsv

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
            Arg::new("replace.tsv")
                .required(true)
                .num_args(1..)
                .index(2)
                .help("Path to replace.tsv"),
        )
        .arg(
            Arg::new("Internal")
                .long("Internal")
                .short('I')
                .action(ArgAction::SetTrue)
                .help("Skip internal labels"),
        )
        .arg(
            Arg::new("Leaf")
                .long("Leaf")
                .short('L')
                .action(ArgAction::SetTrue)
                .help("Skip leaf labels"),
        )
        .arg(
            Arg::new("mode")
                .long("mode")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("label"),
                    builder::PossibleValue::new("taxid"),
                    builder::PossibleValue::new("species"),
                    builder::PossibleValue::new("asis"),
                ])
                .default_value("label")
                .help("Where we place the replaces"),
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
    let mode = args.get_one::<String>("mode").unwrap();

    let skip_internal = args.get_flag("Internal");
    let skip_leaf = args.get_flag("Leaf");

    let mut replace_of: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for rfile in args.get_many::<String>("replace.tsv").unwrap() {
        for line in intspan::read_lines(rfile) {
            let parts: Vec<_> = line.split('\t').collect();

            if parts.len() < 2 {
                continue;
            } else {
                let name = parts.first().unwrap().to_string();
                let replaces = parts
                    .iter()
                    .skip(1)
                    .map(|e| e.to_string())
                    .collect::<Vec<String>>();
                replace_of.insert(name.to_string(), replaces);
            }
        }
    }

    let mut trees = reader::from_file(infile)?;

    for tree in &mut trees {
        let id_of = tree.get_name_id();

        // Collect modifications to avoid mutable borrow issues
        let mut to_modify = Vec::new();

        for (name, id) in &id_of {
            if let Some(replaces) = replace_of.get(name) {
                // Check filters
                if let Some(node) = tree.get_node(*id) {
                    let is_leaf = node.is_leaf();
                    if skip_internal && !is_leaf {
                        continue;
                    }
                    if skip_leaf && is_leaf {
                        continue;
                    }

                    to_modify.push((*id, replaces.clone()));
                }
            }
        }

        // Apply modifications
        for (id, replaces) in to_modify {
            if let Some(node) = tree.get_node_mut(id) {
                let first = replaces.first().unwrap().to_string();
                match mode.as_str() {
                    "label" => node.set_name(first),
                    "taxid" => node.add_property("T", first),
                    "species" => node.add_property("S", first),
                    "asis" => {
                        if first.contains('=') {
                            node.add_property_from_str(&first);
                        } else {
                            node.add_property(&first, "");
                        }
                    }
                    _ => unreachable!(),
                }

                replaces.iter().skip(1).for_each(|e| {
                    if e.contains('=') {
                        node.add_property_from_str(e);
                    } else {
                        node.add_property(e, "");
                    }
                });
            }
        }

        let out_string = tree.to_newick();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}
