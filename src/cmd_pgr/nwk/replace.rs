use clap::{ArgMatches, Command};
use pgr::libs::phylo::tree::Tree;
use std::collections::BTreeMap;
use std::io::Write;

/// Build the clap subcommand for replace.
pub fn make_subcommand() -> Command {
    Command::new("replace")
        .about("Replaces node names or comments in a Newick file")
        .after_help(
            r###"
Replaces node names or appends annotations in a Newick file using a TSV file.

Notes:
* `--replace-tsv` is a tab-separated file with 2 or more columns:
  `<original_name> <replacement> [additional_annotations...]`
* The behavior of the 2nd column (`<replacement>`) depends on `--mode`:
    * `label` (default): Replaces the node name. Empty string removes the name.
    * `taxid`:           Appends as NCBI TaxID (`:T=<replacement>`) in NHX.
    * `species`:         Appends as species name (`:S=<replacement>`) in NHX.
    * `asis`:            Appends verbatim to comments (e.g. `key=val` or `tag`).
* Columns 3+ are ALWAYS appended to the node's comments/properties.
  Key-value pairs (e.g., `color=red`) are stored as properties.
  Simple tags (e.g., `highlight`) are stored as keys with empty values.

Examples:
1. Basic renaming of nodes:
   pgr nwk replace input.nwk --replace-tsv names.tsv > output.nwk

2. Add species and color annotations:
   pgr nwk replace input.nwk --replace-tsv annotations.tsv --mode species

3. Remove node names (2nd column is empty):
   pgr nwk replace input.nwk --replace-tsv remove.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::replace_tsv_arg())
        .arg(crate::cmd_pgr::args::internal_arg())
        .arg(crate::cmd_pgr::args::leaf_arg())
        .arg(crate::cmd_pgr::args::mode_arg(
            "label",
            &["label", "taxid", "species", "asis"],
            "Where we place the replaces",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the replace command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let infile = args.get_one::<String>("infile").unwrap();
    let mode = args.get_one::<String>("mode").unwrap();

    let skip_internal = args.get_flag("internal");
    let skip_leaf = args.get_flag("leaf");

    let mut replace_of: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let rfile = args.get_one::<String>("replace_tsv").unwrap();
    for line in pgr::read_lines(rfile)? {
        let parts: Vec<_> = line.split('\t').collect();

        if parts.len() < 2 {
            log::warn!("skipping malformed line in replace file: {}", line);
            continue;
        }
        let name = parts[0].to_string();
        let replaces = parts
            .iter()
            .skip(1)
            .map(|e| e.to_string())
            .collect::<Vec<String>>();
        replace_of.insert(name, replaces);
    }

    let mut trees = Tree::from_file(infile)?;

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
                let first = replaces
                    .first()
                    .ok_or_else(|| anyhow::anyhow!("no replace values"))?
                    .to_string();
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
                    other => anyhow::bail!("unknown property mode: {}", other),
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
