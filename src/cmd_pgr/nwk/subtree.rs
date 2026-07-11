use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

/// Build the clap subcommand for subtree.
pub fn make_subcommand() -> Command {
    Command::new("subtree")
        .about("Extracts a subtree")
        .after_help(
            r###"
Extracts a subtree (clade) rooted at the lowest common ancestor (LCA) of all
provided nodes.

Notes:
* Node selection:
    * `-n`: Select nodes by exact name.
    * `-l`: Select nodes from a name-list file.
    * `-x`: Select nodes by regular expression.
    * If no selection is provided, no output is generated.
* Monophyly check (`-M`):
    * Ensures the subtree defined by the LCA contains ONLY the selected named
      terminal nodes.
    * Useful to verify if a group is monophyletic.
* Condense (`--condense`):
    * Instead of extracting the subtree, it replaces the subtree with a single
      node in the original tree.
    * The new node inherits the edge length of the subtree root.
    * Added annotations: `member=<count>` and `tri=white`.

Examples:
1. Extract subtree for Human and Chimp:
   pgr nwk subtree tree.nwk -n Human -n Chimp

2. Extract subtree for all nodes matching "Homo":
   pgr nwk subtree tree.nwk -x "^Homo"

3. Condense the Hominini clade (LCA of Homo and Pan) into a single node "Hominini":
   pgr nwk subtree tree.nwk -n Homo -n Pan --condense Hominini

4. Check if a group is monophyletic (prints nothing if not):
   pgr nwk subtree tree.nwk -n Human -n Chimp -n Gorilla -M
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::node_arg())
        .arg(crate::cmd_pgr::args::name_list_arg())
        .arg(crate::cmd_pgr::args::regex_arg())
        .arg(crate::cmd_pgr::args::descendants_arg())
        .arg(crate::cmd_pgr::args::monophyly_arg(
            "Only print the subtree when it's monophyletic",
        ))
        .arg(
            Arg::new("condense")
                .long("condense")
                .short('C')
                .num_args(1)
                .help("Condense the subtree into a single node with this name"),
        )
        .arg(
            Arg::new("context")
                .long("context")
                .short('c')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .default_value("0")
                .help("Extend the subtree by N levels above the LCA"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the subtree command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    let is_monophyly = args.get_flag("monophyly");
    let condense_name = args.get_one::<String>("condense");
    let is_condense = condense_name.is_some();

    let infile = args.get_one::<String>("infile").unwrap();
    let trees = Tree::from_file(infile)?;

    if trees.is_empty() {
        return Ok(());
    }

    let mut trees = trees;

    for tree in &mut trees {
        // IDs matching names
        let ids = super::common::match_names(tree, args)?;

        if ids.is_empty() {
            continue;
        }

        // Find LCA
        let ids_vec: Vec<usize> = ids.iter().cloned().collect();
        let mut sub_root_id = tree.get_lca(&ids_vec)?;

        // Monophyly check
        if is_monophyly {
            let descendants_named = tree.get_named_leaves(sub_root_id);

            if ids != descendants_named {
                if is_condense {
                    let out_string = tree.to_newick();
                    writer.write_fmt(format_args!("{}\n", out_string))?;
                }
                continue;
            }
        }

        // Apply context
        let context_levels = *args.get_one::<usize>("context").unwrap();
        for _ in 0..context_levels {
            if let Some(node) = tree.get_node(sub_root_id) {
                if let Some(parent) = node.parent {
                    sub_root_id = parent;
                } else {
                    break;
                }
            }
        }

        if is_condense {
            // condense_name is Some here (is_condense == condense_name.is_some()).
            let name = condense_name.map_or("", |s| s.as_str());
            if name.is_empty() {
                anyhow::bail!("--condense requires a non-empty name argument");
            }

            tree.condense_subtree(sub_root_id, name, ids.len())?;

            let out_string = tree.to_newick();
            writer.write_fmt(format_args!("{}\n", out_string))?;
        } else {
            // Extract subtree
            let out_string = tree.to_newick_subtree(sub_root_id);
            writer.write_fmt(format_args!("{}\n", out_string))?;
        }
    }

    writer.flush()?;
    Ok(())
}
