use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::phylo::tree::Tree;
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;

/// Build the clap subcommand for reroot.
pub fn make_subcommand() -> Command {
    Command::new("reroot")
        .about("Reroots a tree at a specified node or the longest branch")
        .after_help(
            r###"
Reroots a phylogenetic tree on a specific branch or node.

Notes:
* Target Selection:
    * Default: If NO nodes are specified, reroots at the midpoint of the longest branch.
    * Specified Nodes (`--node` / `-n`):
        * Reroots on the edge leading to the Lowest Common Ancestor (LCA) of specified nodes.
        * Treats specified nodes as the "ingroup" and the rest as the "outgroup".
        * The new root is placed at the midpoint of the parent edge of the LCA.
    * Lax Mode (`--lax` / `-l`):
        * If the LCA of specified nodes is already the root, use the *unspecified* nodes (complement) as the ingroup.
        * Useful for defining an outgroup by exclusion.

* Operations:
    * Reroot (Default): Creates a bifurcating root at the target edge.
    * Deroot (`--deroot` / `-d`): Splices out the ingroup to create a multifurcating root.

* Technical Details:
    * Support Values (`--support-as-labels`):
        * Shifts internal node labels along the rerooting path to maintain split associations.
        * Necessary because rerooting changes edge directions.
    * Topology:
        * The original root's parent edge is merged.
        * Degree-2 nodes created during the process are removed.

Examples:
1. Reroot at the longest branch (default):
   pgr nwk reroot input.nwk

2. Reroot at a specific node (ingroup):
   pgr nwk reroot input.nwk -n Homo

3. Reroot at the LCA of multiple nodes:
   pgr nwk reroot input.nwk -n Homo -n Pan

4. Reroot and preserve support values (internal node labels):
   pgr nwk reroot input.nwk -n Homo --support-as-labels

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::node_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("support_as_labels")
                .long("support-as-labels")
                .action(ArgAction::SetTrue)
                .help("Treat internal node labels as support values and shift them when rerooting"),
        )
        .arg(
            Arg::new("deroot")
                .long("deroot")
                .short('d')
                .action(ArgAction::SetTrue)
                .help("Deroot the tree (create a multifurcating root) (see Notes)"),
        )
        .arg(
            Arg::new("lax")
                .long("lax")
                .short('l')
                .action(ArgAction::SetTrue)
                .help("Lax mode: Use the complement if the specified nodes form the root (see Notes)"),
        )
}

/// Execute the reroot command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let process_support = args.get_flag("support_as_labels");
    let deroot = args.get_flag("deroot");
    let lax = args.get_flag("lax");

    let infile = args.get_one::<String>("infile").unwrap();
    let trees = Tree::from_file(infile)?;

    if trees.len() > 1 {
        log::warn!(
            "file contains {} trees, only the first will be processed",
            trees.len()
        );
    }

    // Process only the first tree for now, or loop if we want to support multi-tree
    // Since arguments are node names, it implies a single tree context or consistent naming.
    // We'll process the first tree to match previous behavior.
    let mut tree = trees
        .into_iter()
        .next()
        .ok_or_else(|| anyhow::anyhow!("no trees found in {}", infile))?;

    if deroot {
        tree.deroot().map_err(|e| anyhow::anyhow!(e))?;
    } else {
        // ids with names
        let id_of: BTreeMap<_, _> = tree.get_name_id();

        // All IDs matched
        let mut ids = BTreeSet::new();
        let user_specified = args.get_many::<String>("node").is_some();
        if let Some(nodes) = args.get_many::<String>("node") {
            for name in nodes {
                if let Some(&id) = id_of.get(name) {
                    ids.insert(id);
                } else {
                    log::warn!("node name not found in tree: {}", name);
                }
            }
        }

        if !ids.is_empty() {
            pgr::libs::phylo::tree::ops::reroot_at_lca(&mut tree, &ids, lax, process_support)?;
        } else if user_specified {
            anyhow::bail!("none of the specified --node names were found in the tree");
        } else {
            pgr::libs::phylo::tree::ops::reroot_at_longest_branch(&mut tree, process_support)?;
        }
    }

    let out_string = tree.to_newick();
    writer.write_all((out_string + "\n").as_ref())?;

    Ok(())
}
