use clap::{ArgMatches, Command};
use pgr::libs::phylo::tree::algo;
use pgr::libs::phylo::tree::query as nwr;
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("prune")
        .about("Removes nodes from a Newick file")
        .after_help(
            r###"
Removes nodes from the Newick tree based on provided labels or patterns.

Notes:
* Target nodes can be specified by name (`--node`), name-list (`--name-list`), or regex (`--regex`).
* With `--invert`, specified nodes (along with their ancestors and descendants) are kept, and everything else is removed.
* Topology changes:
    * If a node removal leaves its parent with only one child, the parent is collapsed (spliced out).
    * Internal nodes that lose all children are also removed.

Examples:
1. Remove specific nodes by name:
   pgr nwk prune input.nwk -n Homo -n Pan

2. Remove nodes using a list in a file:
   pgr nwk prune input.nwk -l remove.txt

3. Keep a clade (e.g., Hominidae) and remove everything else (Invert mode):
   pgr nwk prune input.nwk -i -n Hominidae

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::node_arg())
        .arg(crate::cmd_pgr::args::name_list_arg())
        .arg(crate::cmd_pgr::args::regex_arg())
        .arg(crate::cmd_pgr::args::descendants_arg())
        .arg(crate::cmd_pgr::args::invert_arg_with_help(
            "Invert pruning: keep specified nodes, their ancestors and descendants",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let infile = args.get_one::<String>("infile").unwrap();
    let trees = Tree::from_file(infile)?;

    for mut tree in trees {
        let target_ids = nwr::match_names(&tree, args)?;

        let to_remove: Vec<_> = if args.get_flag("invert") {
            let keep = algo::compute_keep_set(&tree, target_ids.iter().copied());
            match tree.get_root() {
                Some(root) => tree
                    .levelorder(&root)
                    .unwrap_or_default()
                    .into_iter()
                    .filter(|id| !keep.contains(id))
                    .collect(),
                None => Vec::new(),
            }
        } else {
            target_ids.into_iter().collect()
        };

        algo::prune_nodes(&mut tree, to_remove);

        let out_string = tree.to_newick();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}
