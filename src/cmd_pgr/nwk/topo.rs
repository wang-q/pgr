use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("topo")
        .about("Manipulates tree topology and attributes")
        .after_help(
            r###"
Modifies tree topology by optionally removing branch lengths, comments, or labels.

Notes:
* By default, branch lengths and comments are REMOVED.
* Use `--bl` to KEEP branch lengths.
* Use `--comment` to KEEP comments.
* Use `-I` to REMOVE internal labels.
* Use `-L` to REMOVE leaf labels.

Examples:
1. Topology only (remove lengths and comments):
   pgr nwk topo tree.nwk

2. Keep branch lengths but remove comments:
   pgr nwk topo tree.nwk --bl

3. Remove internal node labels (topology only):
   pgr nwk topo tree.nwk -I
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::bl_arg())
        .arg(
            Arg::new("comment")
                .long("comment")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Keep comments"),
        )
        .arg(crate::cmd_pgr::args::internal_arg())
        .arg(crate::cmd_pgr::args::leaf_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let is_bl = args.get_flag("bl");
    let is_comment = args.get_flag("comment");
    let skip_internal = args.get_flag("internal");
    let skip_leaf = args.get_flag("leaf");

    let infile = args.get_one::<String>("infile").unwrap();
    let trees = Tree::from_file(infile)?;

    for mut tree in trees {
        if let Some(root) = tree.get_root() {
            // We use levelorder to traverse all nodes safely
            let ids = tree.levelorder(&root).map_err(anyhow::Error::msg)?;

            for id in ids {
                if let Some(node) = tree.get_node_mut(id) {
                    if !is_bl {
                        node.length = None;
                    }
                    if !is_comment {
                        node.properties = None;
                    }
                    if node.is_leaf() && skip_leaf {
                        node.name = None;
                    }
                    if !node.is_leaf() && skip_internal {
                        node.name = None;
                    }
                }
            }
        }

        let out_string = tree.to_newick();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}
