use clap::*;
use pgr::libs::phylo::tree::Tree;
use std::collections::BTreeMap;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("rename")
        .about("Renames nodes in a Newick file")
        .after_help(
            r###"
Renames nodes in the Newick tree.

Notes:
* For nodes with names, set `--node` to the name.
* For nodes without names (e.g., internal nodes), set `--lca` to a combination
  of two node names, separated by commas.
* The sum of nodes supplied by `--node` and `--lca` should be equal to the number of `--rename`.
* Do not use `--node` and `--lca` alternately.

* This command is not designed to replace names in large batches, but for modifying small amounts
  of data, and therefore does not provide the ability to read a mapping file.
    * `pgr nwk replace` does this kind of jobs.
    * Or use other tools, such as `sed` or `perl`, to accomplish such tasks.

Examples:
1. Rename a named node:
   pgr nwk rename tests/newick/catarrhini.nwk --node Homo --rename Human

2. Rename an internal node via LCA (Hominini is LCA of Homo and Pan):
   pgr nwk rename tests/newick/catarrhini.nwk --lca Homo,Pan --rename CladeX

3. Rename multiple nodes:
   pgr nwk rename tests/newick/catarrhini.nwk \
       --node Homo --rename Human \
       --lca Homo,Pan --rename CladeX
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::node_arg())
        .arg(crate::cmd_pgr::args::lca_arg())
        .arg(
            Arg::new("rename")
                .long("rename")
                .num_args(1)
                .required(true)
                .action(ArgAction::Append)
                .help("New name"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let mut names = vec![];
    if args.contains_id("node") {
        for name in args.get_many::<String>("node").unwrap() {
            names.push(name.to_string());
        }
    }

    let mut lcas = vec![];
    if args.contains_id("lca") {
        for lca in args.get_many::<String>("lca").unwrap() {
            lcas.push(lca.to_string());
        }
    }

    let mut renames = vec![];
    for rename in args.get_many::<String>("rename").unwrap() {
        renames.push(rename.to_string());
    }

    // discard the unnecessary ones
    // make sure renames.len() >= names.len() + lcas.len()
    if names.len() > renames.len() {
        let unnecessary = names.len() - renames.len();
        names.truncate(names.len() - unnecessary);
        // All lcas are unnecessary
        lcas.clear();
    } else if names.len() + lcas.len() > renames.len() {
        let unnecessary = names.len() + lcas.len() - renames.len();
        lcas.truncate(lcas.len() - unnecessary);
    }
    let len_names = names.len();

    let infile = args.get_one::<String>("infile").unwrap();
    let mut trees = Tree::from_file(infile)?;

    for tree in &mut trees {
        // ids with names
        let id_of: BTreeMap<_, _> = tree.get_name_id();

        // all IDs to be modified
        let mut rename_of: BTreeMap<_, _> = BTreeMap::new();

        // ids supplied by --node
        for (i, name) in names.iter().enumerate() {
            if id_of.contains_key(name) {
                let id = id_of.get(name).unwrap();
                let rename = renames.get(i).unwrap();
                rename_of.insert(*id, rename.to_string());
            }
        }

        // ids supplied by --lca
        for (i, lca) in lcas.iter().enumerate() {
            let parts = lca.split(',').map(|e| e.to_string()).collect::<Vec<_>>();
            if parts.len() != 2 {
                continue;
            }

            if parts.iter().all(|e| id_of.contains_key(e)) {
                let id1 = id_of.get(parts.first().unwrap()).unwrap();
                let id2 = id_of.get(parts.last().unwrap()).unwrap();

                if let Ok(x) = tree.get_common_ancestor(id1, id2) {
                    let rename = renames.get(len_names + i).unwrap();
                    rename_of.insert(x, rename.to_string());
                }
            }
        }

        for (k, v) in &rename_of {
            if let Some(node) = tree.get_node_mut(*k) {
                node.name = Some(v.to_string());
            }
        }

        let out_string = tree.to_newick();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}
