use clap::{Arg, ArgAction, ArgGroup, ArgMatches, Command, Id};
use pgr::libs::phylo::tree::{algo, Tree};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("order")
        .about("Orders nodes in a Newick file")
        .after_help(
            r###"
Sorts the children of each node without changing the topology.

Notes:
* Traverses the entire tree in a breadth-first order.
* `--alphanumeric` and `--num-descendants` can be enabled at the same time; sorted first by `--alphanumeric` and then by `--num-descendants`.
* `--name-list` is processed before `--alphanumeric` and `--num-descendants`.
* Sort orders:
    * `--name-list`: By a list of names in the file, one name per line.
    * `--alphanumeric`/`--alphanumeric-rev`: By alphanumeric order of labels.
    * `--num-descendants`/`--num-descendants-rev`: By number of descendants (ladderize).
    * `--deladderize`: Alternate sort direction at each level.

Examples:
1. Sort by number of descendants (ladderize):
   pgr nwk order tree.nwk --num-descendants

2. Sort by alphanumeric order of labels:
   pgr nwk order tree.nwk --alphanumeric

3. Sort by a list of names:
   pgr nwk order tree.nwk --name-list names.txt

4. Sort by alphanumeric order, then by number of descendants (reverse):
   pgr nwk order tree.nwk --alphanumeric --num-descendants-rev

5. De-ladderize (alternate sort direction):
   pgr nwk order tree.nwk --deladderize

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(
            Arg::new("num_descendants")
                .long("num-descendants")
                .action(ArgAction::SetTrue)
                .help("By number of descendants"),
        )
        .arg(
            Arg::new("num_descendants_rev")
                .long("num-descendants-rev")
                .action(ArgAction::SetTrue)
                .help("By number of descendants, reversely"),
        )
        .group(ArgGroup::new("number-of-descendants").args(["num_descendants", "num_descendants_rev"]))
        .arg(
            Arg::new("alphanumeric")
                .long("alphanumeric")
                .action(ArgAction::SetTrue)
                .help("By alphanumeric order of labels"),
        )
        .arg(
            Arg::new("alphanumeric_rev")
                .long("alphanumeric-rev")
                .action(ArgAction::SetTrue)
                .help("By alphanumeric order of labels, reversely"),
        )
        .group(ArgGroup::new("alphanumeric-order").args(["alphanumeric", "alphanumeric_rev"]))
        .arg(
            Arg::new("deladderize")
                .long("deladderize")
                .alias("dl")
                .action(ArgAction::SetTrue)
                .help("De-ladderize (alternate) the tree"),
        )
        .arg(crate::cmd_pgr::args::name_list_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let opt_nd = match args.get_one::<Id>("number-of-descendants") {
        None => "",
        Some(x) => x.as_str(),
    };
    let opt_an = match args.get_one::<Id>("alphanumeric-order") {
        None => "",
        Some(x) => x.as_str(),
    };

    let infile = args.get_one::<String>("infile").unwrap();
    let mut trees = Tree::from_file(infile)?;

    let mut names = vec![];
    if args.contains_id("name_list") {
        let list_file = args.get_one::<String>("name_list").unwrap();
        names = pgr::libs::io::read_names::<Vec<String>>(list_file)?;
    }

    let is_deladderize = args.get_flag("deladderize");

    // Default behavior: if no specific sort order is requested, use alphanumeric
    let default_an = names.is_empty() && opt_an.is_empty() && opt_nd.is_empty() && !is_deladderize;

    for tree in &mut trees {
        if !names.is_empty() {
            algo::sort_by_list(tree, &names);
        }
        if default_an || !opt_an.is_empty() {
            algo::sort_by_name(tree, opt_an == "alphanumeric_rev");
        }
        if !opt_nd.is_empty() {
            algo::ladderize(tree, opt_nd == "num_descendants_rev");
        }
        if is_deladderize {
            algo::deladderize(tree);
        }

        let out_string = tree.to_newick();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}
