use clap::*;
use pgr::libs::phylo::{algo, tree::Tree};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("order")
        .about("Order nodes in a Newick file")
        .after_help(
            r###"
Sorts the children of each node without changing the topology.

Notes:
* Traverses the entire tree in a breadth-first order.
* `--an` and `--nd` can be enabled at the same time; sorted first by `--an` and then by `--nd`.
* `--list` is processed before `--an` and `--nd`.
* Sort orders:
    * `--list`: By a list of names in the file, one name per line.
    * `--an`/`--anr`: By alphanumeric order of labels.
    * `--nd`/`--ndr`: By number of descendants (ladderize).
    * `--deladderize`: Alternate sort direction at each level.

Examples:
1. Sort by number of descendants (ladderize):
   pgr nwk order tree.nwk --nd

2. Sort by alphanumeric order of labels:
   pgr nwk order tree.nwk --an

3. Sort by a list of names:
   pgr nwk order tree.nwk --list names.txt

4. Sort by alphanumeric order, then by number of descendants (reverse):
   pgr nwk order tree.nwk --an --ndr

5. De-ladderize (alternate sort direction):
   pgr nwk order tree.nwk --deladderize

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input filename. [stdin] for standard input"),
        )
        .arg(arg!(--nd  "By number of descendants"))
        .arg(arg!(--ndr "By number of descendants, reversely"))
        .group(ArgGroup::new("number-of-descendants").args(["nd", "ndr"]))
        .arg(arg!(--an  "By alphanumeric order of labels"))
        .arg(arg!(--anr "By alphanumeric order of labels, reversely"))
        .group(ArgGroup::new("alphanumeric").args(["an", "anr"]))
        .arg(
            Arg::new("deladderize")
                .long("deladderize")
                .alias("dl")
                .action(ArgAction::SetTrue)
                .help("De-ladderize (alternate) the tree"),
        )
        .arg(
            Arg::new("list")
                .long("list")
                .short('l')
                .num_args(1)
                .help("Order by a list of names in the file"),
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
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    let opt_nd = match args.get_one::<Id>("number-of-descendants") {
        None => "",
        Some(x) => x.as_str(),
    };
    let opt_an = match args.get_one::<Id>("alphanumeric") {
        None => "",
        Some(x) => x.as_str(),
    };

    let infile = args.get_one::<String>("infile").unwrap();
    let mut trees = Tree::from_file(infile)?;

    let mut names = vec![];
    if args.contains_id("list") {
        let list_file = args.get_one::<String>("list").unwrap();
        names = intspan::read_first_column(list_file);
    }

    let is_deladderize = args.get_flag("deladderize");

    // Default behavior: if no specific sort order is requested, use alphanumeric
    let default_an = names.is_empty() && opt_an.is_empty() && opt_nd.is_empty() && !is_deladderize;

    for tree in &mut trees {
        if !names.is_empty() {
            algo::sort_by_list(tree, &names);
        }
        if default_an || !opt_an.is_empty() {
            algo::sort_by_name(tree, opt_an == "anr");
        }
        if !opt_nd.is_empty() {
            algo::ladderize(tree, opt_nd == "ndr");
        }
        if is_deladderize {
            algo::deladderize(tree);
        }

        let out_string = tree.to_newick();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}
