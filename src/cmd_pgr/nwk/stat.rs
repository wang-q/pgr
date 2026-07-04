use clap::{Arg, ArgMatches, Command};
use pgr::libs::phylo::tree::{stat, Tree};

/// Build the clap subcommand for stat.
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Prints statistics about trees")
        .after_help(
            r###"
Prints information about the trees in the input.

Notes:
* Input format:
    * Newick trees filename or `stdin`.

* Output format:
    * Key-value pairs (TSV, default):
      Type	cladogram
      nodes	18
      leaves	11
      rooted	Yes
      ...
      cherries	5
      sackin	46
      colless	9

    * Tab-separated values (--style line):
      Type	nodes	leaves	rooted	dichotomies	leaf labels	internal labels	cherries	sackin	colless
      cladogram	18	11	Yes	5	11	0	5	46	9

Examples:
1. Default statistics:
   pgr nwk stat data/catarrhini.nwk

2. Output to file:
   pgr nwk stat data/catarrhini.nwk -o stats.tsv
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("style")
                .long("style")
                .value_parser(["col", "line"])
                .default_value("col")
                .help("Output style. [col] for key-value pairs, [line] for TSV"),
        )
}

/// Execute the stat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let infile = args.get_one::<String>("infile").unwrap();
    let style = args.get_one::<String>("style").unwrap();

    let trees = Tree::from_file(infile)?;

    if style == "line" {
        writer.write_fmt(format_args!(
            "Type\tnodes\tleaves\trooted\tdichotomies\tleaf labels\tinternal labels\tcherries\tsackin\tcolless\n"
        ))?;
    }

    for tree in trees {
        let s = stat::tree_summary(&tree);
        let is_rooted = if s.is_rooted { "Yes" } else { "No" };
        let sackin_str = match s.sackin {
            Some(v) => v.to_string(),
            None => "-".to_string(),
        };
        let colless_str = match s.colless {
            Some(v) => v.to_string(),
            None => "-".to_string(),
        };
        let tree_type = s.tree_type.as_str();

        if style == "line" {
            writer.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                tree_type,
                s.nodes,
                s.leaves,
                is_rooted,
                s.dichotomies,
                s.leaf_labels,
                s.internal_labels,
                s.cherries,
                sackin_str,
                colless_str
            ))?;
        } else {
            writer.write_fmt(format_args!("Type\t{}\n", tree_type))?;
            writer.write_fmt(format_args!("nodes\t{}\n", s.nodes))?;
            writer.write_fmt(format_args!("leaves\t{}\n", s.leaves))?;
            writer.write_fmt(format_args!("rooted\t{}\n", is_rooted))?;
            writer.write_fmt(format_args!("dichotomies\t{}\n", s.dichotomies))?;
            writer.write_fmt(format_args!("leaf labels\t{}\n", s.leaf_labels))?;
            writer.write_fmt(format_args!("internal labels\t{}\n", s.internal_labels))?;
            writer.write_fmt(format_args!("cherries\t{}\n", s.cherries))?;
            writer.write_fmt(format_args!("sackin\t{}\n", sackin_str))?;
            writer.write_fmt(format_args!("colless\t{}\n", colless_str))?;
        }
    }

    Ok(())
}
