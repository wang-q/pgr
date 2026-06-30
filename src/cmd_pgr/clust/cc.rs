use clap::*;
use std::io::BufRead;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("cc")
        .about("Connected components clustering (ignoring weights)")
        .after_help(
            r###"
Ignores scores and writes all connected components.

Output formats:
    * cluster: Each line contains points of one cluster.
    * pair: Each line contains a (representative point, cluster member) pair.

Note:
    For the 'pair' format, the representative point is the alphabetically first member of the cluster.

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input file containing pairwise relations (weights ignored) in .tsv format"),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("cluster"),
                    builder::PossibleValue::new("pair"),
                ])
                .default_value("cluster")
                .help("Output format for clustering results"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // 1. Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let opt_format = args.get_one::<String>("format").unwrap();
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap())?;

    //----------------------------
    // 2. Load Graph
    //----------------------------
    let mut names = indexmap::IndexSet::new();

    let mut graph = petgraph::graphmap::UnGraphMap::<_, ()>::new();

    let reader = pgr::reader(infile)?;
    for line in reader.lines().map_while(Result::ok) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 2 {
            names.insert(fields[0].to_string());
            names.insert(fields[1].to_string());
        }

        graph.add_edge(
            names.get_index_of(fields[0]).unwrap(),
            names.get_index_of(fields[1]).unwrap(),
            (),
        );
    }

    //----------------------------
    // 3. Clustering
    //----------------------------
    let mut scc = petgraph::algo::tarjan_scc(&graph);

    //----------------------------
    // 4. Output
    //----------------------------
    let names_vec: Vec<String> = names.iter().cloned().collect();
    let out =
        pgr::libs::clust::format::format_flat_clusters(&mut scc, &names_vec, opt_format, |c| {
            c.first().copied()
        })?;
    writer.write_all(out.as_bytes())?;

    Ok(())
}
