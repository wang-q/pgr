use clap::*;
use pgr::libs::paf::index::PafIndex;

pub fn make_subcommand() -> Command {
    Command::new("query")
        .about("Query PAF index for coordinate projection")
        .after_help(
            r###"
Queries a PAF file for intervals overlapping a target region and projects
them to query coordinates via CIGAR.

Two modes:
* Default: single-hop projection — finds all PAF records whose target
  interval overlaps the query region and lifts coordinates to the
  corresponding query sequence.
* --transitive: multi-hop BFS traversal — iteratively projects through
  intermediate sequences up to --max-depth hops.

Output format:
* Tab-delimited: query_name, query_start, query_end, target_name, target_start, target_end

Notes:
* Input PAF files should contain cg:Z: tags for accurate projection
* Reads from stdin if input file is 'stdin'

Examples:
1. Single-hop projection:
   pgr paf query alignments.paf chr1:1000-5000

2. Transitive BFS up to 3 hops:
   pgr paf query alignments.paf chr1:1000-5000 --transitive --max-depth 3

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input PAF file to query"),
        )
        .arg(
            Arg::new("region")
                .required(true)
                .index(2)
                .help("Target region to query (e.g. chr1:1000-5000)"),
        )
        .arg(
            Arg::new("transitive")
                .long("transitive")
                .short('t')
                .num_args(0)
                .help("Enable transitive BFS traversal"),
        )
        .arg(
            Arg::new("max_depth")
                .long("max-depth")
                .num_args(1)
                .default_value("2")
                .value_parser(clap::value_parser!(u16))
                .help("Maximum BFS depth (0 = unlimited, default: 2)"),
        )
        .arg(
            Arg::new("min_len")
                .long("min-len")
                .num_args(1)
                .default_value("10")
                .value_parser(clap::value_parser!(i32))
                .help("Minimum interval length to propagate (default: 10)"),
        )
        .arg(
            Arg::new("min_dist")
                .long("min-dist")
                .num_args(1)
                .default_value("10")
                .value_parser(clap::value_parser!(i32))
                .help("Minimum distance to merge adjacent intervals (default: 10)"),
        )
}

fn parse_region(s: &str) -> anyhow::Result<(&str, i32, i32)> {
    let parts: Vec<&str> = s.split(':').collect();
    anyhow::ensure!(
        parts.len() == 2,
        "invalid region '{s}': expected name:start-end"
    );
    let name = parts[0];
    let range: Vec<&str> = parts[1].split('-').collect();
    anyhow::ensure!(range.len() == 2, "invalid region '{s}': expected start-end");
    let start: i32 = range[0].parse()?;
    let end: i32 = range[1].parse()?;
    Ok((name, start, end))
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let region_str = args.get_one::<String>("region").unwrap();
    let transitive = args.get_flag("transitive");
    let max_depth = *args.get_one::<u16>("max_depth").unwrap();
    let min_len = *args.get_one::<i32>("min_len").unwrap();
    let min_dist = *args.get_one::<i32>("min_dist").unwrap();

    let (target_name, start, end) = parse_region(region_str)?;

    eprintln!("Building index from {infile}...");
    let reader = pgr::reader(infile);
    let idx = PafIndex::build(reader)?;
    eprintln!(
        "  sequences: {}, targets: {}",
        idx.names.len(),
        idx.num_targets()
    );

    let target_id = match idx.name_to_id(target_name) {
        Some(id) => id,
        None => {
            anyhow::bail!("target '{}' not found in index", target_name);
        }
    };

    let results = if transitive {
        idx.query_transitive_bfs(target_id, start, end, max_depth, min_len, min_dist)
    } else {
        idx.query(target_id, start, end)
    };

    if results.is_empty() {
        eprintln!("No results found.");
    }

    for (query_id, q_iv, t_iv) in &results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            qname, q_iv.first, q_iv.last, tname, t_iv.first, t_iv.last
        );
    }

    Ok(())
}
