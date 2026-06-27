use clap::*;
use pgr::libs::paf::index::PafIndex;
use std::collections::HashSet;
use std::fs;
use std::io::BufRead;

pub fn make_subcommand() -> Command {
    Command::new("query")
        .about("Query PAF index for coordinate projection")
        .after_help(
            r###"
Queries a PAF file or saved index for intervals overlapping a target
region and projects them to query coordinates via CIGAR.

Accepts either a PAF file (built on-the-fly) or a .paf.idx index
(loaded from disk, instant startup).

Two modes:
* Default: single-hop projection — finds all PAF records whose target
  interval overlaps the query region and lifts coordinates to the
  corresponding query sequence.
* --transitive: multi-hop BFS traversal — iteratively projects through
  intermediate sequences up to --max-depth hops.

Output formats:
* Default: query_name, query_start, query_end, target_name, target_start, target_end
* --bed: BED6 format with target annotation in name column
* --paf: standard PAF 12 columns + tags

Notes:
* Input PAF files should contain cg:Z: tags for accurate projection
* Reads from stdin if input file is 'stdin'

Examples:
1. Single-hop projection from a PAF file:
   pgr paf query alignments.paf chr1:1000-5000

2. Query from a saved index (faster):
   pgr paf query alignments.paf.idx chr1:1000-5000

3. Transitive BFS with filters:
   pgr paf query alignments.paf chr1:1000-5000 --transitive --min-identity 0.8 --bed

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input PAF file or .paf.idx index to query"),
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
        .arg(
            Arg::new("min_identity")
                .long("min-identity")
                .num_args(1)
                .default_value("0.0")
                .value_parser(clap::value_parser!(f64))
                .help("Minimum gap-compressed identity (0.0-1.0, default: 0.0)"),
        )
        .arg(
            Arg::new("min_output_len")
                .long("min-output-len")
                .num_args(1)
                .default_value("0")
                .value_parser(clap::value_parser!(i32))
                .help("Minimum output interval length (default: 0 = no filter)"),
        )
        .arg(
            Arg::new("merge_distance")
                .long("merge-distance")
                .num_args(1)
                .default_value("0")
                .value_parser(clap::value_parser!(i32))
                .help("Merge adjacent output intervals within this distance (default: 0 = off)"),
        )
        .arg(
            Arg::new("subset_list")
                .long("subset-sequence-list")
                .num_args(1)
                .help("File with sequence names to include (one per line)"),
        )
        .arg(
            Arg::new("bed")
                .long("bed")
                .num_args(0)
                .help("Output in BED6 format"),
        )
        .arg(
            Arg::new("paf")
                .long("paf")
                .num_args(0)
                .conflicts_with("bed")
                .help("Output in PAF format (12 columns + tags)"),
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

fn load_subset(path: &str) -> anyhow::Result<HashSet<String>> {
    let f = fs::File::open(path)?;
    let mut set = HashSet::new();
    for line in std::io::BufReader::new(f).lines() {
        let line = line?;
        let name = line.trim().to_string();
        if !name.is_empty() {
            set.insert(name);
        }
    }
    Ok(set)
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let region_str = args.get_one::<String>("region").unwrap();
    let transitive = args.get_flag("transitive");
    let max_depth = *args.get_one::<u16>("max_depth").unwrap();
    let min_len = *args.get_one::<i32>("min_len").unwrap();
    let min_dist = *args.get_one::<i32>("min_dist").unwrap();
    let min_identity = *args.get_one::<f64>("min_identity").unwrap();
    let min_output_len = *args.get_one::<i32>("min_output_len").unwrap();
    let merge_distance = *args.get_one::<i32>("merge_distance").unwrap();
    let bed = args.get_flag("bed");
    let paf = args.get_flag("paf");

    let (target_name, start, end) = parse_region(region_str)?;

    let idx = if infile.ends_with(".paf.idx") {
        eprintln!("Loading index from {infile}...");
        PafIndex::load(infile)?
    } else {
        eprintln!("Building index from {infile}...");
        let reader = pgr::reader(infile);
        PafIndex::build(reader)?
    };

    eprintln!(
        "  sequences: {}, targets: {}",
        idx.names.len(),
        idx.num_targets()
    );

    let target_id = match idx.name_to_id(target_name) {
        Some(id) => id,
        None => anyhow::bail!("target '{}' not found in index", target_name),
    };

    let mut results = if transitive {
        idx.query_transitive_bfs(
            target_id,
            start,
            end,
            max_depth,
            min_len,
            min_dist,
            min_identity,
            min_output_len,
            merge_distance,
        )
    } else {
        idx.query(target_id, start, end, min_identity, min_output_len)
    };

    // Subset filter
    if let Some(list_path) = args.get_one::<String>("subset_list") {
        let subset = load_subset(list_path)?;
        results.retain(|(qid, _, _)| {
            let name = idx.id_to_name(*qid).unwrap_or("");
            subset.contains(name)
        });
    }

    if results.is_empty() {
        eprintln!("No results found.");
        return Ok(());
    }

    if bed {
        output_bed(&idx, &results);
    } else if paf {
        output_paf(&idx, &results);
    } else {
        output_default(&idx, &results);
    }

    Ok(())
}

fn output_default(
    idx: &PafIndex,
    results: &[(u32, coitrees::Interval<u32>, coitrees::Interval<u32>)],
) {
    for (query_id, q_iv, t_iv) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            qname, q_iv.first, q_iv.last, tname, t_iv.first, t_iv.last
        );
    }
}

fn output_bed(idx: &PafIndex, results: &[(u32, coitrees::Interval<u32>, coitrees::Interval<u32>)]) {
    for (query_id, q_iv, t_iv) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");
        let name = format!("{}:{}-{}", tname, t_iv.first, t_iv.last);
        println!("{}\t{}\t{}\t{}\t0\t.", qname, q_iv.first, q_iv.last, name);
    }
}

fn output_paf(idx: &PafIndex, results: &[(u32, coitrees::Interval<u32>, coitrees::Interval<u32>)]) {
    for (query_id, q_iv, t_iv) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");
        let block_len = (q_iv.last - q_iv.first).max(1) as u32;
        println!(
            "{}\t0\t{}\t{}\t+\t{}\t0\t{}\t{}\t0\t{}\t255\tgi:f:1.0",
            qname, q_iv.first, q_iv.last, tname, t_iv.first, t_iv.last, block_len
        );
    }
}
