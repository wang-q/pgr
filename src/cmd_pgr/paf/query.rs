use clap::*;
use pgr::libs::chain::record::read_chains;
use pgr::libs::paf::index::{PafIndex, QueryResult};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::BufRead;

/// Add common query arguments to a clap Command.
/// Shared by `paf query`, `paf to-bed`, and `paf to-maf`.
pub fn add_query_args(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("infile")
            .required(true)
            .index(1)
            .help("Input PAF file or .paf.idx index to query"),
    )
    .arg(
        Arg::new("region")
            .index(2)
            .help("Target region to query (e.g. chr1:1000-5000)"),
    )
    .arg(
        Arg::new("bed_regions")
            .long("bed-regions")
            .short('b')
            .num_args(1)
            .help("BED file with multiple regions for batch query (name start end per line)"),
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
            .short('m')
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
        Arg::new("min_degree")
            .long("min-degree")
            .num_args(1)
            .default_value("0")
            .value_parser(clap::value_parser!(usize))
            .help("Minimum distinct query sequences per region (default: 0 = off)"),
    )
    .arg(
        Arg::new("min_chain_length")
            .long("min-chain-length")
            .num_args(1)
            .default_value("0")
            .value_parser(clap::value_parser!(i32))
            .help("Minimum total aligned length per query (default: 0 = off)"),
    )
    .arg(
        Arg::new("subset_list")
            .long("subset-sequence-list")
            .num_args(1)
            .help("File with sequence names to include (one per line)"),
    )
    .arg(
        Arg::new("syntenic_filter")
            .long("syntenic-filter")
            .num_args(1)
            .help("UCSC chain file; drop query results whose query interval is not covered by any chain's query span (chain-level, both target and query name must match)"),
    )
}

/// Add POA scoring arguments to a clap Command.
/// Shared by `paf to-gfa`, `paf to-vcf`, `paf to-fas`, and `paf to-maf`.
pub fn add_poa_args(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("match_score")
            .long("match")
            .num_args(1)
            .default_value("5")
            .value_parser(clap::value_parser!(i32))
            .allow_negative_numbers(true)
            .help("POA match score (default: 5)"),
    )
    .arg(
        Arg::new("mismatch_score")
            .long("mismatch")
            .num_args(1)
            .default_value("-4")
            .value_parser(clap::value_parser!(i32))
            .allow_negative_numbers(true)
            .help("POA mismatch score (default: -4)"),
    )
    .arg(
        Arg::new("gap_open")
            .long("gap-open")
            .num_args(1)
            .default_value("-8")
            .value_parser(clap::value_parser!(i32))
            .allow_negative_numbers(true)
            .help("POA gap open penalty (default: -8)"),
    )
    .arg(
        Arg::new("gap_extend")
            .long("gap-extend")
            .num_args(1)
            .default_value("-6")
            .value_parser(clap::value_parser!(i32))
            .allow_negative_numbers(true)
            .help("POA gap extend penalty (default: -6)"),
    )
}

/// Shared query logic: parse args, build/load index, run queries, apply filters.
/// Returns the index and a list of (region, results) pairs.
#[allow(clippy::type_complexity)]
pub fn run_query(
    args: &ArgMatches,
) -> anyhow::Result<(PafIndex, Vec<((String, i32, i32), Vec<QueryResult>)>)> {
    let infile = args.get_one::<String>("infile").unwrap();
    let region_str = args.get_one::<String>("region");
    let bed_regions_path = args.get_one::<String>("bed_regions");
    let transitive = args.get_flag("transitive");
    let max_depth = *args.get_one::<u16>("max_depth").unwrap();
    let min_len = *args.get_one::<i32>("min_len").unwrap();
    let min_dist = *args.get_one::<i32>("min_dist").unwrap();
    let min_identity = *args.get_one::<f64>("min_identity").unwrap();
    let min_output_len = *args.get_one::<i32>("min_output_len").unwrap();
    let merge_distance = *args.get_one::<i32>("merge_distance").unwrap();
    let min_degree = *args.get_one::<usize>("min_degree").unwrap();
    let min_chain_length = *args.get_one::<i32>("min_chain_length").unwrap();
    let syntenic_filter_path = args.get_one::<String>("syntenic_filter");

    // Region input: exactly one of positional <region> or -b/--bed-regions
    anyhow::ensure!(
        region_str.is_some() || bed_regions_path.is_some(),
        "either positional <region> or -b/--bed-regions must be provided"
    );
    anyhow::ensure!(
        !(region_str.is_some() && bed_regions_path.is_some()),
        "<region> and -b/--bed-regions are mutually exclusive"
    );

    let regions: Vec<(String, i32, i32)> = if let Some(path) = bed_regions_path {
        load_bed_regions(path)?
    } else {
        let (name, start, end) = parse_region(region_str.unwrap())?;
        vec![(name.to_string(), start, end)]
    };

    let idx = if infile.ends_with(".paf.idx") {
        eprintln!("Loading index from {infile}...");
        PafIndex::load(infile)?
    } else {
        eprintln!("Building index from {infile}...");
        // Use build_from_path to enable lazy CIGAR loading for BGZF files.
        PafIndex::build_from_path(infile)?
    };

    eprintln!(
        "  sequences: {}, targets: {}",
        idx.names.len(),
        idx.num_targets()
    );
    if idx.is_lazy() {
        eprintln!("  mode: lazy (BGZF virtual-position CIGAR)");
    }

    let subset = args.get_one::<String>("subset_list").map(|list_path| {
        intspan::read_first_column(list_path)
            .into_iter()
            .collect::<HashSet<String>>()
    });

    // Optional syntenic filter: load UCSC chain file and build
    // (t_name, q_name) -> Vec<(q_start, q_end)> map for chain-level query coverage check.
    let syntenic_map: Option<HashMap<(String, String), Vec<(u64, u64)>>> =
        if let Some(path) = syntenic_filter_path {
            eprintln!("Loading syntenic chains from {path}...");
            let f = fs::File::open(path)?;
            let chains = read_chains(f)?;
            let mut map: HashMap<(String, String), Vec<(u64, u64)>> = HashMap::new();
            for c in &chains {
                let key = (c.header.t_name.clone(), c.header.q_name.clone());
                map.entry(key)
                    .or_default()
                    .push((c.header.q_start, c.header.q_end));
            }
            eprintln!(
                "  loaded {} chains ({} unique name pairs)",
                chains.len(),
                map.len()
            );
            Some(map)
        } else {
            None
        };

    let mut all_results: Vec<((String, i32, i32), Vec<QueryResult>)> = Vec::new();
    let mut total_results = 0usize;

    for (target_name, start, end) in &regions {
        let target_id = match idx.name_to_id(target_name) {
            Some(id) => id,
            None => {
                eprintln!("target '{target_name}' not found in index, skipping");
                continue;
            }
        };

        let mut results = if transitive {
            idx.query_transitive_bfs(
                target_id,
                *start,
                *end,
                max_depth,
                min_len,
                min_dist,
                min_identity,
                min_output_len,
                merge_distance,
            )
        } else {
            idx.query(target_id, *start, *end, min_identity, min_output_len)
        };

        if let Some(ref subset) = subset {
            results.retain(|(qid, _, _, _, _, _, _)| {
                let name = idx.id_to_name(*qid).unwrap_or("");
                subset.contains(name)
            });
        }

        if let Some(ref syntenic) = syntenic_map {
            let before = results.len();
            results.retain(|(qid, qiv, _, _, _, _, _)| {
                let q_name = idx.id_to_name(*qid).unwrap_or("");
                let key = (target_name.clone(), q_name.to_string());
                match syntenic.get(&key) {
                    None => false,
                    Some(spans) => {
                        let qs = qiv.first as u64;
                        let qe = qiv.last as u64;
                        spans.iter().any(|&(cs, ce)| qs < ce && qe > cs)
                    }
                }
            });
            let dropped = before - results.len();
            if dropped > 0 {
                eprintln!("  syntenic-filter: dropped {dropped} non-syntenic results for {target_name}:{start}-{end}");
            }
        }

        if min_chain_length > 0 {
            filter_by_chain_length(&mut results, min_chain_length);
        }

        if min_degree > 0 {
            let distinct: HashSet<u32> =
                results.iter().map(|(qid, _, _, _, _, _, _)| *qid).collect();
            if distinct.len() < min_degree {
                eprintln!(
                    "region {target_name}:{start}-{end} skipped (degree {} < min-degree {min_degree})",
                    distinct.len()
                );
                continue;
            }
        }

        total_results += results.len();
        all_results.push(((target_name.clone(), *start, *end), results));
    }

    if total_results == 0 {
        eprintln!("No results found.");
    } else {
        eprintln!("Total results: {total_results}");
    }

    Ok((idx, all_results))
}

// Parse a region string "name:start-end" (0-based, PAF convention).
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

// Parse BED file (name start end per line, tab-separated). Skips blanks and comments.
fn load_bed_regions(path: &str) -> anyhow::Result<Vec<(String, i32, i32)>> {
    let f = fs::File::open(path)?;
    let mut regions = Vec::new();
    for line in std::io::BufReader::new(f).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        anyhow::ensure!(
            fields.len() >= 3,
            "invalid BED line '{line}': expected at least 3 tab-separated fields"
        );
        let name = fields[0].to_string();
        let start: i32 = fields[1].parse()?;
        let end: i32 = fields[2].parse()?;
        regions.push((name, start, end));
    }
    Ok(regions)
}

// Drop queries whose total aligned length (summed across all result intervals
// for that query_id) is below `min_chain_length`. Operates in place.
fn filter_by_chain_length(results: &mut Vec<QueryResult>, min_chain_length: i32) {
    let mut totals: HashMap<u32, i32> = HashMap::new();
    for (qid, q_iv, _, _, _, _, _) in results.iter() {
        let len = (q_iv.last - q_iv.first).abs();
        *totals.entry(*qid).or_insert(0) += len;
    }
    results.retain(|(qid, _, _, _, _, _, _)| {
        totals.get(qid).copied().unwrap_or(0) >= min_chain_length
    });
}

// Output PAF: 12 columns + gi/bi/cg tags.
fn output_paf(idx: &PafIndex, results: &[QueryResult]) {
    for (query_id, q_iv, t_iv, cigar, _, _, strand) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");
        let block_len = (q_iv.last - q_iv.first).abs().max(1) as u32;
        let matches = pgr::libs::paf::cigar::cigar_stats(cigar).matches;
        let gi = pgr::libs::paf::cigar::gap_compressed_identity(cigar);
        let bi = pgr::libs::paf::cigar::block_identity(cigar);
        let cg = pgr::libs::paf::cigar::format_cigar(cigar);
        let (qs, qe) = if q_iv.first <= q_iv.last {
            (q_iv.first, q_iv.last)
        } else {
            (q_iv.last, q_iv.first)
        };
        let (ts, te) = if t_iv.first <= t_iv.last {
            (t_iv.first, t_iv.last)
        } else {
            (t_iv.last, t_iv.first)
        };
        println!(
            "{}\t0\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\t255\tgi:f:{:.6}\tbi:f:{:.6}\tcg:Z:{}",
            qname, qs, qe, strand, tname, ts, te, matches, block_len, gi, bi, cg
        );
    }
}

pub fn make_subcommand() -> Command {
    let cmd = Command::new("query")
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

Region input (one of):
* Positional <region>: single region (e.g. chr1:1000-5000)
* -b/--bed-regions <file>: BED file with multiple regions (one per line,
  tab-separated `name start end`), enabling batch query

Output: PAF (12 columns + gi/bi/cg tags). For BED or MAF output, use
`pgr paf to-bed` or `pgr paf to-maf` respectively.

Notes:
* Input PAF files should contain cg:Z: tags for accurate projection
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'

Examples:
1. Single-hop projection from a PAF file:
   pgr paf query alignments.paf chr1:1000-5000

2. Query from a saved index (faster):
   pgr paf query alignments.paf.idx chr1:1000-5000

3. Transitive BFS with filters:
   pgr paf query alignments.paf chr1:1000-5000 --transitive --min-identity 0.8

4. Batch query:
   pgr paf query alignments.paf.idx -b regions.bed

"###,
        );
    add_query_args(cmd)
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let (idx, all_results) = run_query(args)?;
    for (_, results) in &all_results {
        output_paf(&idx, results);
    }
    Ok(())
}
