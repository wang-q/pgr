use clap::{Arg, ArgMatches, Command};

use pgr::libs::chain::record::read_chains;
use pgr::libs::paf::fasta::{load_fasta_tsv, validate_tsv_covers_index, FastaStore};
use pgr::libs::paf::index::{PafIndex, QueryResult};
use pgr::libs::paf::msa::orient_interval;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

// Re-export shared POA argument builder and extractor from `cmd_pgr::args`.
pub use crate::cmd_pgr::args::{add_poa_args, get_poa_params};

/// A region paired with its query results.
pub type QueryGroup = ((String, i32, i32), Vec<QueryResult>);

/// Add common query arguments to a clap Command.
/// Shared by `paf query`, `paf to-bed`, `paf to-gfa`, `paf to-vcf`,
/// `paf to-fas`, and `paf to-maf`.
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

/// Add the required `-f/--fasta-tsv` argument.
/// Shared by `paf to-gfa`, `paf to-vcf`, `paf to-fas`, and `paf to-maf`.
pub fn add_fasta_tsv_arg(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("fasta_tsv")
            .long("fasta-tsv")
            .short('f')
            .required(true)
            .num_args(1)
            .help("TSV file: genome_name <tab> bgzf_fasta_path"),
    )
}

/// Add the `--msa` flag for POA-based multi-way output.
/// Shared by `paf to-fas` and `paf to-maf`.
pub fn add_msa_flag(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("msa")
            .long("msa")
            .num_args(0)
            .help("Merge results per region into a multi-way block via POA"),
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
        log::info!("Loading index from {infile}...");
        PafIndex::load(infile)?
    } else {
        log::info!("Building index from {infile}...");
        // Use build_from_path to enable lazy CIGAR loading for BGZF files.
        PafIndex::build_from_path(infile)?
    };

    log::info!(
        "  sequences: {}, targets: {}",
        idx.names.len(),
        idx.num_targets()
    );
    if idx.is_lazy() {
        log::info!("  mode: lazy (BGZF virtual-position CIGAR)");
    }

    let subset = args
        .get_one::<String>("subset_list")
        .map(|list_path| pgr::libs::io::read_names_as_set(list_path))
        .transpose()?;

    // Optional syntenic filter: load UCSC chain file and build
    // (t_name, q_name) -> Vec<(q_start, q_end)> map for chain-level query coverage check.
    let syntenic_map: Option<HashMap<(String, String), Vec<(u64, u64)>>> =
        if let Some(path) = syntenic_filter_path {
            log::info!("Loading syntenic chains from {path}...");
            let chains = read_chains(pgr::reader(path)?)?;
            let mut map: HashMap<(String, String), Vec<(u64, u64)>> = HashMap::new();
            for c in &chains {
                let key = (c.header.t_name.clone(), c.header.q_name.clone());
                map.entry(key)
                    .or_default()
                    .push((c.header.q_start, c.header.q_end));
            }
            log::info!(
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
                log::warn!("target '{target_name}' not found in index, skipping");
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
                log::info!("  syntenic-filter: dropped {dropped} non-syntenic results for {target_name}:{start}-{end}");
            }
        }

        if min_chain_length > 0 {
            filter_by_chain_length(&mut results, min_chain_length);
        }

        if min_degree > 0 {
            let distinct: HashSet<u32> =
                results.iter().map(|(qid, _, _, _, _, _, _)| *qid).collect();
            if distinct.len() < min_degree {
                log::info!(
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
        log::info!("No results found.");
    } else {
        log::info!("Total results: {total_results}");
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

// Parse BED file (name start end per line, tab-separated). Skips blank and comments.
fn load_bed_regions(path: &str) -> anyhow::Result<Vec<(String, i32, i32)>> {
    let reader = pgr::reader(path)?;
    let mut regions = Vec::new();
    for line in reader.lines() {
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
pub fn output_paf(idx: &PafIndex, results: &[QueryResult]) {
    for (query_id, q_iv, t_iv, cigar, _, _, strand) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");
        let block_len = (q_iv.last - q_iv.first).abs().max(1) as u32;
        let matches = pgr::libs::paf::cigar::cigar_stats(cigar).matches;
        let gi = pgr::libs::paf::cigar::gap_compressed_identity(cigar);
        let bi = pgr::libs::paf::cigar::block_identity(cigar);
        let cg = pgr::libs::paf::cigar::format_cigar(cigar);
        let (qs, qe) = orient_interval(q_iv.first, q_iv.last);
        let (ts, te) = orient_interval(t_iv.first, t_iv.last);
        println!(
            "{}\t0\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\t255\tgi:f:{:.6}\tbi:f:{:.6}\tcg:Z:{}",
            qname, qs, qe, strand, tname, ts, te, matches, block_len, gi, bi, cg
        );
    }
}

/// Load fasta TSV, run query, validate TSV covers the index, build FastaStore.
/// Shared by `to-fas`, `to-gfa`, `to-maf`, `to-vcf`.
#[allow(clippy::type_complexity)]
pub fn prepare_query(args: &ArgMatches) -> anyhow::Result<(PafIndex, Vec<QueryGroup>, FastaStore)> {
    let tsv_path = args.get_one::<String>("fasta_tsv").unwrap();
    let seq_to_file = load_fasta_tsv(tsv_path)?;

    let (idx, all_results) = run_query(args)?;

    validate_tsv_covers_index(&seq_to_file, &idx)?;

    let fasta_store = FastaStore::new(&seq_to_file)?;

    Ok((idx, all_results, fasta_store))
}
