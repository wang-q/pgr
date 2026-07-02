//! Query orchestration for the PAF interval-tree index.
//!
//! Decouples query and filter logic from clap argument parsing: [`QueryOptions`]
//! holds raw parameters (parsed by the cmd layer), and [`run_query`] performs
//! index build/load, transitive BFS traversal, and subset / syntenic /
//! chain-length / degree filtering.

use super::cigar;
use super::index::{PafIndex, QueryResult};
use super::msa_build::orient_interval;
use crate::libs::chain::record::read_chains;
use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead, Write};

/// Parameters for a PAF index query, decoupled from clap.
///
/// Either [`region`](Self::region) or [`bed_regions`](Self::bed_regions) must
/// be set (mutually exclusive).
#[derive(Debug, Clone)]
pub struct QueryOptions {
    /// Input PAF file or `.paf.idx` index path.
    pub infile: String,
    /// Optional `name:start-end` region string.
    pub region: Option<String>,
    /// Optional BED file with multiple regions.
    pub bed_regions: Option<String>,
    /// Enable transitive BFS traversal.
    pub transitive: bool,
    /// Maximum BFS depth (0 = unlimited).
    pub max_depth: u16,
    /// Minimum interval length to propagate.
    pub min_len: i32,
    /// Minimum distance to merge adjacent intervals.
    pub min_dist: i32,
    /// Minimum gap-compressed identity (0.0-1.0).
    pub min_identity: f64,
    /// Minimum output interval length.
    pub min_output_len: i32,
    /// Merge adjacent output intervals within this distance.
    pub merge_distance: i32,
    /// Minimum distinct query sequences per region.
    pub min_degree: usize,
    /// Minimum total aligned length per query.
    pub min_chain_length: i32,
    /// Optional file with sequence names to include.
    pub subset_list: Option<String>,
    /// Optional UCSC chain file for syntenic filtering.
    pub syntenic_filter: Option<String>,
}

/// Run a query against a PAF index.
///
/// Builds or loads the index from `opts.infile`, parses regions from
/// `opts.region` or `opts.bed_regions`, runs BFS / direct query, and
/// applies subset / syntenic / chain-length / degree filters.
#[allow(clippy::type_complexity)]
pub fn run_query(
    opts: &QueryOptions,
) -> anyhow::Result<(PafIndex, Vec<((String, i32, i32), Vec<QueryResult>)>)> {
    anyhow::ensure!(
        opts.region.is_some() || opts.bed_regions.is_some(),
        "either <region> or --bed-regions must be provided"
    );
    anyhow::ensure!(
        !(opts.region.is_some() && opts.bed_regions.is_some()),
        "<region> and --bed-regions are mutually exclusive"
    );

    let regions: Vec<(String, i32, i32)> = if let Some(path) = &opts.bed_regions {
        load_bed_regions(path)?
    } else {
        let (name, start, end) = parse_region(opts.region.as_ref().unwrap())?;
        vec![(name.to_string(), start, end)]
    };

    let idx = if opts.infile.ends_with(".paf.idx") {
        log::info!("Loading index from {}...", opts.infile);
        PafIndex::load(&opts.infile)?
    } else {
        log::info!("Building index from {}...", opts.infile);
        PafIndex::build_from_path(&opts.infile)?
    };

    log::info!(
        "  sequences: {}, targets: {}",
        idx.names.len(),
        idx.num_targets()
    );
    if idx.is_lazy() {
        log::info!("  mode: lazy (BGZF virtual-position CIGAR)");
    }

    let subset = opts
        .subset_list
        .as_ref()
        .map(|list_path| {
            crate::libs::io::read_names::<std::collections::HashSet<String>>(list_path)
        })
        .transpose()?;

    // Optional syntenic filter: load UCSC chain file and build
    // (t_name, q_name) -> Vec<(q_start, q_end)> map for chain-level query coverage check.
    let syntenic_map: Option<HashMap<(String, String), Vec<(u64, u64)>>> =
        if let Some(path) = &opts.syntenic_filter {
            log::info!("Loading syntenic chains from {path}...");
            let chains = read_chains(crate::reader(path)?)?;
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

        let mut results = if opts.transitive {
            idx.query_transitive_bfs(
                target_id,
                *start,
                *end,
                opts.max_depth,
                opts.min_len,
                opts.min_dist,
                opts.min_identity,
                opts.min_output_len,
                opts.merge_distance,
            )
        } else {
            idx.query(
                target_id,
                *start,
                *end,
                opts.min_identity,
                opts.min_output_len,
            )
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
                log::info!(
                    "  syntenic-filter: dropped {dropped} non-syntenic results for {target_name}:{start}-{end}"
                );
            }
        }

        if opts.min_chain_length > 0 {
            filter_by_chain_length(&mut results, opts.min_chain_length);
        }

        if opts.min_degree > 0 {
            let distinct: HashSet<u32> =
                results.iter().map(|(qid, _, _, _, _, _, _)| *qid).collect();
            if distinct.len() < opts.min_degree {
                log::info!(
                    "region {target_name}:{start}-{end} skipped (degree {} < min-degree {})",
                    distinct.len(),
                    opts.min_degree
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

/// Parse a region string `name:start-end` (0-based, PAF convention).
pub fn parse_region(s: &str) -> anyhow::Result<(&str, i32, i32)> {
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

/// Parse a BED file (name start end per line, tab-separated). Skips blank and comments.
pub fn load_bed_regions(path: &str) -> anyhow::Result<Vec<(String, i32, i32)>> {
    let reader = crate::reader(path)?;
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

/// Drop queries whose total aligned length (summed across all result intervals
/// for that query_id) is below `min_chain_length`. Operates in place.
pub fn filter_by_chain_length(results: &mut Vec<QueryResult>, min_chain_length: i32) {
    let mut totals: HashMap<u32, i32> = HashMap::new();
    for (qid, q_iv, _, _, _, _, _) in results.iter() {
        let len = (q_iv.last - q_iv.first).abs();
        *totals.entry(*qid).or_insert(0) += len;
    }
    results.retain(|(qid, _, _, _, _, _, _)| {
        totals.get(qid).copied().unwrap_or(0) >= min_chain_length
    });
}

/// A region paired with its query results.
pub type QueryGroup = ((String, i32, i32), Vec<QueryResult>);

/// Write `results` as PAF records (12 columns + `gi`/`bi`/`cg` tags) to `writer`.
///
/// Sequence names are resolved via `idx.id_to_name`. The `cg` tag is built from
/// the per-result CIGAR; `gi` and `bi` are gap-compressed and block identities.
pub fn output_paf<W: Write>(
    writer: &mut W,
    idx: &PafIndex,
    results: &[QueryResult],
) -> io::Result<()> {
    for (query_id, q_iv, t_iv, cigar, _, _, strand) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");
        let block_len = (q_iv.last - q_iv.first).abs().max(1) as u32;
        let matches = cigar::cigar_stats(cigar).matches;
        let gi = cigar::gap_compressed_identity(cigar);
        let bi = cigar::block_identity(cigar);
        let cg = cigar::format_cigar(cigar);
        let (qs, qe) = orient_interval(q_iv.first, q_iv.last);
        let (ts, te) = orient_interval(t_iv.first, t_iv.last);
        writeln!(
            writer,
            "{}\t0\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\t255\tgi:f:{:.6}\tbi:f:{:.6}\tcg:Z:{}",
            qname, qs, qe, strand, tname, ts, te, matches, block_len, gi, bi, cg
        )?;
    }
    Ok(())
}
