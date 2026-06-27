use clap::*;
use pgr::libs::paf::cigar::CigarOp;
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

Region input (one of):
* Positional <region>: single region (e.g. chr1:1000-5000)
* -b/--bed-regions <file>: BED file with multiple regions (one per line,
  tab-separated `name start end`), enabling batch query

Output formats (-o):
* paf (default): PAF 12 columns + tags (gi/bi/cg)
* bed: BED3 (name start end), most pipe-friendly

Notes:
* Input PAF files should contain cg:Z: tags for accurate projection
* Reads from stdin if input file is 'stdin'

Examples:
1. Single-hop projection from a PAF file:
   pgr paf query alignments.paf chr1:1000-5000

2. Query from a saved index (faster):
   pgr paf query alignments.paf.idx chr1:1000-5000

3. Transitive BFS with filters:
   pgr paf query alignments.paf chr1:1000-5000 --transitive --min-identity 0.8

4. Batch query with BED output:
   pgr paf query alignments.paf.idx -b regions.bed -o bed

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
            Arg::new("output_format")
                .long("output")
                .short('o')
                .num_args(1)
                .default_value("paf")
                .value_parser(["paf", "bed"])
                .help("Output format: paf (default) or bed"),
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

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let region_str = args.get_one::<String>("region");
    let bed_regions_path = args.get_one::<String>("bed_regions");
    let output_format = args.get_one::<String>("output_format").unwrap();
    let transitive = args.get_flag("transitive");
    let max_depth = *args.get_one::<u16>("max_depth").unwrap();
    let min_len = *args.get_one::<i32>("min_len").unwrap();
    let min_dist = *args.get_one::<i32>("min_dist").unwrap();
    let min_identity = *args.get_one::<f64>("min_identity").unwrap();
    let min_output_len = *args.get_one::<i32>("min_output_len").unwrap();
    let merge_distance = *args.get_one::<i32>("merge_distance").unwrap();

    // Region input: exactly one of positional <region> or -b/--bed-regions
    anyhow::ensure!(
        region_str.is_some() || bed_regions_path.is_some(),
        "either positional <region> or -b/--bed-regions must be provided"
    );
    anyhow::ensure!(
        !(region_str.is_some() && bed_regions_path.is_some()),
        "<region> and -b/--bed-regions are mutually exclusive"
    );

    // Collect regions to query (single or batch)
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
        let reader = pgr::reader(infile);
        PafIndex::build(reader)?
    };

    eprintln!(
        "  sequences: {}, targets: {}",
        idx.names.len(),
        idx.num_targets()
    );

    // Subset filter
    let subset = if let Some(list_path) = args.get_one::<String>("subset_list") {
        Some(load_subset(list_path)?)
    } else {
        None
    };

    let mut total_results = 0usize;
    let use_bed = output_format == "bed";

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
            results.retain(|(qid, _, _, _)| {
                let name = idx.id_to_name(*qid).unwrap_or("");
                subset.contains(name)
            });
        }

        if use_bed {
            output_bed(&idx, &results);
        } else {
            output_paf(&idx, &results);
        }
        total_results += results.len();
    }

    if total_results == 0 {
        eprintln!("No results found.");
    } else {
        eprintln!("Total results: {total_results}");
    }

    Ok(())
}

fn output_bed(
    idx: &PafIndex,
    results: &[(
        u32,
        coitrees::Interval<u32>,
        coitrees::Interval<u32>,
        Vec<CigarOp>,
    )],
) {
    for (query_id, q_iv, _t_iv, _cigar) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let (qs, qe) = if q_iv.first <= q_iv.last {
            (q_iv.first, q_iv.last)
        } else {
            (q_iv.last, q_iv.first)
        };
        println!("{qname}\t{qs}\t{qe}");
    }
}

fn output_paf(
    idx: &PafIndex,
    results: &[(
        u32,
        coitrees::Interval<u32>,
        coitrees::Interval<u32>,
        Vec<CigarOp>,
    )],
) {
    for (query_id, q_iv, t_iv, cigar) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");
        let block_len = (q_iv.last - q_iv.first).abs().max(1) as u32;
        let matches = pgr::libs::paf::cigar::cigar_stats(cigar).matches;
        let gi = pgr::libs::paf::cigar::gap_compressed_identity(cigar);
        let bi = pgr::libs::paf::cigar::block_identity(cigar);
        let cg = pgr::libs::paf::cigar::format_cigar(cigar);
        let strand = if q_iv.first <= q_iv.last { '+' } else { '-' };
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
