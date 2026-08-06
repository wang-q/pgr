//! `pgr align fill` — 2D gap fill between pgi anchors with LASTZ.

use anyhow::Result;
use clap::{value_parser, Arg, ArgMatches, Command};
use cmd_lib::run_cmd;
use pgr::libs::fmt::psl::Psl;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use super::common;

/// Build the clap subcommand for fill.
pub fn make_subcommand() -> Command {
    Command::new("fill")
        .about("Fills 2D gaps between pgi anchors with LASTZ (syntenic)")
        .after_help(
            r###"
Combines the speed of `align pgi` with the sensitivity of LASTZ (the
FastGA-gapfill idea): pgi produces coarse anchors, LASTZ fills the colinear
gaps between consecutive same-strand anchors, and the two PSL sets are emitted
together.

Only syntenic (colinear, same-strand) searches are supported. Feed the combined
PSL to `pgr pl chainnet --syn` for the final chain/net/axt/maf.

* `align pgi` runs first (or reuse an existing PSL with --avail-psl, e.g.
  produced by pgi, FastGA or minimap2).
* For each pair of adjacent, non-overlapping, same-strand anchors whose target
  and query gaps both fall in [--min-gap, --max-gap], a bounding box is built
  overlapping the anchors by --overlap bp (LASTZ seeding buffer) and LASTZ is
  run on the extracted sub-sequences.
* LASTZ output (LAV) is converted to PSL and the subrange coordinates are
  lifted back to genomic coordinates (the `pgr psl lift` logic).
* No dedup is done here: the anchors and the LASTZ records are written together
  and `pgr pl chainnet` handles the overlap/merge.

Notes:
* `lastz` must be installed and available in PATH.
* Converts the inputs to 2bit in a tempdir and extracts each box with
  `pgr 2bit range` (random access — no whole genome in memory).
* --overlap defaults to 1 kb (paper's box overlap); --min-gap defaults to
  100 bp (smaller gaps are left to pgi); --max-gap optionally skips longer
  gaps (novel-sequence regions); by default no gap is skipped.

Examples:
1. Default gap fill:
   pgr align fill ref.fa query.fa -o out.psl
2. Reuse an existing PSL (pgi, FastGA, minimap2...):
   pgr align fill ref.fa query.fa --avail-psl anchors.psl -o out.psl
3. Larger seeding buffer and a close-species preset:
   pgr align fill ref.fa query.fa --preset set01 --overlap 2000 -o out.psl
"###,
        )
        .arg(
            Arg::new("target")
                .index(1)
                .required(true)
                .help("Target genome (FASTA or .2bit)"),
        )
        .arg(
            Arg::new("query")
                .index(2)
                .required(true)
                .help("Query genome (FASTA or .2bit)"),
        )
        .arg(
            Arg::new("avail_psl")
                .long("avail-psl")
                .help("Precomputed PSL anchor file (skips the internal align pgi)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("preset")
                .long("preset")
                .value_parser(clap::builder::PossibleValuesParser::new(
                    pgr::libs::lastz::preset_names(),
                ))
                .help("Use a predefined LASTZ parameter set (set01..set07)"),
        )
        .arg(
            Arg::new("overlap")
                .long("overlap")
                .default_value("1000")
                .value_parser(value_parser!(i32))
                .help("Box expansion beyond the gap (bp), LASTZ seeding buffer"),
        )
        .arg(
            Arg::new("min_gap")
                .long("min-gap")
                .default_value("100")
                .value_parser(value_parser!(i32))
                .help("Shortest gap to fill (bp); smaller gaps are left to pgi"),
        )
        .arg(
            Arg::new("max_gap")
                .long("max-gap")
                .value_parser(value_parser!(i32))
                .help("Skip gaps longer than this (bp); default: no limit"),
        )
        .arg(
            Arg::new("query_depth")
                .long("query-depth")
                .default_value("50")
                .value_parser(value_parser!(usize))
                .help("Query depth threshold for LASTZ"),
        )
        .arg(
            Arg::new("lastz_args")
                .long("lastz-args")
                .help("Additional arguments passed directly to LASTZ (overrides preset)"),
        )
        .arg(crate::cmd_pgr::args::parallel_arg_with_default("8"))
}

/// Execute the fill command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    if which::which("lastz").is_err() {
        anyhow::bail!("lastz not found in PATH. Please install lastz first.");
    }
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut protected = vec![
        args.get_one::<String>("target").unwrap().as_str(),
        args.get_one::<String>("query").unwrap().as_str(),
    ];
    if let Some(avail_psl) = args.get_one::<String>("avail_psl") {
        protected.push(avail_psl.as_str());
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, protected.iter().copied())?;

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_align_fill_")?;
    let pgr = ctx.pgr.clone();

    run_cmd!(info "==> Absolute paths")?;
    let abs_target = ctx.abs_path(args.get_one::<String>("target").unwrap())?;
    let abs_query = ctx.abs_path(args.get_one::<String>("query").unwrap())?;
    let abs_outfile = pgr::libs::pl::abs_path_or_stdout(outfile)?;
    let abs_avail_psl = match args.get_one::<String>("avail_psl") {
        Some(p) => Some(ctx.abs_path(p)?),
        None => None,
    };

    let _cwd_guard = ctx.enter()?;

    // 1. Anchors: reuse an existing PSL or run `align pgi`.
    let avail_psl_path = match &abs_avail_psl {
        Some(p) => p.clone(),
        None => {
            run_cmd!(info "==> align pgi")?;
            run_cmd!(${pgr} align pgi ${abs_target} ${abs_query} -o pgi.psl)?;
            "pgi.psl".to_string()
        }
    };

    let avail_psl = common::read_psl(&avail_psl_path)?;
    log::info!("available anchors: {}", avail_psl.len());

    // 2. Convert to 2bit for random-access box extraction.
    run_cmd!(info "==> fa to-2bit")?;
    let abs_target_2bit = common::to_2bit(&pgr, &abs_target, "target.2bit")?;
    let abs_query_2bit = common::to_2bit(&pgr, &abs_query, "query.2bit")?;

    // 3. Gap fill: emit the anchors plus the LASTZ records.
    let opts = FillOptions {
        overlap: *args.get_one::<i32>("overlap").unwrap(),
        min_gap: *args.get_one::<i32>("min_gap").unwrap(),
        max_gap: args.get_one::<i32>("max_gap").copied(),
        preset: args.get_one::<String>("preset").cloned(),
        query_depth: *args.get_one::<usize>("query_depth").unwrap(),
        lastz_args: args.get_one::<String>("lastz_args").cloned(),
        parallel: *args.get_one::<usize>("parallel").unwrap(),
    };
    let combined = run_fill(
        &pgr,
        &abs_target_2bit,
        &abs_query_2bit,
        &avail_psl,
        &opts,
        &std::env::current_dir()?,
    )?;

    // 4. Write the combined PSL (anchors + LASTZ gap-fill).
    run_cmd!(info "==> Write combined PSL")?;
    let mut writer = pgr::libs::io::writer(&abs_outfile)?;
    for p in &combined {
        p.write_to(&mut writer)?;
    }
    log::info!(
        "wrote {} PSL records ({} anchors + {} lastz) to {}",
        combined.len(),
        avail_psl.len(),
        combined.len().saturating_sub(avail_psl.len()),
        abs_outfile
    );
    Ok(())
}

// --- Gap-fill logic ---

/// A 2D bounding box between two colinear anchors, in 0-based half-open
/// forward-strand coordinates on both genomes.
#[derive(Debug, Clone)]
struct Box {
    t_name: String,
    t_start: i32,
    t_end: i32,
    q_name: String,
    q_start: i32,
    q_end: i32,
}

/// Options controlling the gap-fill run.
#[derive(Debug, Clone)]
struct FillOptions {
    /// Box overlaps the anchor by this many bp on each side (LASTZ seeding buffer).
    overlap: i32,
    /// Gaps shorter than this (bp) are left to pgi.
    min_gap: i32,
    /// Gaps longer than this (bp) are skipped (likely novel sequence).
    max_gap: Option<i32>,
    /// LASTZ preset (set01..set07); None uses lastz defaults.
    preset: Option<String>,
    /// Query depth threshold passed to lastz.
    query_depth: usize,
    /// Extra LASTZ arguments appended after the preset.
    lastz_args: Option<String>,
    /// Number of parallel lastz jobs.
    parallel: usize,
}

/// Compute the fill boxes between consecutive colinear same-strand anchors.
///
/// Anchors are the PSL records (any aligner: pgi, FastGA, minimap2...),
/// grouped by (target, query, strand) and sorted by target start. A box is
/// emitted for each pair of adjacent, non-overlapping, colinear anchors whose
/// target *and* query gaps both fall in `[min_gap, max_gap]`; the box is
/// expanded beyond the gap by `overlap` bp on each side so LASTZ can seed
/// across the true homology boundary.
fn compute_boxes(psls: &[Psl], overlap: i32, min_gap: i32, max_gap: i32) -> Vec<Box> {
    let mut groups: BTreeMap<(String, String, String), Vec<&Psl>> = BTreeMap::new();
    for p in psls {
        groups
            .entry((p.t_name.clone(), p.q_name.clone(), p.strand.clone()))
            .or_default()
            .push(p);
    }

    let mut boxes = Vec::new();
    for ((t_name, q_name, strand), mut anchors) in groups {
        anchors.sort_by_key(|p| (p.t_start, p.t_end));
        for pair in anchors.windows(2) {
            let a = pair[0];
            let b = pair[1];
            // Non-overlapping, ascending target order.
            if a.t_end >= b.t_start {
                continue;
            }
            let t_gap_start = a.t_end;
            let t_gap_end = b.t_start;
            // Query gap depends on strand: '+' grows with t, '-' shrinks.
            let (q_gap_start, q_gap_end) = match strand.as_str() {
                "+" => {
                    if a.q_end >= b.q_start {
                        continue;
                    }
                    (a.q_end, b.q_start)
                }
                "-" => {
                    if b.q_end >= a.q_start {
                        continue;
                    }
                    (b.q_end, a.q_start)
                }
                _ => continue,
            };
            let t_gap_len = t_gap_end - t_gap_start;
            let q_gap_len = q_gap_end - q_gap_start;
            if t_gap_len < min_gap || t_gap_len > max_gap {
                continue;
            }
            if q_gap_len < min_gap || q_gap_len > max_gap {
                continue;
            }
            boxes.push(Box {
                t_name: t_name.clone(),
                t_start: (t_gap_start - overlap).max(0),
                t_end: t_gap_end + overlap,
                q_name: q_name.clone(),
                q_start: (q_gap_start - overlap).max(0),
                q_end: q_gap_end + overlap,
            });
        }
    }
    boxes
}

/// Run the gap-fill pipeline over two genomes.
fn run_fill(
    pgr: &str,
    target_2bit: &str,
    query_2bit: &str,
    avail_psl: &[Psl],
    opts: &FillOptions,
    tempdir: &Path,
) -> Result<Vec<Psl>> {
    let mut sizes = common::read_2bit_sizes(pgr, target_2bit)?;
    sizes.extend(common::read_2bit_sizes(pgr, query_2bit)?);

    let mut boxes = compute_boxes(
        avail_psl,
        opts.overlap,
        opts.min_gap,
        opts.max_gap.unwrap_or(i32::MAX),
    );
    // Clamp box coordinates to contig bounds so the extracted subrange is
    // valid; otherwise the lift would reject the record (subrange end > real
    // size) and drop the alignment.
    for b in &mut boxes {
        let t_size = sizes.get(&b.t_name).copied().unwrap_or(0);
        let q_size = sizes.get(&b.q_name).copied().unwrap_or(0);
        b.t_start = b.t_start.clamp(0, t_size);
        b.t_end = b.t_end.clamp(b.t_start, t_size);
        b.q_start = b.q_start.clamp(0, q_size);
        b.q_end = b.q_end.clamp(b.q_start, q_size);
    }
    boxes.retain(|b| b.t_end > b.t_start && b.q_end > b.q_start);
    if boxes.is_empty() {
        return Ok(avail_psl.to_vec());
    }

    let (common_args, _matrix_handle) = common::build_common_args(
        opts.preset.as_deref(),
        opts.query_depth,
        opts.lastz_args.as_deref(),
    )?;

    let t_extract = std::time::Instant::now();
    std::fs::create_dir_all(tempdir)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.parallel)
        .build()?;
    let jobs: Vec<Result<(PathBuf, PathBuf, PathBuf)>> = pool.install(|| {
        use rayon::prelude::*;
        boxes
            .par_iter()
            .enumerate()
            .map(|(i, b)| {
                let t_fa = tempdir.join(format!("t_{i}.fa"));
                let q_fa = tempdir.join(format!("q_{i}.fa"));
                let lav = tempdir.join(format!("box_{i}.lav"));
                common::extract_2bit_range(pgr, target_2bit, &b.t_name, b.t_start, b.t_end, &t_fa)?;
                common::extract_2bit_range(pgr, query_2bit, &b.q_name, b.q_start, b.q_end, &q_fa)?;
                Ok((t_fa, q_fa, lav))
            })
            .collect()
    });
    let jobs: Vec<(PathBuf, PathBuf, PathBuf)> = jobs.into_iter().collect::<Result<_>>()?;
    log::info!(
        "extracted {} boxes in {:?}",
        boxes.len(),
        t_extract.elapsed()
    );

    let t_lastz = std::time::Instant::now();
    let lastz_psl = common::run_lastz_jobs(&jobs, &common_args, &sizes, opts.parallel)?;
    log::info!("lastz + convert in {:?}", t_lastz.elapsed());

    let mut out = avail_psl.to_vec();
    out.extend(lastz_psl);
    Ok(out)
}

// --- Tests ---

/// A minimal PSL helper for tests (target/query starts/ends and strand only).
#[cfg(test)]
fn tpsl(
    t_name: &str,
    t_start: i32,
    t_end: i32,
    q_name: &str,
    q_start: i32,
    q_end: i32,
    strand: &str,
) -> Psl {
    Psl {
        t_name: t_name.into(),
        t_start,
        t_end,
        q_name: q_name.into(),
        q_start,
        q_end,
        strand: strand.into(),
        ..Default::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compute_boxes_plus_strand_emits_gap_boxes() {
        let anchors = vec![
            tpsl("t", 0, 100, "q", 0, 100, "+"),
            tpsl("t", 500, 600, "q", 500, 600, "+"),
            tpsl("t", 1000, 1100, "q", 1000, 1100, "+"),
        ];
        let boxes = compute_boxes(&anchors, 10, 50, 1000);
        // Two pairs: [0,100]-[500,600] and [500,600]-[1000,1100].
        let coords: Vec<(i32, i32, i32, i32)> = boxes
            .iter()
            .map(|b| (b.t_start, b.t_end, b.q_start, b.q_end))
            .collect();
        assert_eq!(coords, vec![(90, 510, 90, 510), (590, 1010, 590, 1010)]);
    }

    #[test]
    fn compute_boxes_skips_overlapping_and_out_of_range_gaps() {
        // Overlapping anchors: no box.
        let anchors = vec![
            tpsl("t", 0, 500, "q", 0, 500, "+"),
            tpsl("t", 400, 900, "q", 400, 900, "+"),
        ];
        assert!(compute_boxes(&anchors, 10, 50, 1000).is_empty());
        // Gap below min_gap: no box.
        let anchors = vec![
            tpsl("t", 0, 100, "q", 0, 100, "+"),
            tpsl("t", 120, 220, "q", 120, 220, "+"),
        ];
        assert!(compute_boxes(&anchors, 10, 50, 1000).is_empty());
        // Gap above max_gap: no box.
        let anchors = vec![
            tpsl("t", 0, 100, "q", 0, 100, "+"),
            tpsl("t", 5000, 5100, "q", 5000, 5100, "+"),
        ];
        assert!(compute_boxes(&anchors, 10, 50, 1000).is_empty());
    }

    #[test]
    fn compute_boxes_minus_strand_query_gap_is_reversed() {
        // On the '-' strand the query order is reversed: the earlier target
        // anchor has the later query span.
        let a = tpsl("t", 0, 100, "q", 900, 1000, "-");
        let b = tpsl("t", 500, 600, "q", 300, 400, "-");
        let boxes = compute_boxes(&[a, b], 10, 50, 1000);
        assert_eq!(boxes.len(), 1);
        // Query gap is [400, 900] (forward strand), box q -> [390, 910].
        assert_eq!(boxes[0].q_start, 390);
        assert_eq!(boxes[0].q_end, 910);
    }
}
