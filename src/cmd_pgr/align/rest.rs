//! `pgr align rest` — fill the whole-genome complement of pgi anchors with
//! LASTZ, with the query side trimmed the same way (no 2D coordinates).

use anyhow::Result;
use clap::{value_parser, Arg, ArgMatches, Command};
use cmd_lib::run_cmd;
use pgr::libs::ds::IntSpan;
use pgr::libs::fmt::psl::Psl;
use pgr::libs::runlist::{span_op, SpanOp};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use super::common;

/// Build the clap subcommand for rest.
pub fn make_subcommand() -> Command {
    Command::new("rest")
        .about("Searches beyond the given PSL anchors with LASTZ, filling the rest of the genome (syntenic)")
        .after_help(
            r###"
The name `rest` means "the rest of the genome". Given a PSL of anchors (from
`pgr align pgi` or `--avail-psl`), this command looks **beyond those
alignments**: it computes the target-side regions the anchors do not cover
(trim -> excise small anchors -> whole-genome holes) and tries again to find
homology for them with LASTZ. The query side gets the same treatment, and a
k-mer prefilter pairs target/query holes so only likely pairs are aligned.
`pgr align fill` only fills the gaps between anchors; `rest` fills everything
else.

Only syntenic (colinear, same-strand) searches are supported. Feed the combined
PSL to `pgr pl chainnet --syn` for the final chain/net/axt/maf.

* `align pgi` runs first (or reuse an existing PSL with --avail-psl, e.g.
  produced by pgi, FastGA or minimap2).
* Each side is handled independently with 1D runlist operations (no 2D
  coordinates): PSL spans -> trim -> excise -> whole-genome holes.
* The query holes are extracted as one multi-sequence FASTA and every target
  hole is aligned against it (LASTZ supports multiple query sequences).
* LASTZ output (LAV) is converted to PSL and the subrange coordinates are
  lifted back to genomic coordinates (the `pgr psl lift` logic).
* No dedup is done here: the anchors and the LASTZ records are written together
  and `pgr pl chainnet` handles the overlap/merge.

Notes:
* `lastz` must be installed and available in PATH.
* Converts the inputs to 2bit in a tempdir and extracts each hole with
  `pgr 2bit range` (random access — no whole genome in memory).
* --trim shrinks every anchor span by the given bp before excising;
  --min-anchor drops anchors shorter than the given bp (both before the
  complement, like the `rept s-kmer` pipeline); --max-gap optionally skips
  holes longer than the given bp; by default no hole is skipped.

Examples:
1. Default complement fill:
   pgr align rest ref.fa query.fa -o out.psl
2. Reuse an existing PSL:
   pgr align rest ref.fa query.fa --avail-psl anchors.psl -o out.psl
3. Drop tiny anchors and skip novel-sequence holes:
   pgr align rest ref.fa query.fa --min-anchor 1000 --max-gap 100000 -o out.psl
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
            Arg::new("trim")
                .long("trim")
                .default_value("500")
                .value_parser(value_parser!(i32))
                .help("Shrink each anchor span by this many bp on both ends"),
        )
        .arg(
            Arg::new("min_anchor")
                .long("min-anchor")
                .default_value("500")
                .value_parser(value_parser!(i32))
                .help("Excise anchors shorter than this (bp) before the complement"),
        )
        .arg(
            Arg::new("max_gap")
                .long("max-gap")
                .value_parser(value_parser!(i32))
                .help("Skip holes longer than this (bp); default: no limit"),
        )
        .arg(
            Arg::new("sampler")
                .long("sampler")
                .default_value("syncmer")
                .value_parser(["syncmer", "minimizer", "none"])
                .help("Fragment prefilter sampler (syncmer, minimizer or none = every target hole against the full query-holes set)"),
        )
        .arg(
            Arg::new("smer")
                .long("smer")
                .default_value("17")
                .value_parser(value_parser!(usize))
                .help("Syncmer s-mer length for the prefilter (17 = fast; 15 = higher coverage)"),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .default_value("17")
                .value_parser(value_parser!(usize))
                .help("Minimizer k-mer length for the prefilter"),
        )
        .arg(
            Arg::new("window")
                .long("window")
                .default_value("5")
                .value_parser(value_parser!(usize))
                .help("Prefilter sampler window size"),
        )
        .arg(
            Arg::new("min_shared")
                .long("min-shared")
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Minimum shared sampled k-mers to pair a target/query hole"),
        )
        .arg(
            Arg::new("unmatched")
                .long("unmatched")
                .default_value("skip")
                .value_parser(["skip", "full"])
                .help("How to treat target holes with no paired query hole: skip, or align against the full query-holes set (full)"),
        )
        .arg(
            Arg::new("top_k")
                .long("top-k")
                .value_parser(value_parser!(usize))
                .help("Keep only the top-K query holes per target hole by shared k-mers; default: no limit"),
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

/// Execute the rest command.
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

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_align_rest_")?;
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

    // 2. Convert to 2bit for random-access hole extraction.
    run_cmd!(info "==> fa to-2bit")?;
    let abs_target_2bit = common::to_2bit(&pgr, &abs_target, "target.2bit")?;
    let abs_query_2bit = common::to_2bit(&pgr, &abs_query, "query.2bit")?;

    // 3. Complement fill.
    let opts = RestOptions {
        trim: *args.get_one::<i32>("trim").unwrap(),
        min_anchor: *args.get_one::<i32>("min_anchor").unwrap(),
        max_gap: args.get_one::<i32>("max_gap").copied(),
        sampler: args.get_one::<String>("sampler").unwrap().clone(),
        smer: *args.get_one::<usize>("smer").unwrap(),
        kmer: *args.get_one::<usize>("kmer").unwrap(),
        window: *args.get_one::<usize>("window").unwrap(),
        min_shared: *args.get_one::<usize>("min_shared").unwrap(),
        unmatched: args.get_one::<String>("unmatched").unwrap().clone(),
        top_k: args.get_one::<usize>("top_k").copied(),
        preset: args.get_one::<String>("preset").cloned(),
        query_depth: *args.get_one::<usize>("query_depth").unwrap(),
        lastz_args: args.get_one::<String>("lastz_args").cloned(),
        parallel: *args.get_one::<usize>("parallel").unwrap(),
    };
    let combined = run_rest(
        &pgr,
        &abs_target_2bit,
        &abs_query_2bit,
        &avail_psl,
        &opts,
        &std::env::current_dir()?,
    )?;

    // 4. Write the combined PSL (anchors + LASTZ complement records).
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

// --- Complement-fill logic ---

/// A target-side hole, in 0-based half-open coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Hole {
    t_name: String,
    t_start: i32,
    t_end: i32,
}

/// Options controlling the complement-fill run.
#[derive(Debug, Clone)]
struct RestOptions {
    /// Shrink every anchor span by this many bp on both ends.
    trim: i32,
    /// Excise anchors shorter than this (bp) before the complement.
    min_anchor: i32,
    /// Holes longer than this (bp) are skipped; None means no limit.
    max_gap: Option<i32>,
    /// Fragment prefilter sampler: "syncmer" or "minimizer".
    sampler: String,
    /// Syncmer s-mer length (used when sampler = "syncmer").
    smer: usize,
    /// Minimizer k-mer length (used when sampler = "minimizer").
    kmer: usize,
    /// Prefilter sampler window size.
    window: usize,
    /// Minimum shared sampled k-mers to pair a target/query hole.
    min_shared: usize,
    /// Unpaired target holes: "skip" (default) or "full" (align against the
    /// whole query-holes set).
    unmatched: String,
    /// Keep only the top-K query holes per target hole by shared k-mers.
    top_k: Option<usize>,
    /// LASTZ preset (set01..set07); None uses lastz defaults.
    preset: Option<String>,
    /// Query depth threshold passed to lastz.
    query_depth: usize,
    /// Extra LASTZ arguments appended after the preset.
    lastz_args: Option<String>,
    /// Number of parallel lastz jobs.
    parallel: usize,
}

/// Build a per-contig IntSpan from the anchor PSL spans.
///
/// PSL coordinates are 0-based half-open; IntSpan is 1-based inclusive, so a
/// span `[start, end)` becomes `(start + 1, end)`.
fn psl_to_intspan(psls: &[Psl], is_target: bool) -> BTreeMap<String, IntSpan> {
    let mut set: BTreeMap<String, IntSpan> = BTreeMap::new();
    for p in psls {
        let (name, start, end) = if is_target {
            (&p.t_name, p.t_start, p.t_end)
        } else {
            (&p.q_name, p.q_start, p.q_end)
        };
        if end <= start {
            continue;
        }
        set.entry(name.clone())
            .or_default()
            .add_pair(start + 1, end);
    }
    set
}

/// Apply trim -> excise -> whole-genome complement on one side.
///
/// Returns the holes as 0-based half-open intervals, sorted by contig then
/// start. `sizes` supplies the contig lengths used for the whole-genome
/// complement (the `IntSpan::holes` builtin only covers internal gaps).
fn compute_rest_holes(
    set: &BTreeMap<String, IntSpan>,
    sizes: &BTreeMap<String, i32>,
    trim: i32,
    min_anchor: i32,
    max_gap: Option<i32>,
) -> Vec<Hole> {
    let trimmed = span_op(set, SpanOp::Trim, trim);
    let excised = span_op(&trimmed, SpanOp::Excise, min_anchor);

    let mut holes = Vec::new();
    for (chr, is) in &excised {
        let size = sizes.get(chr).copied().unwrap_or(0);
        if size <= 0 {
            continue;
        }
        let full = IntSpan::from_pair(1, size);
        let complement = full.diff(is);
        for (s, e) in complement.spans() {
            let start0 = s - 1;
            let end0 = e;
            if let Some(mg) = max_gap {
                if end0 - start0 > mg {
                    continue;
                }
            }
            holes.push(Hole {
                t_name: chr.clone(),
                t_start: start0,
                t_end: end0,
            });
        }
    }
    holes.sort_by(|a, b| (&a.t_name, a.t_start).cmp(&(&b.t_name, b.t_start)));
    holes
}

/// Extract every hole into an individual single-record FASTA in parallel,
/// returning the file paths in hole order.
fn extract_holes(
    pgr: &str,
    two_bit: &str,
    holes: &[Hole],
    prefix: &str,
    tempdir: &Path,
    parallel: usize,
) -> Result<Vec<PathBuf>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(parallel)
        .build()?;
    let paths: Vec<Result<PathBuf>> = pool.install(|| {
        use rayon::prelude::*;
        holes
            .par_iter()
            .enumerate()
            .map(|(i, h)| {
                let fa = tempdir.join(format!("{prefix}_{i}.fa"));
                common::extract_2bit_range(pgr, two_bit, &h.t_name, h.t_start, h.t_end, &fa)?;
                Ok(fa)
            })
            .collect()
    });
    paths.into_iter().collect()
}

/// Read a single-record FASTA file into its sequence bytes (header dropped).
fn read_fasta(path: &Path) -> Result<Vec<u8>> {
    let text = std::fs::read_to_string(path)?;
    let mut seq = Vec::new();
    for line in text.lines().skip(1) {
        seq.extend_from_slice(line.trim().as_bytes());
    }
    Ok(seq)
}

/// Sample a hole sequence into a k-mer hash set (syncmer or minimizer).
fn sample_hole(path: &Path, opts: &RestOptions) -> Result<rapidhash::RapidHashSet<u64>> {
    let seq = read_fasta(path)?;
    match opts.sampler.as_str() {
        "minimizer" => pgr::libs::hash::seq_mins(&seq, "rapid", opts.kmer, opts.window),
        _ => pgr::libs::syncmer::seq_syncmer_set(
            &seq,
            &pgr::libs::syncmer::SyncmerParams {
                smer: opts.smer,
                window: opts.window,
                seed: 7,
            },
            false,
        ),
    }
}

/// Pair target/query holes by shared sampled k-mers.
///
/// A pair `(i, j)` is kept when `t_sets[i] ∩ q_sets[j]` has at least
/// `min_shared` elements; unpaired target holes are skipped (likely
/// strain-specific sequence with no candidate homology on the query side).
fn pair_holes(
    t_sets: &[rapidhash::RapidHashSet<u64>],
    q_sets: &[rapidhash::RapidHashSet<u64>],
    min_shared: usize,
    top_k: Option<usize>,
) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for (i, ts) in t_sets.iter().enumerate() {
        let mut scored: Vec<(usize, usize)> = q_sets
            .iter()
            .enumerate()
            .map(|(j, qs)| (j, ts.intersection(qs).count()))
            .filter(|&(_, c)| c >= min_shared)
            .collect();
        scored.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        if let Some(k) = top_k {
            scored.truncate(k);
        }
        for (j, _) in scored {
            pairs.push((i, j));
        }
    }
    pairs
}

/// Run the complement-fill pipeline over two genomes.
fn run_rest(
    pgr: &str,
    target_2bit: &str,
    query_2bit: &str,
    avail_psl: &[Psl],
    opts: &RestOptions,
    tempdir: &Path,
) -> Result<Vec<Psl>> {
    let t_sizes = common::read_2bit_sizes(pgr, target_2bit)?;
    let q_sizes = common::read_2bit_sizes(pgr, query_2bit)?;
    let mut sizes = t_sizes.clone();
    sizes.extend(q_sizes.clone());

    // 1D per-side pipelines (no 2D coordinates).
    let t_set = psl_to_intspan(avail_psl, true);
    let q_set = psl_to_intspan(avail_psl, false);
    let t_holes = compute_rest_holes(&t_set, &t_sizes, opts.trim, opts.min_anchor, opts.max_gap);
    let q_holes = compute_rest_holes(&q_set, &q_sizes, opts.trim, opts.min_anchor, opts.max_gap);
    log::info!(
        "rest holes: {} target / {} query (trim={}, min-anchor={})",
        t_holes.len(),
        q_holes.len(),
        opts.trim,
        opts.min_anchor
    );
    if t_holes.is_empty() || q_holes.is_empty() {
        return Ok(avail_psl.to_vec());
    }

    let (common_args, _matrix_handle) = common::build_common_args(
        opts.preset.as_deref(),
        opts.query_depth,
        opts.lastz_args.as_deref(),
    )?;

    std::fs::create_dir_all(tempdir)?;
    // Extract both hole sets in parallel (individual files per hole).
    let t_extract = std::time::Instant::now();
    let t_fas = extract_holes(pgr, target_2bit, &t_holes, "t", tempdir, opts.parallel)?;
    let q_fas = extract_holes(pgr, query_2bit, &q_holes, "q", tempdir, opts.parallel)?;
    // Concatenate the query holes into one multi-sequence file, used when
    // unpaired target holes fall back to the full query-holes set.
    let q_all = tempdir.join("query_holes.fa");
    {
        use std::io::Write;
        let mut writer = pgr::libs::io::writer(q_all.to_str().unwrap_or(""))?;
        for qf in &q_fas {
            let text = std::fs::read_to_string(qf)?;
            writer.write_all(text.as_bytes())?;
        }
        writer.flush()?;
    }
    log::info!(
        "extracted {} target + {} query holes in {:?}",
        t_fas.len(),
        q_fas.len(),
        t_extract.elapsed()
    );

    let mut jobs: Vec<(PathBuf, PathBuf, PathBuf)> = Vec::new();
    if opts.sampler == "none" {
        // Full path: every target hole against the merged query-holes set.
        jobs.extend(t_fas.iter().enumerate().map(|(i, tf)| {
            (
                tf.clone(),
                q_all.clone(),
                tempdir.join(format!("hole_full_{i}.lav")),
            )
        }));
    } else {
        // Prefilter: sample every hole and pair by shared k-mers.
        let t_sample = std::time::Instant::now();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(opts.parallel)
            .build()?;
        let t_sets: Vec<Result<rapidhash::RapidHashSet<u64>>> = pool.install(|| {
            use rayon::prelude::*;
            t_fas.par_iter().map(|p| sample_hole(p, opts)).collect()
        });
        let t_sets: Vec<_> = t_sets.into_iter().collect::<Result<_>>()?;
        let q_sets: Vec<Result<rapidhash::RapidHashSet<u64>>> = pool.install(|| {
            use rayon::prelude::*;
            q_fas.par_iter().map(|p| sample_hole(p, opts)).collect()
        });
        let q_sets: Vec<_> = q_sets.into_iter().collect::<Result<_>>()?;
        let pairs = pair_holes(&t_sets, &q_sets, opts.min_shared, opts.top_k);
        log::info!(
            "prefilter ({} sampler, min-shared={}): {} pairs in {:?}",
            opts.sampler,
            opts.min_shared,
            pairs.len(),
            t_sample.elapsed()
        );
        if pairs.is_empty() && opts.unmatched != "full" {
            return Ok(avail_psl.to_vec());
        }
        jobs.extend(pairs.iter().enumerate().map(|(k, &(i, j))| {
            (
                t_fas[i].clone(),
                q_fas[j].clone(),
                tempdir.join(format!("hole_{k}.lav")),
            )
        }));
        if opts.unmatched == "full" {
            let mut matched: std::collections::HashSet<usize> = std::collections::HashSet::new();
            for &(i, _) in &pairs {
                matched.insert(i);
            }
            let mut k = pairs.len();
            for (i, _) in t_holes.iter().enumerate() {
                if matched.contains(&i) {
                    continue;
                }
                let t_fa = t_fas[i].clone();
                let lav = tempdir.join(format!("hole_unmatched_{k}.lav"));
                jobs.push((t_fa, q_all.clone(), lav));
                k += 1;
            }
            log::info!(
                "unmatched target holes (full fallback): {}",
                k - pairs.len()
            );
        }
    }

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
fn tpsl(t_name: &str, t_start: i32, t_end: i32, q_name: &str, q_start: i32, q_end: i32) -> Psl {
    Psl {
        t_name: t_name.into(),
        t_start,
        t_end,
        q_name: q_name.into(),
        q_start,
        q_end,
        strand: "+".into(),
        ..Default::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sizes_map(pairs: &[(&str, i32)]) -> BTreeMap<String, i32> {
        pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    #[test]
    fn psl_to_intspan_converts_to_1based() {
        let psls = vec![tpsl("t", 10, 20, "q", 30, 40)];
        let set = psl_to_intspan(&psls, true);
        assert_eq!(set["t"].to_string(), "11-20");
        let set = psl_to_intspan(&psls, false);
        assert_eq!(set["q"].to_string(), "31-40");
    }

    #[test]
    fn compute_rest_holes_whole_genome_complement_includes_ends() {
        // Anchors [100, 200) and [500, 600) on a 1000 bp contig; trim 0 and
        // no excise -> holes [1..99], [201..499], [601..1000] in 0-based
        // half-open: [0,100), [200,500), [600,1000).
        let psls = vec![
            tpsl("t", 100, 200, "q", 100, 200),
            tpsl("t", 500, 600, "q", 500, 600),
        ];
        let set = psl_to_intspan(&psls, true);
        let sizes = sizes_map(&[("t", 1000)]);
        let holes = compute_rest_holes(&set, &sizes, 0, 0, None);
        let coords: Vec<(String, i32, i32)> = holes
            .iter()
            .map(|h| (h.t_name.clone(), h.t_start, h.t_end))
            .collect();
        assert_eq!(
            coords,
            vec![
                ("t".to_string(), 0, 100),
                ("t".to_string(), 200, 500),
                ("t".to_string(), 600, 1000),
            ]
        );
    }

    #[test]
    fn compute_rest_holes_trim_and_excise_apply_before_complement() {
        // Trim 40 bp off each anchor end, then excise anchors shorter than
        // 30 bp. Anchor [100,200) trims to [140,160) (20 bp -> excised);
        // anchor [500,600) trims to [540,560) (20 bp -> excised). With no
        // anchors left, the complement is the whole contig.
        let psls = vec![
            tpsl("t", 100, 200, "q", 100, 200),
            tpsl("t", 500, 600, "q", 500, 600),
        ];
        let set = psl_to_intspan(&psls, true);
        let sizes = sizes_map(&[("t", 1000)]);
        let holes = compute_rest_holes(&set, &sizes, 40, 30, None);
        assert_eq!(holes.len(), 1);
        assert_eq!((holes[0].t_start, holes[0].t_end), (0, 1000));
    }

    #[test]
    fn compute_rest_holes_max_gap_skips_long_holes() {
        let psls = vec![tpsl("t", 100, 200, "q", 100, 200)];
        let set = psl_to_intspan(&psls, true);
        let sizes = sizes_map(&[("t", 10_000)]);
        // Tail hole ~9800 bp > 1000 -> skipped entirely.
        let holes = compute_rest_holes(&set, &sizes, 0, 0, Some(1000));
        assert_eq!(holes.len(), 1);
        assert_eq!((holes[0].t_start, holes[0].t_end), (0, 100));
    }

    #[test]
    fn compute_rest_holes_oversized_trim_drops_the_anchor() {
        // Trim 200 bp off a 100 bp anchor: inset leaves an empty span, the
        // anchor disappears and the complement is the whole contig.
        let psls = vec![tpsl("t", 100, 200, "q", 100, 200)];
        let set = psl_to_intspan(&psls, true);
        let sizes = sizes_map(&[("t", 1000)]);
        let holes = compute_rest_holes(&set, &sizes, 200, 0, None);
        assert_eq!(holes.len(), 1);
        assert_eq!((holes[0].t_start, holes[0].t_end), (0, 1000));
    }
}
