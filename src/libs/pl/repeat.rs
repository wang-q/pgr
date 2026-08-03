//! Repeat-identification pipeline drivers (FastK → Profex → spanr).

use cmd_lib::run_cmd;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

/// Run `Profex -z genome` per chromosome and write `.rg` files.
///
/// For each chromosome, runs `Profex -z genome <sn>` writing `prof.<sn>.txt`,
/// then scans lines with `re_prof` capturing `start` and `end` (1-based inclusive
/// in output). If `min_depth` is set and the regex has a `depth` capture group,
/// entries with depth below the threshold are skipped. Returns the list of
/// `prof.<sn>.rg` file names.
pub fn run_profex_per_chr(
    chrs: &[String],
    re_prof: &regex::Regex,
    min_depth: Option<usize>,
) -> anyhow::Result<Vec<String>> {
    let mut rg_files = vec![];
    for (i, chr) in chrs.iter().enumerate() {
        let sn = i + 1;
        run_cmd!(
            Profex -z genome ${sn} > prof.${sn}.txt
        )?;

        let reader = crate::reader(&format!("prof.{}.txt", sn))?;

        let rg_file = format!("prof.{}.rg", sn);
        let mut writer = crate::writer(&rg_file)?;

        for line in std::io::BufReader::new(reader)
            .lines()
            .map_while(Result::ok)
        {
            let Some(caps) = re_prof.captures(&line) else {
                continue;
            };

            if let Some(min_d) = min_depth {
                if let Some(depth_str) = caps.name("depth") {
                    let depth: usize = depth_str.as_str().parse()?;
                    if depth < min_d {
                        continue;
                    }
                }
            }

            let start = caps["start"].parse::<usize>()? + 1;
            let end = caps["end"].parse::<usize>()? + 1;

            writer.write_fmt(format_args!("{}:{}-{}\n", chr, start, end))?;
        }
        rg_files.push(rg_file);
    }
    Ok(rg_files)
}

/// Options for the shared repeat-identification pipeline (ir/rept).
pub struct RepeatOpts {
    /// Absolute path to the `pgr` executable.
    pub pgr: String,
    /// Absolute path to the genome FASTA.
    pub abs_infile: String,
    /// Absolute path to the output (or `stdout`).
    pub abs_outfile: String,
    pub opt_kmer: usize,
    pub opt_fk: usize,
    pub opt_min: usize,
    pub opt_ff: usize,
    /// For `ir`: absolute path to the repeat database. `None` for `rept`.
    pub abs_repeat: Option<String>,
    /// Keep the FastK repeat table (`repeat.ktab`) next to the library for
    /// reuse on later runs (`--keep-index`).
    pub keep_index: bool,
    /// Profex output regex (captures `start`/`end`, optionally `depth`).
    pub re_prof: regex::Regex,
    /// Minimum depth filter; `None` to skip. `Some(2)` for `rept`.
    pub min_depth: Option<usize>,
}

/// Run the shared FastK → Profex → spanr repeat pipeline.
///
/// When `opts.abs_repeat` is set, runs FastK twice (repeat + genome with
/// `-p:repeat`); otherwise runs FastK once on the genome (`-p`). Then
/// generates `chr.sizes`, runs Profex per chromosome, and finally the
/// spanr cover/fill/excise/fill pipeline.
pub fn run_repeat_pipeline(opts: &RepeatOpts) -> anyhow::Result<()> {
    let pgr = &opts.pgr;
    let abs_infile = &opts.abs_infile;
    let opt_kmer = opts.opt_kmer;

    if let Some(abs_repeat) = &opts.abs_repeat {
        // Cache the FastK table built from the repeat library next to the
        // library (`<lib>.repeat.k<k>.ktab`) when `keep_index` is set, and
        // reuse it on later runs as long as the library has not changed.
        let cache_prefix = format!("{}.repeat.k{}", abs_repeat, opt_kmer);
        if opts.keep_index && cache_is_fresh(abs_repeat, &cache_prefix) {
            run_cmd!(info "==> FastK on genome (reused repeat table)")?;
            run_cmd!(
                FastK -p:${cache_prefix} -k${opt_kmer} -Ngenome ${abs_infile}
            )?;
        } else {
            run_cmd!(info "==> FastK on repeat")?;
            run_cmd!(
                FastK -t -k${opt_kmer} -Nrepeat ${abs_repeat}
            )?;
            if opts.keep_index {
                let cache_path = format!("{}.ktab", cache_prefix);
                if let Err(e) = save_repeat_cache("repeat", &cache_prefix) {
                    log::warn!("failed to cache repeat table at {}: {}", cache_path, e);
                }
            }
            run_cmd!(info "==> FastK on genome")?;
            run_cmd!(
                FastK -p:repeat -k${opt_kmer} -Ngenome ${abs_infile}
            )?;
        }
    } else {
        run_cmd!(info "==> FastK")?;
        run_cmd!(
            FastK -p -k${opt_kmer} -Ngenome ${abs_infile}
        )?;
    }

    run_cmd!(info "==> Process each chromosome")?;
    run_cmd!(
        ${pgr} fa size ${abs_infile} -o chr.sizes
    )?;
    let chrs = crate::libs::io::read_names::<Vec<String>>("chr.sizes")?;

    let rg_files = run_profex_per_chr(&chrs, &opts.re_prof, opts.min_depth)?;

    run_repeat_spanr_pipeline(
        &rg_files,
        opts.opt_fk,
        opts.opt_min,
        opts.opt_ff,
        &opts.abs_outfile,
    )?;

    Ok(())
}

/// Options for the alignment-based repeat pipeline (`pgr rept e-align`).
pub struct AlignRepeatOpts {
    /// Absolute path to the `pgr` executable.
    pub pgr: String,
    /// Absolute path to the repeat library FASTA (query).
    pub abs_repeat: String,
    /// Absolute path to the genome FASTA (reference).
    pub abs_infile: String,
    /// Absolute path to the output (or `stdout`).
    pub abs_outfile: String,
    /// Keep the built `.pgi` indexes next to the inputs for reuse.
    pub keep_index: bool,
    pub kmer: usize,
    pub smer: usize,
    pub window: usize,
    pub freq: usize,
    pub min_span: usize,
    pub max_gap: usize,
    pub band: usize,
    pub merge_gap: usize,
    pub min_shared: usize,
    pub workflow: String,
    /// Minimum alignment identity (fraction of aligned bases matching).
    pub min_identity: f64,
    /// Minimum length of repetitive fragments (bp).
    pub min_len: usize,
    /// Fill holes between repetitive fragments (bp).
    pub fill_fragment: usize,
    /// Number of threads for the alignment.
    pub parallel: usize,
}

/// Run the `pgr align pgi` → PSL filter → spanr repeat pipeline.
///
/// The genome is the reference (PSL target) and the repeat library is the
/// query. Alignment blocks are filtered by identity and target-span length,
/// written as target-side `.rg`, then merged with the spanr pipeline.
pub fn run_align_repeat_pipeline(opts: &AlignRepeatOpts) -> anyhow::Result<()> {
    let pgr = &opts.pgr;
    let abs_infile = &opts.abs_infile;
    let abs_repeat = &opts.abs_repeat;
    let kmer = opts.kmer;
    let smer = opts.smer;
    let window = opts.window;
    let freq = opts.freq;
    let min_span = opts.min_span;
    let max_gap = opts.max_gap;
    let band = opts.band;
    let merge_gap = opts.merge_gap;
    let min_shared = opts.min_shared;
    let workflow = &opts.workflow;
    let parallel = opts.parallel;
    let keep_args = if opts.keep_index { "--keep-index" } else { "" };
    run_cmd!(info "==> Align repeats vs genome")?;
    run_cmd!(
        ${pgr} align pgi ${abs_infile} ${abs_repeat}
            -k ${kmer} --smer ${smer} --window ${window}
            -f ${freq} -c ${min_span} -s ${max_gap}
            --band ${band} --merge-gap ${merge_gap}
            --min-shared ${min_shared} --workflow ${workflow}
            -p ${parallel} ${keep_args} -o hits.psl
    )?;

    run_cmd!(info "==> Filter alignments")?;
    let reader = crate::reader("hits.psl")?;
    let mut writer = crate::writer("hits.rg")?;
    for line in std::io::BufReader::new(reader)
        .lines()
        .map_while(Result::ok)
    {
        if line.is_empty() {
            continue;
        }
        let Ok(psl) = line.parse::<crate::libs::fmt::psl::Psl>() else {
            continue;
        };
        let span = (psl.t_end - psl.t_start) as usize;
        if (psl.ident() as f64) < opts.min_identity || span < opts.min_len {
            continue;
        }
        writer.write_fmt(format_args!(
            "{}:{}-{}\n",
            psl.t_name,
            psl.t_start + 1,
            psl.t_end
        ))?;
    }
    drop(writer);

    run_repeat_spanr_pipeline(
        &["hits.rg".to_string()],
        0,
        opts.min_len,
        opts.fill_fragment,
        &opts.abs_outfile,
    )?;

    Ok(())
}

/// True when a complete cache for `cache_prefix` exists and is not older than
/// `lib` (library unchanged). Completeness is guaranteed by the `.complete`
/// marker written after all table files are copied.
fn cache_is_fresh(lib: &str, cache_prefix: &str) -> bool {
    let Ok(lib_meta) = std::fs::metadata(lib) else {
        return false;
    };
    for suffix in [".ktab", ".complete"] {
        let Ok(cache_meta) = std::fs::metadata(format!("{}{}", cache_prefix, suffix)) else {
            return false;
        };
        if !cache_meta.is_file() {
            return false;
        }
        if let (Ok(lib_m), Ok(cache_m)) = (lib_meta.modified(), cache_meta.modified()) {
            if cache_m < lib_m {
                return false;
            }
        }
    }
    true
}

/// Copy the freshly built FastK table (`<prefix>.ktab` plus its hidden part
/// files `.<prefix>.ktab.N`) to the cache path, renaming the prefix to the
/// cache prefix basename. `-p:` needs both the main table and the part files.
fn save_repeat_cache(src_prefix: &str, cache_prefix: &str) -> anyhow::Result<()> {
    let cache_path = format!("{}.ktab", cache_prefix);
    atomic_copy(&format!("{}.ktab", src_prefix), &cache_path)?;

    let base = Path::new(cache_prefix)
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow::anyhow!("invalid cache prefix: {}", cache_prefix))?;
    let dir = Path::new(cache_prefix)
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));

    for entry in std::fs::read_dir(".")?.flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        if let Some(rest) = name.strip_prefix(&format!(".{}.ktab.", src_prefix)) {
            let dst_name = format!(".{}.ktab.{}", base, rest);
            atomic_copy(&name, &dir.join(&dst_name).display().to_string())?;
        }
    }
    // Mark the cache complete only after every table file is in place.
    atomic_copy(
        &format!("{}.ktab", src_prefix),
        &format!("{}.complete", cache_prefix),
    )?;
    Ok(())
}

/// Atomically copy `src` to `dst` (write a temp file, then rename).
fn atomic_copy(src: &str, dst: &str) -> anyhow::Result<()> {
    let tmp = format!("{}.tmp.{}", dst, std::process::id());
    std::fs::copy(src, &tmp)?;
    std::fs::rename(&tmp, dst)?;
    Ok(())
}

/// Run the spanr cover → fill → excise → fill pipeline on `rg_files`.
pub fn run_repeat_spanr_pipeline(
    rg_files: &[String],
    fk: usize,
    min: usize,
    ff: usize,
    abs_outfile: &str,
) -> anyhow::Result<()> {
    run_cmd!(info "==> Outputs")?;
    run_cmd!(
        spanr cover $[rg_files] |
            spanr span --op fill -n ${fk} stdin |
            spanr span --op excise -n ${min} stdin |
            spanr span --op fill -n ${ff} stdin -o ${abs_outfile}
    )?;
    Ok(())
}

/// Parse a TRF `.dat` file and write `chr:start-end` lines to `writer`.
///
/// Each TRF row has at least 15 whitespace-separated fields; the first two
/// are 1-based start and end coordinates. Rows with fewer fields are skipped
/// with a `log::debug!` message.
pub fn parse_trf_output<R: BufRead, W: Write>(
    reader: R,
    chr: &str,
    writer: &mut W,
) -> anyhow::Result<()> {
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split_ascii_whitespace().collect();
        if fields.len() < 15 {
            log::debug!("skipping short TRF line: {}", line);
            continue;
        }

        let start = fields[0].parse::<usize>()?;
        let end = fields[1].parse::<usize>()?;

        writer.write_fmt(format_args!("{}:{}-{}\n", chr, start, end))?;
    }
    Ok(())
}
