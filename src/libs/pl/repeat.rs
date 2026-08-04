//! Repeat-identification pipeline drivers (FastK → Profex → runlist).

use cmd_lib::run_cmd;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

/// True when any base of the input FASTA is lowercase (soft-masked).
///
/// `pgr fa masked` also reports N/gap regions, so it cannot be used to
/// detect soft-masking specifically; scan the sequences directly instead.
fn has_soft_mask(infile: &str) -> anyhow::Result<bool> {
    let mut reader = crate::libs::fmt::fa::reader(infile)?;
    for result in reader.records() {
        let rec = result?;
        if rec
            .sequence()
            .as_ref()
            .iter()
            .any(|b| b.is_ascii_lowercase())
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Run `Profex -z genome` per chromosome and write `.rg` files.
///
/// For each chromosome, runs `Profex -z genome <sn>` writing `prof.<sn>.txt`,
/// then scans lines with `re_prof` capturing `start` and `end`. Profex prints
/// the 0-based k-mer start of each run and closes it with the 1-based inclusive
/// end (start + run length + kmer - 1), so the `.rg` output is 1-based
/// inclusive with `start + 1` and `end` as-is. If `min_depth` is set and the
/// regex has a `depth` capture group, entries with depth below the threshold
/// are skipped. `Profex -z` never closes the final run of a read (its end and
/// depth are omitted); when no depth threshold is applied (e.g. e-kmer) the
/// run is closed with the chromosome length from `lens`, and with a threshold
/// (e.g. s-kmer) it is conservatively dropped since its depth is unknown.
/// Returns the list of `prof.<sn>.rg` file names.
pub fn run_profex_per_chr(
    chrs: &[String],
    lens: &[usize],
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
        let mut tail_start: Option<usize> = None;

        for line in std::io::BufReader::new(reader).lines() {
            let line = line?;
            let Some(caps) = re_prof.captures(&line) else {
                // The final run of a read is printed as a bare start.
                if let Ok(start) = line.trim().parse::<usize>() {
                    tail_start = Some(start);
                }
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
            let end = caps["end"].parse::<usize>()?;

            writer.write_fmt(format_args!("{}:{}-{}\n", chr, start, end))?;
        }

        if let Some(start) = tail_start {
            if min_depth.is_none() {
                writer.write_fmt(format_args!("{}:{}-{}\n", chr, start + 1, lens[i]))?;
            }
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

/// Run the shared FastK → Profex → runlist repeat pipeline.
///
/// When `opts.abs_repeat` is set, runs FastK twice (repeat + genome with
/// `-p:repeat`); otherwise runs FastK once on the genome (`-p`). Then
/// generates `chr.sizes`, runs Profex per chromosome, and finally the
/// internal cover/fill/excise/fill runlist pipeline.
pub fn run_repeat_pipeline(opts: &RepeatOpts) -> anyhow::Result<()> {
    let pgr = &opts.pgr;
    let abs_infile = &opts.abs_infile;
    let opt_kmer = opts.opt_kmer;
    // FastK's block-level sort files go to a fixed global dir by default
    // (/tmp), so concurrent or repeated runs clobber each other's partial
    // tables (observed as SIGSEGV or corrupted profiles). Point -P at the
    // pipeline tempdir (the current working directory after enter()).
    let sort_dir = std::env::current_dir()?.display().to_string();

    if let Some(abs_repeat) = &opts.abs_repeat {
        // Cache the FastK table built from the repeat library next to the
        // library (`<lib>.repeat.k<k>.ktab`) when `keep_index` is set, and
        // reuse it on later runs as long as the library has not changed.
        let cache_prefix = format!("{}.repeat.k{}", abs_repeat, opt_kmer);
        if opts.keep_index && cache_is_fresh(abs_repeat, &cache_prefix) {
            run_cmd!(info "==> FastK on genome (reused repeat table)")?;
            run_cmd!(
                FastK -p:${cache_prefix} -k${opt_kmer} -Ngenome -P${sort_dir} ${abs_infile}
            )?;
        } else {
            run_cmd!(info "==> FastK on repeat")?;
            run_cmd!(
                FastK -t -k${opt_kmer} -Nrepeat -P${sort_dir} ${abs_repeat}
            )?;
            if opts.keep_index {
                let cache_path = format!("{}.ktab", cache_prefix);
                if let Err(e) = save_repeat_cache("repeat", &cache_prefix) {
                    log::warn!("failed to cache repeat table at {}: {}", cache_path, e);
                }
            }
            run_cmd!(info "==> FastK on genome")?;
            run_cmd!(
                FastK -p:repeat -k${opt_kmer} -Ngenome -P${sort_dir} ${abs_infile}
            )?;
        }
    } else {
        run_cmd!(info "==> FastK")?;
        run_cmd!(
            FastK -p -k${opt_kmer} -Ngenome -P${sort_dir} ${abs_infile}
        )?;
    }

    run_cmd!(info "==> Process each chromosome")?;
    run_cmd!(
        ${pgr} fa size ${abs_infile} -o chr.sizes
    )?;
    let mut chrs = Vec::new();
    let mut lens = Vec::new();
    for line in crate::libs::io::read_lines("chr.sizes")? {
        let mut fields = line.split_whitespace();
        if let (Some(name), Some(len)) = (fields.next(), fields.next()) {
            chrs.push(name.to_string());
            lens.push(len.parse()?);
        }
    }

    // The runlist parser truncates dotted contig names (e.g. `NC_000913.1`
    // -> `1`) at the last '.', so map real names to dot-free placeholders
    // and restore them after the runlist pass.
    let mut name_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut safe_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let safe_chrs: Vec<String> = chrs
        .iter()
        .map(|c| {
            let s = format!("c{}", name_map.len() + 1);
            name_map.insert(c.clone(), s.clone());
            safe_map.insert(s.clone(), c.clone());
            s
        })
        .collect();

    let rg_files = run_profex_per_chr(&safe_chrs, &lens, &opts.re_prof, opts.min_depth)?;

    if count_rg_lines(&rg_files)? == 0 {
        // No repetitive intervals: emit an empty runlist directly.
        let empty = b"{}\n";
        if opts.abs_outfile == "stdout" {
            std::io::stdout().write_all(empty)?;
        } else {
            std::fs::write(&opts.abs_outfile, empty)?;
        }
        return Ok(());
    }

    run_repeat_runlist_pipeline(
        &rg_files,
        opts.opt_fk,
        opts.opt_min,
        opts.opt_ff,
        "out.json",
    )?;

    // Restore the real contig names in the runlist json.
    let mut val: serde_json::Value = serde_json::from_slice(&std::fs::read("out.json")?)?;
    if let Some(obj) = val.as_object_mut() {
        let old = std::mem::take(obj);
        for (k, v) in old {
            // Drop the empty marker `-` so the runlist stays clean.
            if v.as_str() == Some("-") {
                continue;
            }
            obj.insert(safe_map.get(&k).cloned().unwrap_or(k), v);
        }
    }
    let out_bytes = serde_json::to_vec_pretty(&val)?;
    if opts.abs_outfile == "stdout" {
        let mut w = crate::writer("stdout")?;
        w.write_all(&out_bytes)?;
        w.write_all(b"\n")?;
    } else {
        std::fs::write(&opts.abs_outfile, out_bytes)?;
    }

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

/// Run the `pgr align pgi` → PSL filter → runlist repeat pipeline.
///
/// The genome is the reference (PSL target) and the repeat library is the
/// query. Alignment blocks are filtered by identity and target-span length,
/// written as target-side `.rg`, then merged with the runlist pipeline.
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

    // Soft-masked (lowercase) repeats fragment pgi's chain extension, so the
    // alignment pass massively underestimates coverage. Detect and warn
    // instead of silently returning bad numbers.
    if has_soft_mask(abs_infile)? {
        log::warn!(
            "input genome contains soft-masked (lowercase) regions; e-align \
             results will be underestimated, consider uppercasing first \
             (`tr a-z A-Z`)"
        );
    }

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
    // The runlist parser truncates dotted contig names (e.g. `NC_000913.1`
    // -> `1`) at the last '.', so map real names to dot-free placeholders
    // and restore them after the runlist pass.
    let mut name_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut safe_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut n_rg = 0usize;
    for line in std::io::BufReader::new(reader).lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let psl = match line.parse::<crate::libs::fmt::psl::Psl>() {
            Ok(p) => p,
            Err(e) => {
                log::warn!("skipping unparseable psl line: {}: {}", line, e);
                continue;
            }
        };
        // Guard against a malformed record with t_end < t_start (a negative
        // difference would wrap into a huge span and pass the length filter).
        // i64 arithmetic: extreme PSL coordinates would overflow the i32
        // subtraction (e.g. t_start = i32::MIN, t_end = i32::MAX).
        let span = (psl.t_end as i64 - psl.t_start as i64).max(0) as usize;
        if (psl.ident() as f64) < opts.min_identity || span < opts.min_len {
            continue;
        }
        let safe = match name_map.get(&psl.t_name) {
            Some(s) => s.clone(),
            None => {
                let s = format!("c{}", name_map.len() + 1);
                name_map.insert(psl.t_name.clone(), s.clone());
                safe_map.insert(s.clone(), psl.t_name.clone());
                s
            }
        };
        writer.write_fmt(format_args!("{}:{}-{}\n", safe, psl.t_start + 1, psl.t_end))?;
        n_rg += 1;
    }
    drop(writer);

    if n_rg == 0 {
        // No alignments survived the filters: emit an empty runlist directly.
        let empty = b"{}\n";
        if opts.abs_outfile == "stdout" {
            std::io::stdout().write_all(empty)?;
        } else {
            std::fs::write(&opts.abs_outfile, empty)?;
        }
        return Ok(());
    }

    run_repeat_runlist_pipeline(
        &["hits.rg".to_string()],
        0,
        opts.min_len,
        opts.fill_fragment,
        "out.json",
    )?;

    // Restore the real contig names in the runlist json.
    let mut val: serde_json::Value = serde_json::from_slice(&std::fs::read("out.json")?)?;
    if let Some(obj) = val.as_object_mut() {
        let old = std::mem::take(obj);
        for (k, v) in old {
            // Drop the empty marker `-` so the runlist stays clean.
            if v.as_str() == Some("-") {
                continue;
            }
            obj.insert(safe_map.get(&k).cloned().unwrap_or(k), v);
        }
    }
    let out_bytes = serde_json::to_vec_pretty(&val)?;
    if opts.abs_outfile == "stdout" {
        let mut w = crate::writer("stdout")?;
        w.write_all(&out_bytes)?;
        w.write_all(b"\n")?;
    } else {
        std::fs::write(&opts.abs_outfile, out_bytes)?;
    }

    Ok(())
}

/// Options for the self-alignment repeat pipeline (`pgr rept s-align`).
pub struct SelfAlignOpts {
    /// Absolute path to the `pgr` executable.
    pub pgr: String,
    /// Absolute path to the genome FASTA.
    pub abs_infile: String,
    /// Absolute path to the output (or `stdout`).
    pub abs_outfile: String,
    /// Overlapping window length (bp).
    pub window: usize,
    /// Window step size (bp).
    pub step: usize,
    /// Split window output into chunks of N records.
    pub chunk_records: usize,
    /// lastz preset name.
    pub preset: String,
    /// Number of threads for the alignment.
    pub parallel: usize,
    /// Minimum alignment depth for a region to be kept.
    pub min_depth: usize,
}

/// Run the Cactus-style self-alignment repeat pipeline (`pgr-repeat.sh`):
/// window the genome, align the windows back to the genome with lastz, lift
/// to genomic coordinates, and keep regions whose alignment depth exceeds a
/// threshold (baseline 2x from 50%-overlap windows; >= 4 means >= 2 copies).
pub fn run_self_align_pipeline(opts: &SelfAlignOpts) -> anyhow::Result<()> {
    let pgr = &opts.pgr;
    let abs_infile = &opts.abs_infile;
    let abs_outfile = &opts.abs_outfile;
    let window = opts.window;
    let step = opts.step;
    let chunk_records = opts.chunk_records;
    let preset = &opts.preset;
    let parallel = opts.parallel;
    let min_depth = opts.min_depth;

    // Soft-masked (lowercase) repeats are skipped by lastz, so the pass
    // underestimates coverage; detect and warn instead of silent bad data.
    if has_soft_mask(abs_infile)? {
        log::warn!(
            "input genome contains soft-masked (lowercase) regions; self \
             alignment results will be underestimated, consider uppercasing \
             first (`tr a-z A-Z`)"
        );
    }

    run_cmd!(info "==> Windowing")?;
    std::fs::create_dir_all("fragments")?;
    run_cmd!(
        ${pgr} fa window ${abs_infile} -w ${window} --step ${step}
            --chunk-records ${chunk_records} -o fragments/fragments.fa
    )?;

    run_cmd!(info "==> Split genome by name")?;
    run_cmd!(
        ${pgr} fa split name ${abs_infile} -o genome
    )?;

    run_cmd!(info "==> Align windows to genome (lastz)")?;
    run_cmd!(
        ${pgr} align lastz genome fragments --preset ${preset}
            --parallel ${parallel} -o lastz_out
    )?;

    run_cmd!(info "==> Convert LAV to PSL")?;
    let lav_files = crate::libs::io::list_files_ext("lastz_out", "lav");
    for lav in &lav_files {
        run_cmd!(${pgr} lav to-psl ${lav} >> fragments.psl)?;
    }

    run_cmd!(info "==> Lift to genomic coordinates")?;
    run_cmd!(
        ${pgr} fa size ${abs_infile} -o chrom.sizes
    )?;
    run_cmd!(
        ${pgr} psl lift fragments.psl --q-sizes chrom.sizes -o lifted.psl
    )?;

    run_cmd!(info "==> Extract ranges")?;
    run_cmd!(
        ${pgr} psl to-range lifted.psl -o coverage.rg
    )?;

    if count_rg_lines(&["coverage.rg".to_string()])? == 0 {
        let empty = b"{}\n";
        if abs_outfile == "stdout" {
            std::io::stdout().write_all(empty)?;
        } else {
            std::fs::write(abs_outfile, empty)?;
        }
        return Ok(());
    }

    // The runlist parser truncates dotted contig names (e.g. `NC_000913.1`
    // -> `1`) at the last '.', so map real names to dot-free placeholders
    // and restore them after the runlist pass (same convention as the other
    // runlist pipelines).
    let chrs = crate::libs::io::read_names::<Vec<String>>("chrom.sizes")?;
    let mut name_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut safe_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    for c in &chrs {
        let s = format!("c{}", name_map.len() + 1);
        name_map.insert(c.clone(), s.clone());
        safe_map.insert(s, c.clone());
    }
    let reader = crate::reader("coverage.rg")?;
    let mut writer = crate::writer("coverage.safe.rg")?;
    for line in std::io::BufReader::new(reader).lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let (name, rest) = line
            .split_once(':')
            .ok_or_else(|| anyhow::anyhow!("invalid .rg line: {}", line))?;
        let safe = name_map
            .get(name)
            .cloned()
            .unwrap_or_else(|| name.to_string());
        writer.write_fmt(format_args!("{}:{}\n", safe, rest))?;
    }
    drop(writer);

    run_cmd!(info "==> Coverage")?;
    let reader = crate::reader("coverage.safe.rg")?;
    let iv_of = crate::libs::runlist::rg_to_intervals(reader)?;
    let mut set: std::collections::BTreeMap<String, crate::libs::ds::IntSpan> =
        std::collections::BTreeMap::new();
    for (chr, ivs) in &iv_of {
        set.insert(
            chr.clone(),
            crate::libs::runlist::depth_at_least(ivs, min_depth as u32),
        );
    }
    let json = crate::libs::ds::intspan::set2json(&set);
    std::fs::write("out.json", serde_json::to_vec_pretty(&json)?)?;

    // Restore the real contig names in the runlist json.
    let mut val: serde_json::Value = serde_json::from_slice(&std::fs::read("out.json")?)?;
    if let Some(obj) = val.as_object_mut() {
        let old = std::mem::take(obj);
        for (k, v) in old {
            // Drop the empty marker `-` so the runlist stays clean.
            if v.as_str() == Some("-") {
                continue;
            }
            obj.insert(safe_map.get(&k).cloned().unwrap_or(k), v);
        }
    }
    let out_bytes = serde_json::to_vec_pretty(&val)?;
    if abs_outfile == "stdout" {
        let mut w = crate::writer("stdout")?;
        w.write_all(&out_bytes)?;
        w.write_all(b"\n")?;
    } else {
        std::fs::write(abs_outfile, out_bytes)?;
    }

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

/// Run the cover → fill → excise → fill pipeline on `rg_files`.
pub fn run_repeat_runlist_pipeline(
    rg_files: &[String],
    fk: usize,
    min: usize,
    ff: usize,
    abs_outfile: &str,
) -> anyhow::Result<()> {
    run_cmd!(info "==> Outputs")?;
    let set = crate::libs::runlist::rg_files_to_set(rg_files)?;
    // The original spanr pipeline ran `spanr span` three times; folding them
    // into sequential passes on the merged set gives identical results.
    let set = crate::libs::runlist::span_op(&set, crate::libs::runlist::SpanOp::Fill, fk as i32);
    let set = crate::libs::runlist::span_op(&set, crate::libs::runlist::SpanOp::Excise, min as i32);
    let set = crate::libs::runlist::span_op(&set, crate::libs::runlist::SpanOp::Fill, ff as i32);
    let mut res = std::collections::BTreeMap::new();
    res.insert("__single__".to_string(), set);
    crate::libs::runlist::write_sets(abs_outfile, &res)?;
    Ok(())
}

/// Count the total number of `chr:start-end` lines across the given `.rg`
/// files.
pub fn count_rg_lines(rg_files: &[String]) -> anyhow::Result<usize> {
    let mut n = 0usize;
    for rg in rg_files {
        let reader = crate::reader(rg)?;
        for line in std::io::BufReader::new(reader).lines() {
            if !line?.trim().is_empty() {
                n += 1;
            }
        }
    }
    Ok(n)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn soft_mask_detection_ignores_n_gaps() {
        // Regression: the s-align/e-align soft-mask warning used `pgr fa
        // masked`, which reports N/gap regions too, so a genome with N runs
        // but no lowercase bases warned about lowercase soft-masking.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("g.fa");
        std::fs::write(
            &path,
            format!(
                ">chr\n{}\n",
                "ACGT".repeat(50) + &"N".repeat(100) + &"ACGT".repeat(50)
            ),
        )
        .unwrap();
        assert!(
            !has_soft_mask(path.to_str().unwrap()).unwrap(),
            "N gaps must not count as soft-masking"
        );

        let lower = dir.path().join("lower.fa");
        std::fs::write(
            &lower,
            format!(">chr\n{}\n", "ACGT".repeat(50) + &"acgt".repeat(10)),
        )
        .unwrap();
        assert!(has_soft_mask(lower.to_str().unwrap()).unwrap());
    }
}
