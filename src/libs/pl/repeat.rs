//! Repeat-identification pipeline drivers (FastK → Profex → runlist).

use cmd_lib::run_cmd;
use rayon::prelude::*;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

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

/// True when the FASTA contains at least one record with a non-empty
/// sequence. FastK segfaults on sequence-less inputs (e.g. an empty repeat
/// library), so the pipeline rejects them up front with a friendly error.
fn has_sequences(infile: &str) -> anyhow::Result<bool> {
    let mut reader = crate::libs::fmt::fa::reader(infile)?;
    for result in reader.records() {
        if !result?.sequence().is_empty() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Find an executable by name in `$PATH`, or `None` when absent.
fn find_in_path(name: &str) -> Option<PathBuf> {
    std::env::var_os("PATH").and_then(|paths| {
        std::env::split_paths(&paths)
            .map(|dir| dir.join(name))
            .find(|p| p.is_file())
    })
}

/// Read the first FASTA record's sequence bytes.
fn read_fasta_sequence(path: &str) -> anyhow::Result<Vec<u8>> {
    let mut reader = crate::libs::fmt::fa::reader(path)?;
    if let Some(result) = reader.records().next() {
        let rec = result?;
        return Ok(rec.sequence().as_ref().to_vec());
    }
    anyhow::bail!("no FASTA records in {path}")
}

/// RepeatMasker `SimpleBatcher` fragment layout for one sequence: a list of
/// (0-based start, length) fragments with `overlap` bp between neighbours.
fn rm_batches(len: usize, frag: usize, overlap: usize) -> Vec<(usize, usize)> {
    if len <= frag {
        return vec![(0, len)];
    }
    let mut divisor = 2usize;
    while (len + (divisor - 1) * overlap) / divisor > frag {
        divisor += 1;
    }
    let mut size = (len + (divisor - 1) * overlap) / divisor;
    let mut batches = Vec::new();
    for i in 0..divisor - 1 {
        let start = i * (size - overlap);
        batches.push((start, size));
    }
    let rem = (len + (divisor - 1) * overlap) % divisor;
    size += rem;
    batches.push((len - size, size));
    batches
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
    anyhow::ensure!(
        has_sequences(abs_infile)?,
        "input genome FASTA has no sequences: {}",
        abs_infile
    );
    if let Some(abs_repeat) = &opts.abs_repeat {
        anyhow::ensure!(
            has_sequences(abs_repeat)?,
            "repeat library FASTA has no sequences: {}",
            abs_repeat
        );
    }
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
    pub min_shared: usize,
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
    let min_shared = opts.min_shared;
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
            -f ${freq} --min-shared ${min_shared}
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

/// Options for the RMBlast repeat pipeline (`pgr rept rmblast`), replicating
/// RepeatMasker 4.2.4's `-lib` search recipe.
pub struct RmblastOpts {
    /// Absolute path to the `pgr` executable.
    pub pgr: String,
    /// Absolute path to the repeat library FASTA (`.gz` accepted).
    pub abs_repeat: String,
    /// Absolute path to the genome FASTA (`.gz` accepted).
    pub abs_infile: String,
    /// Absolute path to the output (or `stdout`).
    pub abs_outfile: String,
    /// RepeatMasker `-cutoff` (default 225, passed to rmblastn unchanged).
    pub cutoff: i32,
    /// rmblastn `-word_size` (default 9).
    pub word_size: usize,
    /// Fixed GC percentage for matrix selection; `None` = per chromosome.
    pub matrix_gc: Option<i64>,
    /// Shortest fragment kept after merging (0 = RepeatMasker raw hits).
    pub min_len: usize,
    /// Fill holes between fragments (0 = RepeatMasker raw hits).
    pub fill_fragment: usize,
    /// Total threads across rmblastn processes (4 per process, like RM).
    pub parallel: usize,
    /// RepeatMasker `-frag`: max fragment length before splitting (0 = no
    /// fragmentation, whole chromosome per job).
    pub frag: usize,
    /// Directory containing `makeblastdb` / `rmblastn` (default `$PATH`).
    pub rmblast_dir: Option<PathBuf>,
}

/// Run the RMBlast repeat pipeline: split the genome by chromosome, build a
/// makeblastdb from the library, search each chromosome with the RepeatMasker
/// `general_search_parameters` recipe, and write a runlist JSON.
pub fn run_rmblast_repeat_pipeline(opts: &RmblastOpts) -> anyhow::Result<()> {
    let pgr = &opts.pgr;
    let cwd = std::env::current_dir()?;

    if has_soft_mask(&opts.abs_infile)? {
        log::warn!(
            "input genome contains soft-masked (lowercase) regions; rmblastn \
             will skip them, consider uppercasing first (`tr a-z A-Z`)"
        );
    }

    // Resolve the RMBlast binaries: `--rmblast-dir` overrides, otherwise
    // fall back to `$PATH` with a friendly error when missing.
    let makeblastdb = match &opts.rmblast_dir {
        Some(dir) => {
            let bin = dir.join("makeblastdb");
            anyhow::ensure!(
                bin.is_file(),
                "makeblastdb not found in {}; is the RMBlast directory correct?",
                dir.display()
            );
            bin
        }
        None => find_in_path("makeblastdb").ok_or_else(|| {
            anyhow::anyhow!(
                "makeblastdb not found in $PATH; install RMBlast or pass --rmblast-dir <dir>"
            )
        })?,
    };
    let rmblastn = match &opts.rmblast_dir {
        Some(dir) => {
            let bin = dir.join("rmblastn");
            anyhow::ensure!(
                bin.is_file(),
                "rmblastn not found in {}; is the RMBlast directory correct?",
                dir.display()
            );
            bin
        }
        None => find_in_path("rmblastn").ok_or_else(|| {
            anyhow::anyhow!(
                "rmblastn not found in $PATH; install RMBlast or pass --rmblast-dir <dir>"
            )
        })?,
    };

    // makeblastdb does not read gzipped FASTA; normalize the library to a
    // plain file in the tempdir (also keeps the db files out of the user dir).
    run_cmd!(info "==> Prepare repeat library")?;
    let lib_fa = "repeats.fa";
    if opts.abs_repeat.ends_with(".gz") {
        let mut reader = crate::reader(&opts.abs_repeat)?;
        let mut writer = crate::writer(lib_fa)?;
        std::io::copy(&mut reader, &mut writer)?;
    } else {
        std::fs::copy(&opts.abs_repeat, lib_fa)?;
    }

    run_cmd!(info "==> Build library database")?;
    let db_status = std::process::Command::new(&makeblastdb)
        .args(["-dbtype", "nucl", "-in", lib_fa, "-out", lib_fa])
        .output()?;
    if !db_status.status.success() {
        anyhow::bail!(
            "makeblastdb failed: {}",
            String::from_utf8_lossy(&db_status.stderr)
        );
    }

    // rmblastn resolves `-matrix <name>` through `BLASTMAT`; write all GC
    // matrices once and point every job at the same directory.
    run_cmd!(info "==> Write scoring matrices")?;
    let matrices_dir = cwd.join("matrices");
    std::fs::create_dir_all(&matrices_dir)?;
    for name in crate::libs::rmblast::MATRIX_GC_NAMES {
        let content = crate::libs::rmblast::matrix_content(name).expect("matrix list mismatch");
        std::fs::write(matrices_dir.join(format!("20p{name}.matrix")), content)?;
    }

    // Split the genome into per-chromosome files (sanitized names, .fa.gz ok).
    run_cmd!(info "==> Split genome by chromosomes")?;
    let abs_infile = &opts.abs_infile;
    run_cmd!(${pgr} fa size ${abs_infile} -o chr.sizes)?;
    let chrs = crate::libs::io::read_names::<Vec<String>>("chr.sizes")?;
    run_cmd!(${pgr} fa split name ${abs_infile} -o .)?;

    // The runlist parser truncates dotted contig names (e.g. `NC_000913.1` ->
    // `1`) at the last '.', so map real names to dot-free placeholders and
    // restore them after the runlist pass (same as trf / e-align).
    let mut safe_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();

    // Fragment each chromosome with RepeatMasker's `SimpleBatcher` layout
    // (fragmentLen/overlap), write one query file per batch, and pick the
    // GC-keyed matrix per batch (RepeatMasker uses per-batch GC when the
    // batch is a single sequence longer than 2000 bp, else 43).
    run_cmd!(info "==> Fragment genome (RepeatMasker batching)")?;
    let overlap = 2000usize;
    let mut jobs: Vec<(String, String, String, usize, &'static str)> = Vec::new();
    for (i, chr) in chrs.iter().enumerate() {
        let safe = format!("c{}", i + 1);
        safe_map.insert(safe.clone(), chr.clone());
        let chr_file = format!("{}.fa", crate::libs::io::sanitize_filename(chr));
        let seq = read_fasta_sequence(&chr_file)?;
        let batches = rm_batches(seq.len(), opts.frag, overlap);
        for (j, (start, len)) in batches.iter().enumerate() {
            let idx = jobs.len();
            let query = format!("batch.{idx}.fa");
            let out = format!("hits.{idx}.out");
            let mut w = crate::writer(&query)?;
            writeln!(w, ">c{}frag-{}", i + 1, j + 1)?;
            for chunk in seq[*start..*start + len].chunks(60) {
                w.write_all(chunk)?;
                w.write_all(b"\n")?;
            }
            drop(w);
            let matrix_name = match opts.matrix_gc {
                Some(gc) => crate::libs::rmblast::matrix_name_for_gc(gc),
                None if batches.len() == 1 && *len <= 2000 => "43g",
                None => crate::libs::rmblast::matrix_name_for_gc(crate::libs::rmblast::gc_bytes(
                    &seq[*start..*start + len],
                )),
            };
            jobs.push((safe.clone(), query, out, *start, matrix_name));
        }
    }

    let minscore = crate::libs::rmblast::effective_minscore(opts.cutoff);
    let n_jobs = jobs.len();
    let failures = AtomicUsize::new(0);
    let first_err: Mutex<Option<String>> = Mutex::new(None);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads((opts.parallel / 4).max(1))
        .build()?;

    run_cmd!(info "==> Search repeats vs genome (rmblastn)")?;
    let rg_results: Vec<anyhow::Result<Option<String>>> = pool.install(|| {
        jobs.par_iter()
            .enumerate()
            .map(
                |(i, (safe, query, out, start, matrix_name))| -> anyhow::Result<Option<String>> {
                    let rg = format!("hits.{i}.rg");
                    let args = crate::libs::rmblast::build_args(
                        lib_fa,
                        query,
                        matrix_name,
                        minscore,
                        opts.word_size,
                        out,
                        4,
                    );
                    let mut cmd = std::process::Command::new(&rmblastn);
                    cmd.args(&args);
                    cmd.env("BLASTMAT", &matrices_dir);
                    log::info!("rmblastn: {cmd:?}");

                    match cmd.output() {
                        Ok(o) if o.status.success() => {
                            let mut writer = crate::writer(&rg)?;
                            let reader = crate::reader(out)?;
                            for line in reader.lines() {
                                let line = line?;
                                if let Some((_, qstart, qend)) =
                                    crate::libs::rmblast::parse_tab_row(&line)
                                {
                                    writer.write_fmt(format_args!(
                                        "{}:{}-{}\n",
                                        safe,
                                        qstart + *start as i64,
                                        qend + *start as i64
                                    ))?;
                                }
                            }
                            Ok(Some(rg))
                        }
                        Ok(o) => {
                            let msg = String::from_utf8_lossy(&o.stderr).trim().to_string();
                            if !msg.is_empty() {
                                let mut guard = first_err.lock().unwrap();
                                if guard.is_none() {
                                    *guard = Some(msg);
                                }
                            }
                            failures.fetch_add(1, Ordering::Relaxed);
                            Ok(None)
                        }
                        Err(e) => {
                            log::error!("failed to spawn rmblastn for {query}: {e}");
                            failures.fetch_add(1, Ordering::Relaxed);
                            Ok(None)
                        }
                    }
                },
            )
            .collect()
    });
    let rg_files: Vec<String> = rg_results
        .into_iter()
        .collect::<anyhow::Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect();

    let n_failures = failures.load(Ordering::Relaxed);
    if n_failures > 0 {
        let detail = first_err.lock().unwrap().clone().unwrap_or_default();
        if detail.is_empty() {
            anyhow::bail!("rmblastn failed for {n_failures} of {n_jobs} jobs");
        } else {
            anyhow::bail!("rmblastn failed for {n_failures} of {n_jobs} jobs: {detail}");
        }
    }

    if rg_files.is_empty() {
        // No hits: emit an empty runlist directly.
        let empty = b"{}\n";
        if opts.abs_outfile == "stdout" {
            std::io::stdout().write_all(empty)?;
        } else {
            std::fs::write(&opts.abs_outfile, empty)?;
        }
        return Ok(());
    }

    run_repeat_runlist_pipeline(&rg_files, 0, opts.min_len, opts.fill_fragment, "out.json")?;

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
        ${pgr} psl to-rg lifted.psl -o coverage.rg
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
        // The rg line is `{contig}:{start}-{end}`; `to-rg` writes the real
        // contig name, which may contain '.' or ':' (e.g. `NC_000913.1`,
        // `chr1:alt`). Parse via the range suffix and rewrite with the
        // dot/colon-free placeholder so the downstream rg parser cannot
        // misread the name separators.
        let Some((contig, start, end)) = crate::libs::fmt::psl::parse_subrange(&line) else {
            log::warn!("skipping unparseable coverage .rg line: {}", line);
            continue;
        };
        let safe = name_map
            .get(&contig)
            .cloned()
            .unwrap_or_else(|| contig.clone());
        writer.write_fmt(format_args!("{}:{}-{}\n", safe, start, end))?;
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
/// marker written after all table files are copied; the marker is an exact
/// copy of the `.ktab` table, so a size mismatch (e.g. a truncated table from
/// an interrupted write) also marks the cache stale instead of letting FastK
/// silently read a corrupt table and report empty/partial repeats.
fn cache_is_fresh(lib: &str, cache_prefix: &str) -> bool {
    let Ok(lib_meta) = std::fs::metadata(lib) else {
        return false;
    };
    let mut ktab_len: Option<u64> = None;
    for suffix in [".ktab", ".complete"] {
        let Ok(cache_meta) = std::fs::metadata(format!("{}{}", cache_prefix, suffix)) else {
            return false;
        };
        if !cache_meta.is_file() {
            return false;
        }
        if suffix == ".ktab" {
            ktab_len = Some(cache_meta.len());
        } else if ktab_len != Some(cache_meta.len()) {
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

    #[test]
    fn truncated_cache_table_is_not_fresh() {
        // Regression: a truncated `.ktab` (e.g. interrupted write) used to
        // pass `cache_is_fresh` (existence + mtime only), so FastK silently
        // read the corrupt table and e-kmer reported empty repeats. The
        // `.complete` marker is an exact copy of `.ktab`, so a size mismatch
        // must mark the cache stale and trigger a rebuild.
        let dir = tempfile::tempdir().unwrap();
        let lib = dir.path().join("lib.fa");
        let prefix = dir.path().join("lib.fa.repeat.k17");
        std::fs::write(&lib, ">seq\nACGT\n").unwrap();
        let full = vec![b'x'; 4096];
        std::fs::write(format!("{}.ktab", prefix.display()), &full).unwrap();
        std::fs::write(format!("{}.complete", prefix.display()), &full).unwrap();
        assert!(
            cache_is_fresh(lib.to_str().unwrap(), prefix.to_str().unwrap()),
            "intact cache must be fresh"
        );

        std::fs::write(format!("{}.ktab", prefix.display()), &full[..100]).unwrap();
        assert!(
            !cache_is_fresh(lib.to_str().unwrap(), prefix.to_str().unwrap()),
            "truncated .ktab must not be fresh"
        );
    }

    #[test]
    fn sequence_less_fasta_is_detected() {
        // Regression: a sequence-less repeat library made FastK segfault
        // ("terminated by signal: 11"); the pipeline now rejects it up front.
        let dir = tempfile::tempdir().unwrap();
        let empty = dir.path().join("empty.fa");
        std::fs::write(&empty, ">only_header\n").unwrap();
        assert!(!has_sequences(empty.to_str().unwrap()).unwrap());

        let normal = dir.path().join("normal.fa");
        std::fs::write(&normal, ">chr\nACGT\n").unwrap();
        assert!(has_sequences(normal.to_str().unwrap()).unwrap());

        let no_records = dir.path().join("no_records.fa");
        std::fs::write(&no_records, "").unwrap();
        assert!(!has_sequences(no_records.to_str().unwrap()).unwrap());
    }
}
