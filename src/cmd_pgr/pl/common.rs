//! Shared helpers for `pgr pl` pipeline subcommands.

use cmd_lib::*;
use std::io::BufRead;

/// Read chromosome names from a `chr.sizes` file (lines of `<chr>\t<size>`).
///
/// Returns the first column in file order. Supports stdin and `.gz` via
/// [`pgr::libs::io::read_names_as_vec`].
pub fn read_chr_names(sizes_file: &str) -> anyhow::Result<Vec<String>> {
    pgr::libs::io::read_names_as_vec(sizes_file)
}

/// Resolve `path` to an absolute path string. `stdout` is passed through as-is.
pub fn abs_path_or_stdout(path: &str) -> anyhow::Result<String> {
    if path == "stdout" {
        Ok(path.to_string())
    } else {
        Ok(intspan::absolute_path(path)?.display().to_string())
    }
}

/// Run `Profex -z genome` per chromosome and write `.rg` files.
///
/// For each chromosome, runs `Profex -z genome <sn>` writing `prof.<sn>.txt`,
/// then scans lines with `re_prof` capturing `start` and `end` (1-based inclusive
/// in output). If `min_depth` is set and the regex has a `depth` capture group,
/// entries with depth below the threshold are skipped. Returns the list of
/// `prof.<sn>.rg` file names.
#[allow(unused_variables)]
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

        let reader = pgr::reader(&format!("prof.{}.txt", sn))?;

        let rg_file = format!("prof.{}.rg", sn);
        let mut writer = pgr::writer(&rg_file)?;

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

/// Shared pipeline context: current dir, pgr executable, and tempdir.
///
/// Created at the start of a pipeline; call [`PipelineCtx::enter`] to switch
/// into the tempdir and [`PipelineCtx::leave`] to restore the original cwd.
pub struct PipelineCtx {
    /// Original working directory, restored by `leave()`.
    pub curdir: std::path::PathBuf,
    /// Absolute path to the current `pgr` executable.
    pub pgr: String,
    /// Owned tempdir; dropped when the ctx is dropped.
    pub tempdir: tempfile::TempDir,
}

impl PipelineCtx {
    /// Create a new context with a tempdir using `prefix` (e.g. `"pgr_rm_"`).
    ///
    /// Prints the `==> Paths` info block.
    pub fn new(prefix: &str) -> anyhow::Result<Self> {
        let curdir = std::env::current_dir()?;
        let pgr = std::env::current_exe()?.display().to_string();
        let tempdir = tempfile::Builder::new().prefix(prefix).tempdir()?;
        let tempdir_str = tempdir.path().to_str().unwrap();

        run_cmd!(info "==> Paths")?;
        run_cmd!(info "    \"pgr\"     = ${pgr}")?;
        run_cmd!(info "    \"curdir\"  = ${curdir:?}")?;
        run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

        Ok(Self {
            curdir,
            pgr,
            tempdir,
        })
    }

    /// Resolve `p` to an absolute path string.
    pub fn abs_path(&self, p: &str) -> anyhow::Result<String> {
        Ok(intspan::absolute_path(p)?.display().to_string())
    }

    /// Switch the current working directory into the tempdir.
    pub fn enter(&self) -> anyhow::Result<()> {
        let tempdir_str = self.tempdir.path().to_str().unwrap();
        run_cmd!(info "==> Switch to tempdir")?;
        std::env::set_current_dir(tempdir_str)?;
        Ok(())
    }

    /// Restore the original working directory.
    pub fn leave(&self) -> anyhow::Result<()> {
        std::env::set_current_dir(&self.curdir)?;
        Ok(())
    }
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
        run_cmd!(info "==> FastK on repeat")?;
        run_cmd!(
            FastK -t -k${opt_kmer} -Nrepeat ${abs_repeat}
        )?;
        run_cmd!(info "==> FastK on genome")?;
        run_cmd!(
            FastK -p:repeat -k${opt_kmer} -Ngenome ${abs_infile}
        )?;
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
    let chrs = read_chr_names("chr.sizes")?;

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
