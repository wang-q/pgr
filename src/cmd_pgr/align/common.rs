//! Shared LASTZ pipeline helpers for the `align fill` / `align rest`
//! subcommands: 2bit conversion, PSL/sizes readers, subrange extraction,
//! LASTZ invocation and LAV → PSL → lift conversion.

use anyhow::Result;
use pgr::libs::fmt::psl::Psl;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// Path to a 2bit version of `path`: reuse a sibling `.2bit` (e.g. `ref.fa`
/// -> `ref.2bit`) when present, else convert the input with `pgr fa to-2bit`.
pub(super) fn to_2bit(pgr: &str, path: &str, out_2bit: &str) -> anyhow::Result<String> {
    if std::path::Path::new(path)
        .extension()
        .and_then(|e| e.to_str())
        == Some("2bit")
    {
        return Ok(path.to_string());
    }
    let sibling = sibling_2bit_path(std::path::Path::new(path));
    if sibling.is_file() {
        log::info!("reusing sibling 2bit {}", sibling.display());
        return Ok(sibling.to_string_lossy().into_owned());
    }
    let status = std::process::Command::new(pgr)
        .args(["fa", "to-2bit", path, "-o", out_2bit])
        .status()?;
    if !status.success() {
        anyhow::bail!("`pgr fa to-2bit {path}` failed with status {status}");
    }
    Ok(out_2bit.to_string())
}

/// Sibling 2bit path for a genome input: `ref.fa` -> `ref.2bit`,
/// `ref.fa.gz` -> `ref.fa.2bit` (same sibling convention as the `.pgi` index).
fn sibling_2bit_path(path: &std::path::Path) -> PathBuf {
    let mut p = path.to_path_buf();
    p.set_extension("2bit");
    p
}

/// Read a PSL file into records.
pub(super) fn read_psl(path: &str) -> anyhow::Result<Vec<Psl>> {
    let reader = pgr::libs::io::reader(path)?;
    let mut out = Vec::new();
    for p in pgr::libs::fmt::psl::iter_psl(reader) {
        out.push(p?);
    }
    Ok(out)
}

/// Read `name<TAB>len` contig sizes from `pgr 2bit size`.
pub(super) fn read_2bit_sizes(pgr: &str, two_bit: &str) -> Result<BTreeMap<String, i32>> {
    let out = std::process::Command::new(pgr)
        .args(["2bit", "size", two_bit, "-o", "stdout"])
        .output()?;
    if !out.status.success() {
        anyhow::bail!("`pgr 2bit size {two_bit}` failed");
    }
    let mut sizes = BTreeMap::new();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let mut it = line.split_whitespace();
        if let (Some(name), Some(len)) = (it.next(), it.next()) {
            if let Ok(len) = len.parse::<i32>() {
                sizes.insert(name.to_string(), len);
            }
        }
    }
    Ok(sizes)
}

/// Extract a subrange `[start, end)` (0-based half-open) from a 2bit file with
/// `pgr 2bit range`, writing a single-record FASTA (header `chr:start+1-end`).
pub(super) fn extract_2bit_range(
    pgr: &str,
    two_bit: &str,
    chr: &str,
    start: i32,
    end: i32,
    out: &Path,
) -> Result<()> {
    let range = format!("{chr}:{}-{}", start + 1, end);
    let status = std::process::Command::new(pgr)
        .args(["2bit", "range", two_bit, &range, "-o"])
        .arg(out)
        .status()?;
    if !status.success() {
        anyhow::bail!("`pgr 2bit range {two_bit} {range}` failed with status {status}");
    }
    if !out.exists() {
        anyhow::bail!("`pgr 2bit range` produced no output for {range}");
    }
    Ok(())
}

/// Invoke LASTZ on one target/query pair, writing the LAV output.
pub(super) fn run_lastz_pair(
    target_fa: &Path,
    query_fa: &Path,
    common_args: &[String],
    out_lav: &Path,
) -> Result<()> {
    let t_arg = format!("{}[nameparse=darkspace]", target_fa.display());
    let q_arg = format!("{}[nameparse=darkspace]", query_fa.display());
    let mut cmd = std::process::Command::new("lastz");
    cmd.arg(&t_arg).arg(&q_arg);
    for a in common_args {
        cmd.arg(a);
    }
    cmd.arg(format!("--output={}", out_lav.display()));
    let out = cmd.output()?;
    if !out.status.success() {
        let msg = String::from_utf8_lossy(&out.stderr).trim().to_string();
        anyhow::bail!(
            "lastz failed for {} vs {}: {}",
            target_fa.display(),
            query_fa.display(),
            msg
        );
    }
    Ok(())
}

/// Run LASTZ on every (target, query, LAV-output) job in parallel, then
/// convert each LAV to PSL and lift the subrange coordinates back to genomic
/// coordinates (same logic as `pgr psl lift`, reused in-process).
pub(super) fn run_lastz_jobs(
    jobs: &[(PathBuf, PathBuf, PathBuf)],
    common_args: &[String],
    sizes: &BTreeMap<String, i32>,
    parallel: usize,
) -> Result<Vec<Psl>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(parallel)
        .build()?;
    pool.install(|| {
        use rayon::prelude::*;
        jobs.par_iter()
            .try_for_each(|(t_fa, q_fa, lav)| run_lastz_pair(t_fa, q_fa, common_args, lav))
    })?;

    // Convert + lift every LAV in parallel.
    let converted: Vec<Result<Vec<Psl>>> = pool.install(|| {
        use rayon::prelude::*;
        jobs.par_iter()
            .map(|(_, _, lav_path)| convert_lav(lav_path, sizes))
            .collect()
    });
    let mut out = Vec::new();
    for c in converted {
        out.extend(c?);
    }
    Ok(out)
}

/// Convert one LAV file to PSL records and lift the subrange coordinates.
fn convert_lav(lav_path: &Path, sizes: &BTreeMap<String, i32>) -> Result<Vec<Psl>> {
    let mut out = Vec::new();
    if !lav_path.exists() {
        return Ok(out);
    }
    let mut buf = Vec::new();
    let reader = pgr::libs::io::reader(lav_path.to_str().unwrap_or(""))?;
    pgr::libs::lav::lav_to_psl(reader, &mut buf, None, false)?;
    for p in pgr::libs::fmt::psl::iter_psl(std::io::Cursor::new(buf)) {
        let mut p = p?;
        if !p.lift_query(sizes) {
            log::warn!("failed to lift query: {}", p.q_name);
        }
        if !p.lift_target(sizes) {
            log::warn!("failed to lift target: {}", p.t_name);
        }
        out.push(p);
    }
    Ok(out)
}

/// Build the common LASTZ arguments shared by fill/rest (preset + query depth
/// + user overrides).
pub(super) fn build_common_args(
    preset: Option<&str>,
    query_depth: usize,
    lastz_args: Option<&str>,
) -> Result<(Vec<String>, Option<tempfile::NamedTempFile>)> {
    let (mut args, matrix_handle) = pgr::libs::lastz::build_common_args(preset, query_depth)?;
    if let Some(extra) = lastz_args {
        for arg in extra.split_whitespace() {
            args.push(arg.to_string());
        }
    }
    Ok((args, matrix_handle))
}
