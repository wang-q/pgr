//! LASTZ-based putative SD detection (migration design: notes/design/sd.md).

use crate::libs::fmt::lav::lav_to_psl;
use crate::libs::fmt::psl::{parse_or_warn, Psl};
use crate::libs::lastz::{build_common_args, run_lastz, RunLastzOptions};
use anyhow::Context;
use std::io::BufRead;

/// Options for LASTZ-based SD search.
pub struct SearchLastzOptions {
    /// lastz preset (set01..set07); defaults to set01 when `None`.
    pub preset: Option<String>,
    /// lastz query-depth threshold (`--querydepth=keep,nowarn:N`).
    pub query_depth: usize,
    /// Minimum alignment block length in bp (T2T-CHM13 SD standard: 1000).
    pub min_len: u32,
    /// Minimum block identity, 0.0-1.0 (T2T-CHM13 SD standard: 0.90).
    pub min_identity: f64,
    /// Worker threads for the lastz batch.
    pub parallel: usize,
}

impl Default for SearchLastzOptions {
    fn default() -> Self {
        Self {
            preset: Some("set01".to_string()),
            query_depth: 50,
            min_len: 1000,
            min_identity: 0.90,
            parallel: 4,
        }
    }
}

/// Run lastz for a target/query pair, convert LAV to PSL, and filter hits by
/// `min_len` / `min_identity`. Returns the surviving PSL records.
///
/// `workdir` receives the intermediate `.lav` files; it must exist and be
/// writable. Chaining/refinement is intentionally NOT done here — the caller
/// passes the PSL through UCSC chain/net (see the migration design).
pub fn lastz_to_hits(
    target: &str,
    query: &str,
    is_self: bool,
    workdir: &str,
    opts: &SearchLastzOptions,
) -> anyhow::Result<Vec<Psl>> {
    anyhow::ensure!(
        !super::is_pgi_input(target) && !super::is_pgi_input(query),
        "sd search/cross needs genome FASTA (plain or .gz), not a .pgi index"
    );
    if which::which("lastz").is_err() {
        anyhow::bail!("lastz not found in PATH. Please install lastz first.");
    }

    let mut target_files = crate::libs::fmt::fa::find_fasta_files(target);
    target_files.sort();
    if target_files.is_empty() {
        anyhow::bail!("no FASTA files found in {target}");
    }
    let mut query_files = crate::libs::fmt::fa::find_fasta_files(query);
    query_files.sort();
    if query_files.is_empty() {
        anyhow::bail!("no FASTA files found in {query}");
    }

    // lastz cannot read gzipped input; decompress .gz files into the workdir.
    let (plain_target, plain_query) = decompress_target_query(&target_files, &query_files, workdir)?;

    let (common_args, _matrix_handle) =
        build_common_args(opts.preset.as_deref(), opts.query_depth)?;
    let run_opts = RunLastzOptions {
        depth: opts.query_depth,
        is_self,
        common_args,
        output_dir: workdir.to_string(),
        parallel: opts.parallel,
    };
    run_lastz(plain_target, plain_query, run_opts)?;

    // Convert every LAV in the workdir to PSL and apply the SD filters.
    let mut hits = Vec::new();
    // Sort the LAV files so the output PSL order is deterministic across
    // runs (read_dir order is filesystem-dependent).
    let mut lav_files: Vec<std::path::PathBuf> = std::fs::read_dir(workdir)?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().and_then(|e| e.to_str()) == Some("lav"))
        .collect();
    lav_files.sort();
    for path in lav_files {
        let reader = crate::libs::io::reader(&path.to_string_lossy())
            .with_context(|| format!("failed to open LAV {}", path.display()))?;
        let mut psl_bytes = Vec::new();
        lav_to_psl(reader, &mut psl_bytes, None, false)?;
        for line in std::io::Cursor::new(psl_bytes).lines() {
            let line = line?;
            if line.trim_start().starts_with('#') {
                continue;
            }
            if let Some(psl) = parse_or_warn(&line, false)? {
                if super::passes_sd_filters(&psl, opts.min_len, opts.min_identity) {
                    hits.push(psl);
                }
            }
        }
    }
    Ok(hits)
}

/// Run `lastz --self` on a genome and filter the resulting PSL hits.
pub fn search_lastz(
    genome: &str,
    workdir: &str,
    opts: &SearchLastzOptions,
) -> anyhow::Result<Vec<Psl>> {
    lastz_to_hits(genome, genome, true, workdir, opts)
}

/// Decompress any `.gz` FASTA inputs into `workdir`, returning plain-text
/// paths for the target and query lists.
///
/// A single used-basename set spans both lists so a cross-mode run cannot have
/// the target and query overwrite each other's decompressed file when they
/// share a basename (`a/sample.fa.gz` vs `b/sample.fa.gz` would both map to
/// `sample.plain.fa` otherwise, and the query would silently replace the
/// target). Each unique input path is decompressed once, so self-mode
/// (target == query lists) still maps every file to the same plain path on
/// both sides and `run_lastz`'s `--self` detection keeps working.
fn decompress_target_query(
    target_files: &[std::path::PathBuf],
    query_files: &[std::path::PathBuf],
    workdir: &str,
) -> anyhow::Result<(Vec<std::path::PathBuf>, Vec<std::path::PathBuf>)> {
    let mut plain: std::collections::HashMap<std::path::PathBuf, std::path::PathBuf> =
        std::collections::HashMap::new();
    let mut used: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut assign = |f: &std::path::PathBuf| -> anyhow::Result<std::path::PathBuf> {
        if let Some(p) = plain.get(f) {
            return Ok(p.clone());
        }
        let out = if f.extension().and_then(|e| e.to_str()) == Some("gz") {
            let base = crate::libs::io::get_basename(&f.to_string_lossy()).unwrap_or_default();
            let out = if used.insert(base.clone()) {
                std::path::Path::new(workdir).join(format!("{base}.plain.fa"))
            } else {
                std::path::Path::new(workdir).join(format!("{base}.{}.plain.fa", plain.len()))
            };
            let mut reader = crate::libs::io::reader(&f.to_string_lossy())
                .with_context(|| format!("failed to open {}", f.display()))?;
            let mut writer = std::io::BufWriter::new(std::fs::File::create(&out)?);
            std::io::copy(&mut reader, &mut writer)?;
            out
        } else {
            f.clone()
        };
        plain.insert(f.clone(), out.clone());
        Ok(out)
    };
    let mut plain_target = Vec::with_capacity(target_files.len());
    for f in target_files {
        plain_target.push(assign(f)?);
    }
    let mut plain_query = Vec::with_capacity(query_files.len());
    for f in query_files {
        plain_query.push(assign(f)?);
    }
    Ok((plain_target, plain_query))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::sd::{passes_sd_filters, psl_block_len, psl_identity};

    fn psl(m: u32, mm: u32, ins: u32) -> Psl {
        Psl {
            match_count: m,
            mismatch_count: mm,
            q_base_insert: ins as i32,
            ..Default::default()
        }
    }

    #[test]
    fn block_len_and_identity() {
        let p = psl(900, 50, 50);
        assert_eq!(psl_block_len(&p), 1000);
        assert!((psl_identity(&p) - 0.90).abs() < 1e-9);
    }

    #[test]
    fn identity_excludes_inserts() {
        let p = psl(80, 0, 20);
        assert_eq!(psl_block_len(&p), 100);
        assert!((psl_identity(&p) - 0.80).abs() < 1e-9);
    }

    #[test]
    fn sd_filters_respect_length_and_identity() {
        let good = psl(900, 50, 50); // len 1000, id 0.90
        assert!(passes_sd_filters(&good, 1000, 0.90));
        let too_short = psl(800, 100, 50); // len 950
        assert!(!passes_sd_filters(&too_short, 1000, 0.90));
        let low_id = psl(800, 150, 50); // len 1000, id 0.80
        assert!(!passes_sd_filters(&low_id, 1000, 0.90));
    }

    #[test]
    fn decompress_colliding_basenames_stay_distinct() {
        // Regression: two `.fa.gz` files with the same basename in different
        // directories decompressed to the same flat output path, silently
        // overwriting the first. Each input must keep its own plain file.
        let dir = tempfile::TempDir::new().unwrap();
        let workdir = dir.path().join("out");
        std::fs::create_dir_all(&workdir).unwrap();
        for (sub, name, seq) in [
            ("a", "dup", "ACGTACGTACGTACGT"),
            ("b", "dup", "TTTTCCCCAAAAGGGG"),
        ] {
            let subdir = dir.path().join(sub);
            std::fs::create_dir_all(&subdir).unwrap();
            let mut gz = flate2::write::GzEncoder::new(
                std::fs::File::create(subdir.join(format!("{name}.fa.gz"))).unwrap(),
                flate2::Compression::default(),
            );
            use std::io::Write;
            write!(gz, ">{name}\n{seq}\n").unwrap();
            gz.finish().unwrap();
        }
        let files = crate::libs::fmt::fa::find_fasta_files(dir.path());
        assert_eq!(files.len(), 2);
        let is_a = |p: &std::path::PathBuf| p.to_str().unwrap().contains("/a/");
        let a = files.iter().find(|p| is_a(p)).unwrap().clone();
        let b = files.iter().find(|p| !is_a(p)).unwrap().clone();
        // Cross mode: target and query are different files sharing a basename.
        let (target, query) = decompress_target_query(
            std::slice::from_ref(&a),
            std::slice::from_ref(&b),
            workdir.to_str().unwrap(),
        )
        .unwrap();
        assert_eq!(target.len(), 1);
        assert_eq!(query.len(), 1);
        assert_ne!(
            target[0], query[0],
            "cross-mode target/query must stay distinct"
        );
        assert_ne!(
            std::fs::read_to_string(&target[0]).unwrap(),
            std::fs::read_to_string(&query[0]).unwrap(),
            "each input must keep its own sequence"
        );
        // Self mode: the same file on both sides maps to the same plain path.
        let (self_t, self_q) = decompress_target_query(
            std::slice::from_ref(&a),
            std::slice::from_ref(&a),
            workdir.to_str().unwrap(),
        )
        .unwrap();
        assert_eq!(self_t, self_q, "self-mode target and query must match");
    }
}
