//! LASTZ-based putative SD detection (migration design: notes/references/biser.md §6.8).

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
    let plain_target = decompress_if_gz(target_files, workdir)?;
    let plain_query = decompress_if_gz(query_files, workdir)?;

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
    for entry in std::fs::read_dir(workdir)? {
        let path = entry?.path();
        if path.extension().and_then(|e| e.to_str()) != Some("lav") {
            continue;
        }
        let reader = crate::libs::io::reader(&path.to_string_lossy())
            .with_context(|| format!("failed to open LAV {}", path.display()))?;
        let mut psl_bytes = Vec::new();
        lav_to_psl(reader, &mut psl_bytes, None, false)?;
        for line in std::io::Cursor::new(psl_bytes).lines() {
            let line = line?;
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

/// Decompress `.gz` FASTA files into the workdir; returns plain-text paths.
fn decompress_if_gz(
    files: Vec<std::path::PathBuf>,
    workdir: &str,
) -> anyhow::Result<Vec<std::path::PathBuf>> {
    let mut plain = Vec::with_capacity(files.len());
    for f in files {
        let is_gz = f.extension().and_then(|e| e.to_str()) == Some("gz");
        if !is_gz {
            plain.push(f);
            continue;
        }
        let base = crate::libs::io::get_basename(&f.to_string_lossy()).unwrap_or_default();
        let out = std::path::Path::new(workdir).join(format!("{base}.plain.fa"));
        let mut reader = crate::libs::io::reader(&f.to_string_lossy())
            .with_context(|| format!("failed to open {}", f.display()))?;
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&out)?);
        std::io::copy(&mut reader, &mut writer)?;
        plain.push(out);
    }
    Ok(plain)
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
}
