//! Native pgi self-alignment based putative SD detection.

use crate::libs::fmt::psl::{parse_or_warn, Psl};
use crate::libs::pl::PipelineCtx;
use anyhow::Context;
use cmd_lib::run_cmd;
use std::io::BufRead;
use std::path::Path;

/// Options for the native pgi SD search.
pub struct SearchPgiOptions {
    /// Worker threads for the pgi alignment.
    pub parallel: usize,
    /// Minimum alignment block length in bp (T2T-CHM13 SD standard: 1000).
    pub min_len: u32,
    /// Minimum block identity, 0.0-1.0 (T2T-CHM13 SD standard: 0.90).
    pub min_identity: f64,
}

impl Default for SearchPgiOptions {
    fn default() -> Self {
        Self {
            parallel: 4,
            min_len: 1000,
            min_identity: 0.90,
        }
    }
}

/// Run `pgr align pgi` for a target/query pair (or self-alignment), read the
/// PSL output, and keep hits passing the SD filters (`min_len`/`min_identity`).
///
/// `workdir` receives the intermediate `.psl` file; it must exist and be
/// writable. Chaining/refinement is intentionally NOT done here — the caller
/// passes the PSL through `pgr sd align` (mirror of `search_lastz`).
pub fn pgi_to_hits(
    target: &str,
    query: &str,
    is_self: bool,
    workdir: &str,
    opts: &SearchPgiOptions,
) -> anyhow::Result<Vec<Psl>> {
    anyhow::ensure!(
        !super::is_pgi_input(target) && !super::is_pgi_input(query),
        "sd search/cross needs genome FASTA (plain or .gz); a .pgi index \
         aligns without extension sequences and every block would score 0"
    );
    let ctx = PipelineCtx::new("pgr_sd_search_pgi_")?;
    let pgr = ctx.pgr.clone();
    let abs_target = ctx.abs_path(target)?;
    let abs_query = ctx.abs_path(query)?;
    let raw = ctx.abs_path(&Path::new(workdir).join("hits.raw.psl").to_string_lossy())?;
    let _cwd_guard = ctx.enter()?;
    let parallel = opts.parallel;
    if is_self {
        run_cmd!(${pgr} align pgi ${abs_target} -o ${raw} --parallel ${parallel})?;
    } else {
        run_cmd!(${pgr} align pgi ${abs_target} ${abs_query} -o ${raw} --parallel ${parallel})?;
    }

    let mut hits = Vec::new();
    let mut reader = crate::libs::io::reader(&raw)
        .with_context(|| format!("failed to open pgi SD hits {}", raw))?;
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if let Some(psl) = parse_or_warn(line.trim_end(), false)? {
            if super::passes_sd_filters(&psl, opts.min_len, opts.min_identity) {
                hits.push(psl);
            }
        }
    }
    Ok(hits)
}

/// Run `pgr align pgi` self-alignment on a genome and filter the PSL hits.
pub fn search_pgi(
    genome: &str,
    workdir: &str,
    opts: &SearchPgiOptions,
) -> anyhow::Result<Vec<Psl>> {
    pgi_to_hits(genome, genome, true, workdir, opts)
}
