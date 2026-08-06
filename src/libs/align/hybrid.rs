//! Hybrid genome alignment: anchors + LASTZ gap filling.
//!
//! Mirrors the FastGA-gapfill scheme from the FastGA paper: use `align pgi`
//! for fast coarse anchors, then run LASTZ on the colinear gaps between
//! consecutive anchors to recover homology the anchors missed. Only syntenic
//! (colinear, same-strand) search is supported; the anchors and the LASTZ
//! records are written together and the downstream chainnet handles the merge.

use crate::libs::fmt::psl::Psl;
use crate::libs::lastz;
use anyhow::Result;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// A bounding box between two colinear anchors, in 0-based half-open
/// forward-strand coordinates on both genomes.
#[derive(Debug, Clone)]
pub struct Box {
    pub t_name: String,
    pub t_start: i32,
    pub t_end: i32,
    pub q_name: String,
    pub q_start: i32,
    pub q_end: i32,
    pub strand: String,
}

/// Options controlling the hybrid gap-fill run.
#[derive(Debug, Clone)]
pub struct HybridOptions {
    /// Box overlaps the anchor by this many bp on each side (LASTZ seeding buffer).
    pub overlap: i32,
    /// Gaps shorter than this (bp) are left to pgi.
    pub min_gap: i32,
    /// Gaps longer than this (bp) are skipped (likely novel sequence).
    pub max_gap: i32,
    /// LASTZ preset (set01..set07); None uses lastz defaults.
    pub preset: Option<String>,
    /// Query depth threshold passed to lastz.
    pub query_depth: usize,
    /// Extra LASTZ arguments appended after the preset.
    pub lastz_args: Option<String>,
    /// Number of parallel lastz jobs.
    pub parallel: usize,
}

/// Compute the fill boxes between consecutive colinear same-strand anchors.
///
/// Anchors are the PSL records (any aligner: pgi, FastGA, minimap2...),
/// grouped by (target, query, strand) and
/// sorted by target start. A box is emitted for each pair of adjacent,
/// non-overlapping, colinear anchors whose target *and* query gaps both fall
/// in `[min_gap, max_gap]`; the box is expanded beyond the gap by `overlap`
/// bp on each side so LASTZ can seed across the true homology boundary.
pub fn compute_boxes(psls: &[Psl], overlap: i32, min_gap: i32, max_gap: i32) -> Vec<Box> {
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
                strand: strand.clone(),
            });
        }
    }
    boxes
}

/// Run the hybrid gap-fill pipeline over two genomes.
///
/// `target_2bit`/`query_2bit` are 2bit files (converted from the inputs by the
/// caller) used for random-access box extraction (`pgr 2bit range`) and contig
/// sizes (`pgr 2bit size`). `avail_psl` are the anchor records. Returns the
/// anchors plus the LASTZ records, in that order (no dedup).
pub fn run_hybrid(
    pgr: &str,
    target_2bit: &str,
    query_2bit: &str,
    avail_psl: &[Psl],
    opts: &HybridOptions,
    tempdir: &Path,
) -> Result<Vec<Psl>> {
    if which::which("lastz").is_err() {
        anyhow::bail!("lastz not found in PATH. Please install lastz first.");
    }

    let boxes = compute_boxes(avail_psl, opts.overlap, opts.min_gap, opts.max_gap);
    if boxes.is_empty() {
        return Ok(avail_psl.to_vec());
    }

    // Contig sizes, for clamping and coordinate lifting.
    let mut sizes = read_2bit_sizes(pgr, target_2bit)?;
    sizes.extend(read_2bit_sizes(pgr, query_2bit)?);

    // Clamp box coordinates to contig bounds so the extracted subrange is
    // valid; otherwise `lift_query`/`lift_target` would reject the record
    // (subrange end > real size) and drop the alignment.
    let mut boxes = boxes;
    for b in &mut boxes {
        let t_size = sizes.get(&b.t_name).copied().unwrap_or(0);
        let q_size = sizes.get(&b.q_name).copied().unwrap_or(0);
        b.t_start = b.t_start.clamp(0, t_size);
        b.t_end = b.t_end.clamp(b.t_start, t_size);
        b.q_start = b.q_start.clamp(0, q_size);
        b.q_end = b.q_end.clamp(b.q_start, q_size);
        if b.t_end <= b.t_start || b.q_end <= b.q_start {
            log::warn!(
                "empty box after clamping on {} vs {}; skipping",
                b.t_name,
                b.q_name
            );
        }
    }
    boxes.retain(|b| b.t_end > b.t_start && b.q_end > b.q_start);

    // Build LASTZ common args (preset + query depth + user overrides).
    let (mut common_args, _matrix_handle) =
        lastz::build_common_args(opts.preset.as_deref(), opts.query_depth)?;
    if let Some(extra) = &opts.lastz_args {
        for arg in extra.split_whitespace() {
            common_args.push(arg.to_string());
        }
    }

    // Extract each box subrange with `pgr 2bit range` and write one
    // single-sequence FASTA per box for the target and query sides.
    std::fs::create_dir_all(tempdir)?;
    let mut jobs: Vec<(PathBuf, PathBuf, PathBuf)> = Vec::with_capacity(boxes.len());
    for (i, b) in boxes.iter().enumerate() {
        let t_fa = tempdir.join(format!("t_{i}.fa"));
        let q_fa = tempdir.join(format!("q_{i}.fa"));
        let lav = tempdir.join(format!("box_{i}.lav"));
        extract_2bit_range(pgr, target_2bit, &b.t_name, b.t_start, b.t_end, &t_fa)?;
        extract_2bit_range(pgr, query_2bit, &b.q_name, b.q_start, b.q_end, &q_fa)?;
        jobs.push((t_fa, q_fa, lav));
    }

    // Run LASTZ per box in parallel.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.parallel)
        .build()?;
    pool.install(|| {
        use rayon::prelude::*;
        jobs.par_iter()
            .try_for_each(|(t_fa, q_fa, lav)| run_lastz_pair(t_fa, q_fa, &common_args, lav))
    })?;

    // Convert each LAV to PSL and lift the subrange coordinates back to
    // genomic coordinates.
    let mut lastz_psl = Vec::new();
    for (_, _, lav_path) in &jobs {
        if !lav_path.exists() {
            continue;
        }
        let mut buf = Vec::new();
        let reader = crate::reader(lav_path.to_str().unwrap_or(""))?;
        crate::libs::lav::lav_to_psl(reader, &mut buf, None, false)?;
        for p in crate::libs::fmt::psl::iter_psl(std::io::Cursor::new(buf)) {
            let mut p = p?;
            let _ = p.lift_query(&sizes);
            let _ = p.lift_target(&sizes);
            lastz_psl.push(p);
        }
    }

    // Emit the anchors and the LASTZ records together; the downstream
    // chainnet pipeline handles the overlap/merge.
    let mut out = avail_psl.to_vec();
    out.extend(lastz_psl);
    Ok(out)
}

/// Read `name<TAB>len` contig sizes from `pgr 2bit size`.
fn read_2bit_sizes(pgr: &str, two_bit: &str) -> Result<BTreeMap<String, i32>> {
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
fn extract_2bit_range(
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
fn run_lastz_pair(
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

/// A minimal PSL helper for tests (target-starts/ends and strand only).
#[cfg(test)]
fn tpsl(t_name: &str, t_start: i32, t_end: i32, q_name: &str, strand: &str) -> Psl {
    Psl {
        t_name: t_name.into(),
        t_start,
        t_end,
        q_name: q_name.into(),
        q_start: t_start,
        q_end: t_end,
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
            tpsl("t", 0, 100, "q", "+"),
            tpsl("t", 500, 600, "q", "+"),
            tpsl("t", 1000, 1100, "q", "+"),
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
        let anchors = vec![tpsl("t", 0, 500, "q", "+"), tpsl("t", 400, 900, "q", "+")];
        assert!(compute_boxes(&anchors, 10, 50, 1000).is_empty());
        // Gap below min_gap: no box.
        let anchors = vec![tpsl("t", 0, 100, "q", "+"), tpsl("t", 120, 220, "q", "+")];
        assert!(compute_boxes(&anchors, 10, 50, 1000).is_empty());
        // Gap above max_gap: no box.
        let anchors = vec![tpsl("t", 0, 100, "q", "+"), tpsl("t", 5000, 5100, "q", "+")];
        assert!(compute_boxes(&anchors, 10, 50, 1000).is_empty());
    }

    #[test]
    fn compute_boxes_minus_strand_query_gap_is_reversed() {
        // On the '-' strand the query order is reversed: the earlier target
        // anchor has the later query span.
        let anchors = [tpsl("t", 0, 100, "q", "-"), tpsl("t", 500, 600, "q", "-")];
        // For '-' we need matching plus-strand q spans: anchor A q[900,1000],
        // anchor B q[300,400] (B is earlier in query).
        let mut a = anchors[0].clone();
        a.q_start = 900;
        a.q_end = 1000;
        let mut b = anchors[1].clone();
        b.q_start = 300;
        b.q_end = 400;
        let boxes = compute_boxes(&[a, b], 10, 50, 1000);
        assert_eq!(boxes.len(), 1);
        // Query gap is [400, 900] (forward strand), box q -> [390, 910].
        assert_eq!(boxes[0].q_start, 390);
        assert_eq!(boxes[0].q_end, 910);
    }
}
