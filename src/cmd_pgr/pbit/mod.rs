pub mod append;
pub mod append_ref;
pub mod create;
pub mod range;
pub mod some;
pub mod stat;
pub mod to_fa;
pub mod to_paf;

use anyhow::{Context, Result};
use clap::{ArgMatches, Command};
use std::path::PathBuf;

/// One sample specification: (name, fasta path, optional paf path, optional
/// reference spec).
type SampleSpec = (String, String, Option<String>, Option<String>);

/// RAII guard that keeps a `tempfile::NamedTempFile` alive on drop unless
/// disarmed. The temp file is deleted automatically when the guard is
/// dropped; a successful in-place update disarms the guard before renaming.
pub(crate) struct TempFileGuard {
    file: Option<tempfile::NamedTempFile>,
}

impl TempFileGuard {
    pub(crate) fn new(file: tempfile::NamedTempFile) -> Self {
        Self { file: Some(file) }
    }

    /// Keep the temporary file so it can be renamed over the original archive.
    pub(crate) fn disarm(mut self) -> Result<PathBuf> {
        let file = self
            .file
            .take()
            .expect("disarm called on an empty TempFileGuard");
        let (_, path) = file
            .keep()
            .context("failed to keep temp file for in-place rename")?;
        Ok(path)
    }
}

/// Prepare a work path for an in-place archive update: copy to a sibling temp
/// file (in-place) or to `-o` (must differ from the input), returning the
/// path to operate on plus the disarmable guard for in-place mode.
pub(crate) fn stage_work_path(
    infile: &str,
    outfile_opt: Option<&String>,
) -> Result<(String, Option<TempFileGuard>)> {
    match outfile_opt {
        Some(out) => {
            let in_path = std::path::Path::new(infile);
            let out_path = std::path::Path::new(out);
            if in_path == out_path {
                anyhow::bail!("outfile must differ from infile; omit -o for in-place update");
            }
            let same_file = if in_path.exists() && out_path.exists() {
                match (
                    std::fs::canonicalize(in_path),
                    std::fs::canonicalize(out_path),
                ) {
                    (Ok(i), Ok(o)) => i == o,
                    _ => false,
                }
            } else {
                false
            };
            if same_file {
                anyhow::bail!("outfile must differ from infile; omit -o for in-place update");
            }
            std::fs::copy(infile, out)
                .with_context(|| format!("failed to copy {} to {}", infile, out))?;
            Ok((out.clone(), None))
        }
        None => {
            let in_path = std::path::Path::new(infile);
            let parent = in_path
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or_else(|| std::path::Path::new("."));
            let temp_file = tempfile::Builder::new()
                .suffix(".pbit.tmp")
                .tempfile_in(parent)
                .with_context(|| {
                    format!(
                        "failed to create temp file for in-place update in {}",
                        parent.display()
                    )
                })?;
            let tmp_path = temp_file.path().to_path_buf();
            std::fs::copy(infile, &tmp_path)
                .with_context(|| "failed to stage temp file for in-place update")?;
            let guard = TempFileGuard::new(temp_file);
            Ok((tmp_path.to_string_lossy().into_owned(), Some(guard)))
        }
    }
}

/// Read a TSV of `sample_name<TAB>fasta_path[<TAB>paf_path][<TAB>ref_name]`
/// lines (3rd/4th columns optional; the 4th selects the reference genome).
pub(crate) fn read_name_tsv(path: &str) -> Result<Vec<SampleSpec>> {
    let lines = pgr::libs::io::read_lines(path)
        .with_context(|| format!("failed to read name TSV: {}", path))?;
    let mut out = Vec::new();
    for (line_no, line) in lines.iter().enumerate() {
        let line_no = line_no + 1;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        let name = parts[0].trim().to_string();
        let fasta_path = parts
            .get(1)
            .ok_or_else(|| anyhow::anyhow!("line {}: missing FASTA path: {}", line_no, trimmed))?
            .trim()
            .to_string();
        let paf_path = parts
            .get(2)
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty());
        let ref_spec = parts
            .get(3)
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty());
        if name.is_empty() || fasta_path.is_empty() {
            anyhow::bail!("line {}: empty name or path: {}", line_no, trimmed);
        }
        out.push((name, fasta_path, paf_path, ref_spec));
    }
    Ok(out)
}

/// Collect `(sample_name, fasta_path, paf_path_opt, ref_spec_opt)` from
/// `--name` TSV or `-i`/`--paf`.
pub(crate) fn collect_samples_from_args(args: &ArgMatches) -> Result<Vec<SampleSpec>> {
    let has_name = args.get_one::<String>("name").is_some();
    let has_infiles = args.get_many::<String>("infiles").is_some();
    let has_paf = args.get_many::<String>("paf").is_some();
    if has_name && has_paf {
        anyhow::bail!(
            "--name and --paf are mutually exclusive (use --name TSV with 3rd column for PAF)"
        );
    }
    if has_name && has_infiles {
        anyhow::bail!("--name and -i/--infile are mutually exclusive");
    }

    let samples: Vec<SampleSpec> = if let Some(name_tsv) = args.get_one::<String>("name") {
        read_name_tsv(name_tsv)?
    } else {
        let infiles = args
            .get_many::<String>("infiles")
            .ok_or_else(|| anyhow::anyhow!("no input files: provide -i or --name"))?;
        let pafs: Vec<String> = args
            .get_many::<String>("paf")
            .map(|v| v.cloned().collect())
            .unwrap_or_default();
        if !pafs.is_empty() && pafs.len() != infiles.len() {
            anyhow::bail!(
                "--paf count ({}) does not match -i count ({})",
                pafs.len(),
                infiles.len()
            );
        }
        let mut pairs = Vec::new();
        for (i, path) in infiles.enumerate() {
            let name = pgr::libs::io::basename_or_err(path)?;
            let paf = pafs.get(i).cloned();
            pairs.push((name, path.clone(), paf, None));
        }
        pairs
    };

    if samples.is_empty() {
        anyhow::bail!("no sample FASTA files provided");
    }

    // 2026-08-09 (2026-08-09): PAF is mandatory — every sample must come with a
    // PAF (via `--paf` or the TSV 3rd column); the no-PAF compression path
    // is retired. An empty PAF file is allowed (all segments fall back to
    // LZ-diff/Raw).
    for (name, _, paf, _) in &samples {
        if paf.is_none() {
            anyhow::bail!(
                "sample '{}' has no PAF: --paf is required (or the 3rd column \
                 of --name TSV); pass an empty PAF file to skip CIGAR encoding",
                name
            );
        }
    }

    // Reject duplicate sample names within this command. Sample names derived
    // from `-i` basenames collapse distinct files (e.g. `sample.1.fa` and
    // `sample.2.fa` both become `sample`), and `append_sample` would silently
    // merge their segments into one sample, corrupting the archive on extract.
    let mut seen = std::collections::HashSet::new();
    for (name, _, _, _) in &samples {
        if !seen.insert(name.as_str()) {
            anyhow::bail!(
                "duplicate sample name '{}'; sample names must be distinct \
                 (use --name to assign explicit names)",
                name
            );
        }
    }

    Ok(samples)
}

/// Resolve a user reference spec (index or name) against the `-r` list.
pub(crate) fn resolve_ref_id(spec: Option<&str>, ref_fastas: &[&str]) -> Result<u32> {
    let Some(spec) = spec else {
        if ref_fastas.len() > 1 {
            log::warn!(
                "no reference specified for sample; defaulting to reference 0 of {} \
                 (route samples with the 4th TSV column or --name)",
                ref_fastas.len()
            );
        }
        return Ok(0);
    };
    if let Ok(n) = spec.parse::<usize>() {
        if n < ref_fastas.len() {
            return Ok(n as u32);
        }
        anyhow::bail!(
            "reference index {} out of range ({} references)",
            n,
            ref_fastas.len()
        );
    }
    let names: Vec<String> = ref_fastas
        .iter()
        .map(|p| pgr::libs::io::get_basename(p).unwrap_or_else(|| p.to_string()))
        .collect();
    names
        .iter()
        .position(|n| n == spec)
        .map(|i| i as u32)
        .ok_or_else(|| anyhow::anyhow!("reference '{}' not found (available: {names:?})", spec))
}

/// Build the clap subcommand for pbit.
pub fn make_subcommand() -> Command {
    Command::new("pbit")
        .about("Manages pbit (population 2bit + delta) files")
        .after_help(
            r###"Subcommand groups:

* build:     create / append / append-ref
* info:      stat
* subset:    range / some
* transform: to-fa

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(create::make_subcommand())
        .subcommand(append::make_subcommand())
        .subcommand(append_ref::make_subcommand())
        .subcommand(range::make_subcommand())
        .subcommand(some::make_subcommand())
        .subcommand(stat::make_subcommand())
        .subcommand(to_fa::make_subcommand())
        .subcommand(to_paf::make_subcommand())
}

/// Execute the pbit command.
pub fn execute(args: &ArgMatches) -> Result<()> {
    match args.subcommand() {
        Some(("create", sub_matches)) => create::execute(sub_matches),
        Some(("append", sub_matches)) => append::execute(sub_matches),
        Some(("append-ref", sub_matches)) => append_ref::execute(sub_matches),
        Some(("range", sub_matches)) => range::execute(sub_matches),
        Some(("some", sub_matches)) => some::execute(sub_matches),
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        Some(("to-fa", sub_matches)) => to_fa::execute(sub_matches),
        Some(("to-paf", sub_matches)) => to_paf::execute(sub_matches),
        _ => Ok(()),
    }
}
