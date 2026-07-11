use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::pbit::compressor::Compressor;
use std::path::PathBuf;

/// RAII guard that keeps a `tempfile::NamedTempFile` alive on drop unless
/// disarmed. The temp file is deleted automatically when the guard is dropped;
/// a successful in-place append disarms the guard before renaming the file.
struct TempFileGuard {
    file: Option<tempfile::NamedTempFile>,
}

impl TempFileGuard {
    fn new(file: tempfile::NamedTempFile) -> Self {
        Self { file: Some(file) }
    }

    /// Keep the temporary file so it can be renamed over the original archive.
    fn disarm(mut self) -> anyhow::Result<PathBuf> {
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

impl Drop for TempFileGuard {
    fn drop(&mut self) {
        // `NamedTempFile` deletes the underlying file on drop automatically.
        let _ = self.file.take();
    }
}

/// Build the clap subcommand for append.
pub fn make_subcommand() -> Command {
    Command::new("append")
        .about("Appends samples to an existing pbit archive")
        .after_help(
            r###"
This command appends new sample FASTA files to an existing pbit archive.
The reference is already embedded in the archive, so no -r is needed.

When `--paf` is provided, segments covered by PAF alignments are CIGAR-encoded
(replacing LZ-diff); uncovered segments fall back to LZ-diff.

Notes:
* Sample names are derived from input FASTA basenames (use --name to override)
* If -o is omitted, the input archive is modified in place
* If -o is specified, the input archive is copied to the output path first
* Reference and sample FASTA files may be plain text or gzipped (.gz)
* Contigs in sample FASTA that do not match any reference contig are skipped
* Only ACGTN characters are supported; IUPAC degenerate codes are mapped to N
* `--paf` files are paired with `-i` files by order; `--name` and `--paf`
  are mutually exclusive (use the TSV's optional 3rd column for PAF)

Examples:
1. Append a sample in place:
   pgr pbit append archive.pbit -i new_sample.fa

2. Append multiple samples to a new archive:
   pgr pbit append archive.pbit -i s1.fa -i s2.fa -o new_archive.pbit

3. Provide sample names via TSV:
   pgr pbit append archive.pbit --name samples.tsv -o new_archive.pbit

4. CIGAR-driven encoding with PAF:
   pgr pbit append archive.pbit -i sample.fa -p sample.paf -o new_archive.pbit
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Existing pbit archive to append to",
        ))
        .arg(crate::cmd_pgr::args::pbit_infiles_arg())
        .arg(crate::cmd_pgr::args::outfile_arg_optional())
        .arg(crate::cmd_pgr::args::pbit_name_arg())
        .arg(crate::cmd_pgr::args::pbit_paf_arg())
}

/// Execute the append command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let outfile_opt = args.get_one::<String>("outfile");

    let samples = super::collect_samples_from_args(args)?;

    // Determine the working path: copy if -o specified, else in-place via a
    // temp file + atomic rename so a mid-append failure cannot corrupt the
    // input archive.
    let in_place = outfile_opt.is_none();
    let mut temp_guard: Option<TempFileGuard> = None;
    let work_path = match outfile_opt {
        Some(out) => {
            // Refuse to overwrite the source archive in-place: fs::copy would
            // truncate the source before reading, destroying the archive.
            let in_path = std::path::Path::new(infile);
            let out_path = std::path::Path::new(out);
            if in_path == out_path {
                anyhow::bail!("outfile must differ from infile; omit -o for in-place append");
            }
            // Canonicalize when both paths exist; treat failure as "not the
            // same file" rather than aborting the operation.
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
                anyhow::bail!("outfile must differ from infile; omit -o for in-place append");
            }
            std::fs::copy(infile, out)
                .with_context(|| format!("failed to copy {} to {}", infile, out))?;
            out.clone()
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
                        "failed to create temp file for in-place append in {}",
                        parent.display()
                    )
                })?;
            let tmp_path = temp_file.path().to_path_buf();
            std::fs::copy(infile, &tmp_path)
                .with_context(|| "failed to stage temp file for in-place append")?;
            temp_guard = Some(TempFileGuard::new(temp_file));
            tmp_path.to_string_lossy().into_owned()
        }
    };

    let mut comp = Compressor::open_for_append(&work_path)
        .with_context(|| format!("failed to open pbit archive for append: {}", work_path))?;

    let mut cmd_line = format!("pgr pbit append {}", infile);
    if let Some(out) = outfile_opt {
        cmd_line.push_str(&format!(" -o {}", out));
    }
    comp.set_cmd_line(&cmd_line);

    for (name, path, paf_opt) in &samples {
        match paf_opt {
            Some(paf) => comp
                .append_sample_with_paf(name, path, paf)
                .with_context(|| format!("failed to append sample '{}' with PAF", name))?,
            None => comp
                .append_sample(name, path)
                .with_context(|| format!("failed to append sample '{}'", name))?,
        }
    }
    comp.finish().context("failed to finalize pbit archive")?;

    // Atomic in-place replacement: rename temp file over the input archive.
    if in_place {
        // Disarm the guard so the temp file survives the rename.
        let rename_from = if let Some(guard) = temp_guard.take() {
            guard
                .disarm()
                .with_context(|| "failed to prepare temp file for in-place rename")?
        } else {
            PathBuf::from(&work_path)
        };
        std::fs::rename(&rename_from, infile).with_context(|| {
            format!(
                "failed to finalize in-place append: rename {} -> {}",
                rename_from.display(),
                infile
            )
        })?;
    }

    Ok(())
}
