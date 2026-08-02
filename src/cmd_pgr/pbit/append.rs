//! Append samples to an existing pbit archive.

use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::pbit::compressor::Compressor;
use std::path::PathBuf;

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

    let in_place = outfile_opt.is_none();
    let (work_path, mut temp_guard) = super::stage_work_path(infile, outfile_opt)?;

    let mut comp = Compressor::open_for_append(&work_path)
        .with_context(|| format!("failed to open pbit archive for append: {}", work_path))?;

    let mut cmd_line = format!("pgr pbit append {}", infile);
    if let Some(out) = outfile_opt {
        cmd_line.push_str(&format!(" -o {}", out));
    }
    comp.set_cmd_line(&cmd_line);

    for (name, path, paf_opt, ref_spec) in &samples {
        let ref_id = match ref_spec {
            Some(spec) => {
                let names = comp.ref_names();
                if let Ok(n) = spec.parse::<usize>() {
                    anyhow::ensure!(
                        n < names.len(),
                        "reference index {} out of range ({} references)",
                        n,
                        names.len()
                    );
                    n as u32
                } else {
                    names.iter().position(|n| n == spec).ok_or_else(|| {
                        anyhow::anyhow!("reference '{}' not found (available: {names:?})", spec)
                    })? as u32
                }
            }
            None => 0,
        };
        comp.set_cur_ref_id(ref_id);
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
