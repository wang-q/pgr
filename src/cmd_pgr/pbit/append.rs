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

PAF is mandatory (`--paf`, or the TSV 3rd column): segments covered by PAF
alignments are CIGAR-encoded (pure-match segments become zero-cost Identity
references); uncovered segments fall back to LZ-diff/Raw. An empty PAF file
disables CIGAR encoding (all segments use LZ-diff/Raw).

Notes:
* Sample names are derived from input FASTA basenames (use --name to override)
* If -o is omitted, the input archive is modified in place
* If -o is specified, the input archive is copied to the output path first
* Reference and sample FASTA files may be plain text or gzipped (.gz)
* Contigs in sample FASTA that do not match any reference contig are skipped
* Only ACGTN characters are supported; IUPAC degenerate codes are mapped to N
* `--paf` files are paired with `-i` files by order; `--name` and `--paf`
  are mutually exclusive (use the TSV's 3rd column for PAF, which is required)

Examples:
1. Append a sample in place:
   pgr pbit append archive.pbit -i new_sample.fa -p new_sample.paf

2. Append multiple samples to a new archive:
   pgr pbit append archive.pbit -i s1.fa -p s1.paf -i s2.fa -p s2.paf -o new_archive.pbit

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

    // Guard -o against overwriting a sample input file (stage_work_path copies
    // the archive to -o before samples are read).
    if let Some(out) = outfile_opt {
        let mut inputs: Vec<&str> = Vec::new();
        for (_, path, paf_opt, _) in &samples {
            inputs.push(path.as_str());
            if let Some(paf) = paf_opt {
                inputs.push(paf.as_str());
            }
        }
        if let Some(name_tsv) = args.get_one::<String>("name") {
            inputs.push(name_tsv.as_str());
        }
        crate::cmd_pgr::args::ensure_outfile_distinct(out, inputs)?;
    }

    let in_place = outfile_opt.is_none();
    let (work_path, mut temp_guard) = super::stage_work_path(infile, outfile_opt)?;

    let mut comp = Compressor::open_for_append(&work_path)
        .with_context(|| format!("failed to open pbit archive for append: {}", work_path))?;

    let mut cmd_line = format!("pgr pbit append {}", infile);
    if let Some(out) = outfile_opt {
        cmd_line.push_str(&format!(" -o {}", out));
    }
    // Record the sample inputs for provenance (consistent with `create`).
    for (name, path, paf_opt, ref_spec) in &samples {
        cmd_line.push_str(&format!(" -i {}:{}", name, path));
        if let Some(paf) = paf_opt {
            cmd_line.push_str(&format!(" -p {}", paf));
        }
        if let Some(ref_spec) = ref_spec {
            cmd_line.push_str(&format!(" @ref {}", ref_spec));
        }
    }
    comp.set_cmd_line(&cmd_line);

    // Reject appending a sample whose name already exists in the archive;
    // append_sample would silently merge its segments into the existing sample.
    for (name, _, _, _) in &samples {
        if comp.has_sample(name.as_str()) {
            anyhow::bail!(
                "sample '{}' already exists in the archive; append would corrupt it \
                 (use a distinct --name)",
                name
            );
        }
    }

    let num_refs = comp.ref_names().len();
    for (name, path, paf_opt, ref_spec) in &samples {
        let ref_id = match ref_spec {
            None => {
                if num_refs > 1 {
                    log::warn!(
                        "no reference specified for sample '{}'; defaulting to reference 0 of \
                         {} (route samples with the 4th TSV column or --name)",
                        name,
                        num_refs
                    );
                }
                0
            }
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
