//! Append a new reference genome to an existing pbit archive.

use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::pbit::compressor::Compressor;
use std::path::PathBuf;

/// Build the clap subcommand for append-ref.
pub fn make_subcommand() -> Command {
    Command::new("append-ref")
        .about("Appends reference genome(s) to an existing pbit archive")
        .after_help(
            r###"
Adds a new reference genome to an existing pbit archive: its 2bit segment
records are appended after the existing reference records, and the Reference
Index + reference table are rewritten. Existing samples and their deltas are
preserved; new samples can route to the new reference (TSV 4th column).

Notes:
* If -o is omitted, the input archive is modified in place (atomic rename).
* If -o is specified, the input archive is copied to the output path first.

Examples:
1. Append a reference in place:
   pgr pbit append-ref cohort.pbit -r newref.fa
2. Append to a new archive:
   pgr pbit append-ref cohort.pbit -r newref.fa -o cohort2.pbit
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Existing pbit archive to append to",
        ))
        .arg(crate::cmd_pgr::args::pbit_ref_arg())
        .arg(crate::cmd_pgr::args::outfile_arg_optional())
}

/// Execute the append-ref command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let ref_fastas: Vec<&str> = args
        .get_many::<String>("ref")
        .context("missing required argument: --ref")?
        .map(|s| s.as_str())
        .collect();
    let outfile_opt = args.get_one::<String>("outfile");

    // Guard -o against overwriting a reference input file (stage_work_path
    // copies the archive to -o before the new reference is read).
    if let Some(out) = outfile_opt {
        crate::cmd_pgr::args::ensure_outfile_distinct(out, ref_fastas.iter().copied())?;
    }

    let in_place = outfile_opt.is_none();
    let (work_path, mut temp_guard) = super::stage_work_path(infile, outfile_opt)?;

    let mut comp = Compressor::open_for_append(&work_path)
        .with_context(|| format!("failed to open pbit archive: {}", work_path))?;
    let mut cmd_line = format!("pgr pbit append-ref {}", infile);
    if let Some(out) = outfile_opt {
        cmd_line.push_str(&format!(" -o {}", out));
    }
    // Record the reference inputs for provenance.
    for ref_fasta in &ref_fastas {
        cmd_line.push_str(&format!(" -r {}", ref_fasta));
    }
    comp.set_cmd_line(&cmd_line);

    for ref_fasta in &ref_fastas {
        comp.append_reference(ref_fasta)?;
    }
    comp.finish().context("failed to finalize pbit archive")?;

    if in_place {
        let rename_from = if let Some(guard) = temp_guard.take() {
            guard
                .disarm()
                .with_context(|| "failed to prepare temp file for in-place rename")?
        } else {
            PathBuf::from(&work_path)
        };
        std::fs::rename(&rename_from, infile).with_context(|| {
            format!(
                "failed to finalize in-place append-ref: rename {} -> {}",
                rename_from.display(),
                infile
            )
        })?;
    }
    Ok(())
}
