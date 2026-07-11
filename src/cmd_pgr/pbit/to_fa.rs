//! Convert a pbit archive to per-sample FASTA files.

use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::pbit::decompressor::Decompressor;
use std::io::Write;
use std::path::Path;

/// Build the clap subcommand for to-fa.
pub fn make_subcommand() -> Command {
    Command::new("to-fa")
        .about("Converts a pbit archive to per-sample FASTA files")
        .after_help(
            r###"
This command extracts all sample sequences from a pbit archive and writes
them as FASTA files (one per sample) into the specified output directory.

Notes:
* pbit files are binary and require random access (seeking)
* Does not support stdin or gzipped inputs
* One output file per sample: `{outdir}/{sample_name}.fa`
* Sequence lines are wrapped at 60 bases

Examples:
1. Convert pbit to FASTA (one file per sample):
   pgr pbit to-fa input.pbit -o outdir

2. Extract a single sample:
   pgr pbit to-fa input.pbit -o outdir -s sample1
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input pbit file to process",
        ))
        .arg(crate::cmd_pgr::args::outdir_arg_required())
        .arg(crate::cmd_pgr::args::pbit_sample_filter_arg(
            "Extract only this sample (default: all samples)",
        ))
}

/// Execute the to-fa command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let outdir = args
        .get_one::<String>("outdir")
        .context("missing required argument: outdir")?;
    let sample_filter = args.get_one::<String>("sample");

    // Create output directory if it does not exist.
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("failed to create output directory: {}", outdir))?;

    let mut dec = Decompressor::open(input_path)
        .with_context(|| format!("Failed to open pbit file {}", input_path))?;

    let samples: Vec<String> = match sample_filter {
        Some(s) => {
            if !dec.list_samples().contains(&s.as_str()) {
                anyhow::bail!("sample '{}' not found in archive", s);
            }
            vec![s.clone()]
        }
        None => dec.list_samples().into_iter().map(String::from).collect(),
    };

    for sample in &samples {
        // Guard against path traversal: sample names come from the archive
        // (untrusted .pbit input) and must not escape the output directory.
        if sample.contains('/') || sample.contains('\\') || sample == "." || sample == ".." {
            anyhow::bail!(
                "invalid sample name '{}': must not contain path separators or be '.'/'..'",
                sample
            );
        }
        let out_path = Path::new(outdir).join(format!("{}.fa", sample));
        let out_str = out_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("invalid output path: {}", out_path.display()))?;
        let mut writer = pgr::libs::io::writer(out_str)
            .with_context(|| format!("failed to open output file: {}", out_str))?;
        dec.get_sample(sample, &mut writer)
            .with_context(|| format!("failed to extract sample '{}'", sample))?;
        writer.flush()?;
    }

    Ok(())
}
