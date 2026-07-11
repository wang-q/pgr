//! Display statistics about a pbit archive.

use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::pbit::decompressor::Decompressor;
use std::io::Write;

/// Build the clap subcommand for stat.
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Displays statistics about a pbit archive")
        .after_help(
            r###"
This command displays information about a pbit archive. By default it
shows an overview; use flags to list samples, reference contigs, or
sample contigs.

Notes:
* pbit files are binary and require random access (seeking)
* Does not support stdin or gzipped inputs

Examples:
1. Show archive overview:
   pgr pbit stat input.pbit

2. List all samples:
   pgr pbit stat input.pbit --samples

3. List reference contigs (with segment counts):
   pgr pbit stat input.pbit --refs

4. List contigs per sample:
   pgr pbit stat input.pbit --contigs

5. List contigs for a specific sample:
   pgr pbit stat input.pbit --contigs -s sample1
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input pbit file to process",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(crate::cmd_pgr::args::pbit_samples_flag_arg())
        .arg(crate::cmd_pgr::args::pbit_refs_flag_arg())
        .arg(crate::cmd_pgr::args::pbit_contigs_flag_arg())
        .arg(crate::cmd_pgr::args::pbit_sample_filter_arg(
            "Restrict --contigs to this sample",
        ))
}

/// Execute the stat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let show_samples = args.get_flag("samples");
    let show_refs = args.get_flag("refs");
    let show_contigs = args.get_flag("contigs");
    let sample_filter = args.get_one::<String>("sample");

    let dec = Decompressor::open(infile)
        .with_context(|| format!("Failed to open pbit file {}", infile))?;
    let mut writer = pgr::libs::io::writer(outfile)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    // Validate --sample early so a non-existent sample is reported clearly
    // rather than silently producing empty output.
    if let Some(s) = sample_filter {
        if !dec.list_samples().contains(&s.as_str()) {
            anyhow::bail!("sample '{}' not found in archive", s);
        }
    }

    // If no flag is set, show the overview.
    let show_overview = !show_samples && !show_refs && !show_contigs;

    if show_overview || show_samples {
        if show_samples {
            for s in dec.list_samples() {
                writeln!(writer, "{}", s)?;
            }
        } else {
            let header = dec.header();
            writeln!(writer, "File: {}", infile)?;
            writeln!(writer, "Version: {}", header.version)?;
            writeln!(writer, "Segment size: {}", header.segment_size)?;
            writeln!(writer, "K-mer length: {}", header.kmer_len)?;
            writeln!(writer, "Min match length: {}", header.min_match_len)?;
            writeln!(writer, "Reference groups: {}", header.ref_group_count)?;
            writeln!(writer, "Samples: {}", header.sample_count)?;
            // Count unique reference contigs.
            let ref_contigs = dec
                .ref_groups()
                .iter()
                .map(|e| e.contig_name.as_str())
                .collect::<std::collections::HashSet<_>>()
                .len();
            writeln!(writer, "Reference contigs: {}", ref_contigs)?;
        }
    }

    if show_refs {
        // Count segments per reference contig.
        let mut ref_counts: indexmap::IndexMap<&str, usize> = indexmap::IndexMap::new();
        for entry in dec.ref_groups() {
            *ref_counts.entry(entry.contig_name.as_str()).or_default() += 1;
        }
        for (name, count) in ref_counts {
            writeln!(writer, "{}\t{}", name, count)?;
        }
    }

    if show_contigs {
        match sample_filter {
            Some(s) => {
                for c in dec.list_contigs(Some(s)) {
                    writeln!(writer, "{}", c)?;
                }
            }
            None => {
                for sample in dec.list_samples() {
                    for c in dec.list_contigs(Some(sample)) {
                        writeln!(writer, "{}\t{}", sample, c)?;
                    }
                }
            }
        }
    }

    writer.flush()?;
    Ok(())
}
