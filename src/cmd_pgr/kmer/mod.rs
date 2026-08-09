//! `pgr kmer` — k-mer table, profile, and histogram generation.

pub mod hist;
pub mod profile;
pub mod table;

use anyhow::Context;
use clap::{ArgMatches, Command};

/// Build the `pgr kmer` subcommand tree.
pub fn make_subcommand() -> Command {
    Command::new("kmer")
        .about("Analyzes k-mer counts, profiles, and frequency histograms")
        .subcommand_required(true)
        .subcommand(table::make_subcommand())
        .subcommand(profile::make_subcommand())
        .subcommand(hist::make_subcommand())
}

/// Dispatch `pgr kmer` subcommands.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("table", sub)) => table::execute(sub),
        Some(("profile", sub)) => profile::execute(sub),
        Some(("hist", sub)) => hist::execute(sub),
        _ => unreachable!("kmer subcommand match"),
    }
}

/// Read all sequences from a FASTA/FASTQ file (plain, gzipped, or stdin).
pub(crate) fn read_seqs(path: &str) -> anyhow::Result<Vec<Vec<u8>>> {
    let mut reader = pgr::libs::fmt::seq::SeqReader::new(path)
        .with_context(|| format!("failed to open {path}"))?;
    let mut rec = pgr::libs::fmt::seq::SeqRecord::new();
    let mut seqs = Vec::new();
    while reader.read_record(&mut rec)? {
        seqs.push(rec.sequence().to_vec());
    }
    Ok(seqs)
}

/// Resolve `--kmer` against an optional `--table`: a table supplies its own
/// k unless the command line pins one that must match.
pub(crate) fn resolve_k(k_arg: Option<&usize>, table_path: Option<&str>) -> anyhow::Result<usize> {
    let table_k = table_path
        .map(|t| pgr::libs::kmer::count::k_of(std::path::Path::new(t)))
        .transpose()?;
    match (k_arg, table_k) {
        (Some(&k), Some(tk)) => {
            anyhow::ensure!(k == tk, "--kmer {k} does not match table k {tk}");
            Ok(k)
        }
        (Some(&k), None) => Ok(k),
        (None, Some(tk)) => Ok(tk),
        (None, None) => anyhow::bail!("--kmer is required when no --table is given"),
    }
}
