pub mod append;
pub mod create;
pub mod range;
pub mod some;
pub mod stat;
pub mod to_fa;

use anyhow::{Context, Result};
use clap::{ArgMatches, Command};

/// Read a TSV of `sample_name<TAB>fasta_path[<TAB>paf_path]` lines (3rd column optional).
pub(crate) fn read_name_tsv(path: &str) -> Result<Vec<(String, String, Option<String>)>> {
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
        let name = parts
            .first()
            .ok_or_else(|| anyhow::anyhow!("line {}: missing sample name: {}", line_no, trimmed))?
            .trim()
            .to_string();
        let fasta_path = parts
            .get(1)
            .ok_or_else(|| anyhow::anyhow!("line {}: missing FASTA path: {}", line_no, trimmed))?
            .trim()
            .to_string();
        let paf_path = parts
            .get(2)
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty());
        if name.is_empty() || fasta_path.is_empty() {
            anyhow::bail!("line {}: empty name or path: {}", line_no, trimmed);
        }
        out.push((name, fasta_path, paf_path));
    }
    Ok(out)
}

/// Collect `(sample_name, fasta_path, paf_path_opt)` triples from `--name` TSV or `-i`/`--paf`.
pub(crate) fn collect_samples_from_args(
    args: &ArgMatches,
) -> Result<Vec<(String, String, Option<String>)>> {
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

    let samples: Vec<(String, String, Option<String>)> =
        if let Some(name_tsv) = args.get_one::<String>("name") {
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
                pairs.push((name, path.clone(), paf));
            }
            pairs
        };

    if samples.is_empty() {
        anyhow::bail!("no sample FASTA files provided");
    }

    Ok(samples)
}

/// Build the clap subcommand for pbit.
pub fn make_subcommand() -> Command {
    Command::new("pbit")
        .about("Manages pbit (population 2bit + delta) files")
        .after_help(
            r###"Subcommand groups:

* build:     create / append
* info:      stat
* subset:    range / some
* transform: to-fa

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(create::make_subcommand())
        .subcommand(append::make_subcommand())
        .subcommand(range::make_subcommand())
        .subcommand(some::make_subcommand())
        .subcommand(stat::make_subcommand())
        .subcommand(to_fa::make_subcommand())
}

/// Execute the pbit command.
pub fn execute(args: &ArgMatches) -> Result<()> {
    match args.subcommand() {
        Some(("create", sub_matches)) => create::execute(sub_matches),
        Some(("append", sub_matches)) => append::execute(sub_matches),
        Some(("range", sub_matches)) => range::execute(sub_matches),
        Some(("some", sub_matches)) => some::execute(sub_matches),
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        Some(("to-fa", sub_matches)) => to_fa::execute(sub_matches),
        _ => Ok(()),
    }
}
