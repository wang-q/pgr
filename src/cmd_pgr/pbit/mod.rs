pub mod append;
pub mod create;
pub mod range;
pub mod some;
pub mod stat;
pub mod to_fa;

use anyhow::{Context, Result};
use clap::{ArgMatches, Command};

/// Read a TSV file of `sample_name<TAB>fasta_path[<TAB>paf_path]` lines.
/// Empty lines and lines starting with '#' are skipped. The 3rd column is
/// optional; when present and non-empty, it enables CIGAR-driven encoding.
pub(crate) fn read_name_tsv(path: &str) -> Result<Vec<(String, String, Option<String>)>> {
    let lines = pgr::libs::io::read_lines(path)
        .with_context(|| format!("failed to read name TSV: {}", path))?;
    let mut out = Vec::new();
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        let name = parts
            .first()
            .ok_or_else(|| anyhow::anyhow!("missing sample name in line: {}", line))?
            .trim()
            .to_string();
        let fasta_path = parts
            .get(1)
            .ok_or_else(|| anyhow::anyhow!("missing FASTA path in line: {}", line))?
            .trim()
            .to_string();
        let paf_path = parts
            .get(2)
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty());
        if name.is_empty() || fasta_path.is_empty() {
            anyhow::bail!("empty name or path in line: {}", line);
        }
        out.push((name, fasta_path, paf_path));
    }
    Ok(out)
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
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
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
