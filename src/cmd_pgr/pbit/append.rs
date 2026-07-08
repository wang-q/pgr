use anyhow::{Context, Result};
use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::io::basename_or_err;
use pgr::libs::pbit::compressor::Compressor;

/// Build the clap subcommand for append.
pub fn make_subcommand() -> Command {
    Command::new("append")
        .about("Appends samples to an existing pbit archive")
        .after_help(
            r###"
This command appends new sample FASTA files to an existing pbit archive.
The reference is already embedded in the archive, so no -r is needed.

Notes:
* Sample names are derived from input FASTA basenames (use --name to override)
* If -o is omitted, the input archive is modified in place
* If -o is specified, the input archive is copied to the output path first
* Reference and sample FASTA files may be plain text or gzipped (.gz)
* contigs in sample FASTA that do not match any reference contig are skipped
* Only ACGTN characters are supported; IUPAC degenerate codes are mapped to N

Examples:
1. Append a sample in place:
   pgr pbit append archive.pbit -i new_sample.fa

2. Append multiple samples to a new archive:
   pgr pbit append archive.pbit -i s1.fa -i s2.fa -o new_archive.pbit

3. Provide sample names via TSV:
   pgr pbit append archive.pbit --name samples.tsv -o new_archive.pbit
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Existing pbit archive to append to",
        ))
        .arg(
            Arg::new("infiles")
                .long("infile")
                .short('i')
                .required(false)
                .num_args(1)
                .action(ArgAction::Append)
                .help("Sample FASTA file(s) to append (plain or .gz)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_optional())
        .arg(
            Arg::new("name")
                .long("name")
                .num_args(1)
                .help("TSV file of `sample_name<TAB>fasta_path` lines (overrides -i)"),
        )
}

/// Execute the append command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile_opt = args.get_one::<String>("outfile");

    // Collect (sample_name, fasta_path) pairs.
    let samples: Vec<(String, String)> = if let Some(name_tsv) = args.get_one::<String>("name") {
        read_name_tsv(name_tsv)?
    } else {
        let infiles = args
            .get_many::<String>("infiles")
            .ok_or_else(|| anyhow::anyhow!("no input files: provide -i or --name"))?;
        let mut pairs = Vec::new();
        for path in infiles {
            let name = basename_or_err(path)?;
            pairs.push((name, path.clone()));
        }
        pairs
    };

    if samples.is_empty() {
        anyhow::bail!("no sample FASTA files provided");
    }

    // Determine the working path: copy if -o specified, else in-place.
    let work_path = match outfile_opt {
        Some(out) => {
            std::fs::copy(infile, out)
                .with_context(|| format!("failed to copy {} to {}", infile, out))?;
            out.clone()
        }
        None => infile.clone(),
    };

    let mut comp = Compressor::open_for_append(&work_path)
        .with_context(|| format!("failed to open pbit archive for append: {}", work_path))?;
    for (name, path) in &samples {
        comp.append_sample(name, path)
            .with_context(|| format!("failed to append sample '{}'", name))?;
    }
    comp.finish().context("failed to finalize pbit archive")?;

    Ok(())
}

/// Read a TSV file of `sample_name<TAB>fasta_path` lines.
/// Empty lines and lines starting with '#' are skipped.
fn read_name_tsv(path: &str) -> Result<Vec<(String, String)>> {
    let lines = pgr::libs::io::read_lines(path)
        .with_context(|| format!("failed to read name TSV: {}", path))?;
    let mut out = Vec::new();
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut parts = trimmed.splitn(2, '\t');
        let name = parts
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing sample name in line: {}", line))?
            .trim()
            .to_string();
        let fasta_path = parts
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing FASTA path in line: {}", line))?
            .trim()
            .to_string();
        if name.is_empty() || fasta_path.is_empty() {
            anyhow::bail!("empty name or path in line: {}", line);
        }
        out.push((name, fasta_path));
    }
    Ok(out)
}
