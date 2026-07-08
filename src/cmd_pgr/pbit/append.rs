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

When `--paf` is provided, segments covered by PAF alignments are CIGAR-encoded
(replacing LZ-diff); uncovered segments fall back to LZ-diff.

Notes:
* Sample names are derived from input FASTA basenames (use --name to override)
* If -o is omitted, the input archive is modified in place
* If -o is specified, the input archive is copied to the output path first
* Reference and sample FASTA files may be plain text or gzipped (.gz)
* contigs in sample FASTA that do not match any reference contig are skipped
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
            Arg::new("name").long("name").num_args(1).help(
                "TSV file of `sample_name<TAB>fasta_path[<TAB>paf_path]` lines (overrides -i)",
            ),
        )
        .arg(
            Arg::new("paf")
                .long("paf")
                .short('p')
                .num_args(1)
                .action(ArgAction::Append)
                .help("PAF file(s) for CIGAR-driven encoding (paired with -i by order)"),
        )
}

/// Execute the append command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile_opt = args.get_one::<String>("outfile");

    // Mutex: --name and --paf cannot coexist.
    let has_name = args.get_one::<String>("name").is_some();
    let has_paf = args.get_many::<String>("paf").is_some();
    if has_name && has_paf {
        anyhow::bail!(
            "--name and --paf are mutually exclusive (use --name TSV with 3rd column for PAF)"
        );
    }

    // Collect (sample_name, fasta_path, paf_path_opt) triples.
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
                let name = basename_or_err(path)?;
                let paf = pafs.get(i).cloned();
                pairs.push((name, path.clone(), paf));
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

    Ok(())
}

/// Read a TSV file of `sample_name<TAB>fasta_path[<TAB>paf_path]` lines.
/// Empty lines and lines starting with '#' are skipped. The 3rd column is
/// optional; when present and non-empty, it enables CIGAR-driven encoding.
fn read_name_tsv(path: &str) -> Result<Vec<(String, String, Option<String>)>> {
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
