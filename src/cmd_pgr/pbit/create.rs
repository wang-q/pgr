use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::pbit::compressor::Compressor;

/// Build the clap subcommand for create.
pub fn make_subcommand() -> Command {
    Command::new("create")
        .about("Creates a pbit archive from a reference FASTA and sample FASTA files")
        .after_help(
            r###"
This command creates a new pbit archive. The reference FASTA is stored as
standard 2bit records; each sample FASTA is LZ-diff encoded against the
matching reference segment, flate2-compressed, and stored as delta entries.

When `--paf` is provided, segments covered by PAF alignments are CIGAR-encoded
(replacing LZ-diff); uncovered segments fall back to LZ-diff.

Notes:
* Sample names are derived from the input FASTA basenames (use `--name` to
  override with a TSV file of `name<TAB>path[<TAB>paf_path]` lines)
* Reference and sample FASTA files may be plain text or gzipped (.gz)
* contigs in sample FASTA that do not match any reference contig are skipped
* Only ACGTN characters are supported; IUPAC degenerate codes (R, Y, S, W,
  K, M, B, D, H, V) are lossily mapped to N
* `--paf` files are paired with `-i` files by order; `--name` and `--paf`
  are mutually exclusive (use the TSV's optional 3rd column for PAF)

Examples:
1. Create a pbit archive with one sample:
   pgr pbit create -r ref.fa -i sample1.fa -o out.pbit

2. Create with multiple samples:
   pgr pbit create -r ref.fa -i s1.fa -i s2.fa -i s3.fa -o out.pbit

3. Custom segment size and k-mer length:
   pgr pbit create -r ref.fa -i sample.fa -o out.pbit -s 8192 -k 15

4. Provide sample names via a TSV file:
   pgr pbit create -r ref.fa --name samples.tsv -o out.pbit

5. CIGAR-driven encoding with PAF:
   pgr pbit create -r ref.fa -i sample.fa -p sample.paf -o out.pbit
"###,
        )
        .arg(
            Arg::new("ref")
                .long("ref")
                .short('r')
                .required(true)
                .num_args(1)
                .help("Reference FASTA file (plain or .gz)"),
        )
        .arg(
            Arg::new("infiles")
                .long("infile")
                .short('i')
                .required(false)
                .num_args(1)
                .action(ArgAction::Append)
                .help("Sample FASTA file(s) (plain or .gz)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("segment_size")
                .long("segment-size")
                .short('s')
                .num_args(1)
                .default_value("4096")
                .value_parser(clap::value_parser!(usize))
                .help("Reference segment size in bp (default: 4096)"),
        )
        .arg(
            Arg::new("kmer_len")
                .long("kmer-len")
                .short('k')
                .num_args(1)
                .default_value("15")
                .value_parser(clap::value_parser!(usize))
                .help("K-mer length for LZ-diff hashing (default: 15)"),
        )
        .arg(
            Arg::new("min_match_len")
                .long("min-match-len")
                .short('l')
                .num_args(1)
                .default_value("18")
                .value_parser(clap::value_parser!(u32))
                .help("Minimum match length for LZ-diff (default: 18)"),
        )
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

/// Execute the create command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let ref_fasta = args
        .get_one::<String>("ref")
        .context("missing required argument: --ref")?;
    let outfile = args
        .get_one::<String>("outfile")
        .context("missing required argument: --outfile")?;
    let segment_size = *args
        .get_one::<usize>("segment_size")
        .context("missing --segment-size")?;
    let kmer_len = *args
        .get_one::<usize>("kmer_len")
        .context("missing --kmer-len")?;
    let min_match_len = *args
        .get_one::<u32>("min_match_len")
        .context("missing --min-match-len")?;

    anyhow::ensure!(segment_size > 0, "segment-size must be positive");
    anyhow::ensure!(kmer_len > 0, "kmer-len must be positive");
    anyhow::ensure!(min_match_len > 0, "min-match-len must be positive");

    let samples = super::collect_samples_from_args(args)?;

    let mut comp = Compressor::create(outfile, ref_fasta, segment_size, kmer_len, min_match_len)
        .with_context(|| format!("failed to create pbit archive: {}", outfile))?;
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
