use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for qhist.
pub fn make_subcommand() -> Command {
    Command::new("qhist")
        .about("Builds a quality-weighted k-mer histogram (quorum-style)")
        .after_help(
            r###"
Builds the quality-weighted k-mer frequency histogram from FASTQ reads and
writes it in quorum's `histo_mer_database` format: `count n_lowq n_highq`
per non-empty count bin.

A k-mer counts as high quality iff all k bases of its window score at least
the quality threshold; per k-mer the final count is the number of
high-quality occurrences when any exist, otherwise the number of
low-quality occurrences (quorum `hash_with_quality` semantics: low-quality
evidence never raises a high-quality count). Counts are capped at 1000.
Counts are additionally capped by --bits (default 7, quorum's
create_database default: max count 127).

The threshold defaults to the detected Phred offset (+33/+64) plus 5
(quorum's default min-quality offset) and can be pinned with --qual-thresh.

* Supports both plain text and gzipped (.gz) FASTQ files
* Reads from stdin if input file is 'stdin'
* FASTA input is rejected (quality scores are required)

Examples:
1. Auto-detected threshold:
   pgr kmer qhist reads.fq.gz -k 21 -o reads.qhist
2. Pinned threshold (ASCII value; Phred+33 Q10 = 43):
   pgr kmer qhist reads.fq.gz -k 21 -q 43 -o reads.qhist
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTQ file to process",
        ))
        .arg(crate::cmd_pgr::args::kmer_arg_with_default("17"))
        .arg(qual_thresh_arg())
        .arg(bits_arg())
        .arg(crate::cmd_pgr::args::outfile_arg_required())
}

/// Optional `-q/--qual-thresh` argument for quality-weighted commands.
pub fn qual_thresh_arg() -> Arg {
    Arg::new("qual_thresh")
        .long("qual-thresh")
        .short('q')
        .num_args(1)
        .value_parser(clap::value_parser!(u8))
        .help("Quality ASCII threshold (default: detected Phred offset + 5)")
}

/// Optional `-b/--bits` argument for quality-weighted commands.
pub fn bits_arg() -> Arg {
    Arg::new("bits")
        .long("bits")
        .short('b')
        .num_args(1)
        .default_value("7")
        .value_parser(clap::value_parser!(u8))
        .help("Count bits (quorum create_database -b; max count = 2^bits - 1)")
}

/// Execute the qhist command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();
    let k = *args.get_one::<usize>("kmer").unwrap();
    anyhow::ensure!(
        k > 0 && k <= pgr::libs::kmer::key::Kmer::MAX_K,
        "k must be in 1..={}, got {k}",
        pgr::libs::kmer::key::Kmer::MAX_K
    );
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, [infile.as_str()])?;

    let mut reader = pgr::libs::fmt::seq::SeqReader::new(infile)
        .with_context(|| format!("failed to open {infile}"))?;
    let mut rec = pgr::libs::fmt::seq::SeqRecord::new();
    let mut recs = Vec::new();
    while reader.read_record(&mut rec)? {
        recs.push(rec.clone());
    }
    anyhow::ensure!(!recs.is_empty(), "no reads in {infile}");
    anyhow::ensure!(
        recs.iter().all(|r| r.is_fastq()),
        "qhist requires FASTQ reads with quality scores"
    );

    let thresh = match args.get_one::<u8>("qual_thresh").copied() {
        Some(t) => t,
        None => pgr::libs::fq::qual::detect_quality_base(&recs) + 5,
    };
    let bits = *args.get_one::<u8>("bits").unwrap();
    anyhow::ensure!(
        (1..=63).contains(&bits),
        "bits must be in 1..=63, got {bits}"
    );
    let count_cap = (1u64 << bits) - 1;
    let seqs: Vec<Vec<u8>> = recs.iter().map(|r| r.sequence().to_vec()).collect();
    let quals: Vec<Vec<u8>> = recs.iter().map(|r| r.quality_scores().to_vec()).collect();

    let table = pgr::libs::kmer::quality::build_table(&seqs, &quals, k, thresh, count_cap);
    let hist = pgr::libs::kmer::quality::histogram(&table);
    let mut w = pgr::writer(outfile)?;
    pgr::libs::kmer::quality::write_hist(&mut w, &hist)?;
    w.flush()?;
    log::info!(
        "==> Wrote quality-weighted k-mer histogram (threshold {thresh}) to {}",
        outfile
    );
    Ok(())
}
