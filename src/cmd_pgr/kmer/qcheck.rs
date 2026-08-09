use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use rayon::prelude::*;
use std::io::Write;

/// Build the clap subcommand for qcheck.
pub fn make_subcommand() -> Command {
    Command::new("qcheck")
        .about("Flags error-prone reads from quality-weighted k-mers")
        .after_help(
            r###"
Flags reads that quorum would correct or truncate and keeps the rest
untouched. A quality-weighted k-mer table is built from the input reads
first, then each read is checked for quorum's error signals: no high-quality
anchor, a k-mer with no continuation (truncation), or a base that quorum
would substitute (including the Poisson collision test). No corrected
sequence is produced — the read is kept as-is or discarded.

The output is the kept reads as FASTQ; --discard-file additionally writes
the flagged reads so they can be inspected.

* Supports both plain text and gzipped (.gz) FASTQ files
* Reads from stdin if input file is 'stdin'
* FASTA input is rejected (quality scores are required)

Examples:
1. Filter error-prone reads:
   pgr kmer qcheck reads.fq.gz -k 21 -o kept.fq.gz
2. Also write the discarded reads:
   pgr kmer qcheck reads.fq.gz -k 21 -o kept.fq --discard-file bad.fq
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTQ file to process",
        ))
        .arg(crate::cmd_pgr::args::kmer_arg_with_default("17"))
        .arg(super::qhist::qual_thresh_arg())
        .arg(super::qhist::bits_arg())
        .arg(
            Arg::new("skip")
                .long("skip")
                .num_args(1)
                .default_value("0")
                .value_parser(clap::value_parser!(usize))
                .help("Bases to skip before searching for an anchor"),
        )
        .arg(
            Arg::new("good")
                .long("good")
                .num_args(1)
                .default_value("1")
                .value_parser(clap::value_parser!(usize))
                .help("Consecutive anchor k-mers required"),
        )
        .arg(
            Arg::new("anchor_count")
                .long("anchor-count")
                .num_args(1)
                .default_value("1")
                .value_parser(clap::value_parser!(usize))
                .help("Minimum count for a high-quality anchor k-mer"),
        )
        .arg(
            Arg::new("min_count")
                .long("min-count")
                .num_args(1)
                .default_value("1")
                .value_parser(clap::value_parser!(u64))
                .help("Count above which a base is trusted before the cutoff check"),
        )
        .arg(
            Arg::new("cutoff")
                .long("cutoff")
                .num_args(1)
                .default_value("4")
                .value_parser(clap::value_parser!(u64))
                .help("Trusted count for keeping the current base"),
        )
        .arg(
            Arg::new("apriori_error_rate")
                .long("apriori-error-rate")
                .num_args(1)
                .default_value("0.01")
                .value_parser(clap::value_parser!(f64))
                .help("Prior error rate (collision prob = rate / 3)"),
        )
        .arg(
            Arg::new("poisson_threshold")
                .long("poisson-threshold")
                .num_args(1)
                .default_value("1e-06")
                .value_parser(clap::value_parser!(f64))
                .help("Poisson probability threshold"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("discard_file")
                .long("discard-file")
                .num_args(1)
                .help("Write flagged reads to this FASTQ file"),
        )
}

/// Execute the qcheck command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();
    let k = *args.get_one::<usize>("kmer").unwrap();
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
        "qcheck requires FASTQ reads with quality scores"
    );

    let thresh = match args.get_one::<u8>("qual_thresh").copied() {
        Some(t) => t,
        None => pgr::libs::fq::trim::detect_quality_base(&recs) + 5,
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

    let params = pgr::libs::kmer::qcheck::CheckParams {
        k,
        skip: *args.get_one::<usize>("skip").unwrap(),
        good: *args.get_one::<usize>("good").unwrap(),
        anchor_count: *args.get_one::<usize>("anchor_count").unwrap(),
        min_count: *args.get_one::<u64>("min_count").unwrap(),
        cutoff: *args.get_one::<u64>("cutoff").unwrap(),
        collision_prob: *args.get_one::<f64>("apriori_error_rate").unwrap() / 3.0,
        poisson_threshold: *args.get_one::<f64>("poisson_threshold").unwrap(),
    };

    // Per-read checks are independent and the table is read-only: run them
    // in parallel, then emit both outputs in input order.
    let verdicts: Vec<Result<(), pgr::libs::kmer::qcheck::ReadError>> = recs
        .par_iter()
        .map(|rec| pgr::libs::kmer::qcheck::check_read(&table, rec.sequence(), &params))
        .collect();

    let mut w = pgr::writer(outfile)?;
    let mut discard_w = args
        .get_one::<String>("discard_file")
        .map(|f| pgr::writer(f))
        .transpose()?;
    let mut kept = 0usize;
    let mut flagged = 0usize;
    for (rec, verdict) in recs.iter().zip(&verdicts) {
        let seq = rec.sequence();
        match verdict {
            Ok(()) => {
                write_record(&mut w, rec, seq, rec.quality_scores())?;
                kept += 1;
            }
            Err(e) => {
                if let Some(dw) = &mut discard_w {
                    write_record(dw, rec, seq, rec.quality_scores())?;
                }
                log::debug!("flagged {}: {:?}", rec.name(), e);
                flagged += 1;
            }
        }
    }
    w.flush()?;
    if let Some(dw) = &mut discard_w {
        dw.flush()?;
    }
    log::info!(
        "==> Kept {kept} reads, flagged {flagged} (threshold {thresh}) -> {}",
        outfile
    );
    Ok(())
}

/// Writes a FASTQ record keeping the original name and comment.
fn write_record<W: Write>(
    w: &mut W,
    rec: &pgr::libs::fmt::seq::SeqRecord,
    seq: &[u8],
    qual: &[u8],
) -> std::io::Result<()> {
    let comment = rec.comment();
    if comment.is_empty() {
        writeln!(w, "@{}", rec.name())?;
    } else {
        writeln!(w, "@{} {}", rec.name(), comment)?;
    }
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")?;
    Ok(())
}
