use clap::{Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for hist.
pub fn make_subcommand() -> Command {
    Command::new("hist")
        .about("Builds a k-mer frequency histogram (.hist)")
        .after_help(
            r###"
Builds a k-mer frequency histogram from sequences or an existing table and
writes it in the FastK `.hist` binary layout (readable by Histex, KatGC,
and GenomeScope tooling).

Give either a sequence file (histogram is computed on the fly) or --table
to reuse an existing `.pkt` table. Bins are fixed 1..=32767; counts above
the top bin are folded into it, matching FastK semantics.

* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Histogram from reads:
   pgr kmer hist reads.fq.gz -k 21 -o reads.hist
2. Histogram from an existing table:
   pgr kmer hist -t reads.pkt -o reads.hist
"###,
        )
        .arg(
            Arg::new("infile")
                .num_args(1)
                .index(1)
                .required_unless_present("table")
                .help("Input FASTA/FASTQ file to process (unless --table is given)"),
        )
        .arg(super::profile::table_arg())
        .arg(super::profile::kmer_arg())
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("khist_text")
                .long("khist-text")
                .num_args(1)
                .help("Write the kmercountexact-style text histogram (khist.txt)"),
        )
        .arg(
            Arg::new("peaks")
                .long("peaks")
                .num_args(1)
                .help("Write the kmercountexact-style peaks summary (peaks.txt)"),
        )
}

/// Execute the hist command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = args.get_one::<String>("outfile").unwrap();
    let table_path = args.get_one::<String>("table").map(String::as_str);
    let k = super::resolve_k(args.get_one::<usize>("kmer"), table_path)?;

    let table = if let Some(t) = table_path {
        pgr::libs::kmer::count::load(std::path::Path::new(t), k)?
    } else {
        let infile = args.get_one::<String>("infile").unwrap();
        let seqs = super::read_seqs(infile)?;
        pgr::libs::kmer::count::build_table(&seqs, k)?
    };
    let hist = pgr::libs::kmer::hist::from_table(&table);
    pgr::libs::kmer::hist::write(std::path::Path::new(outfile), &hist)?;
    let khist_text = args.get_one::<String>("khist_text");
    let peaks = args.get_one::<String>("peaks");
    if khist_text.is_some() || peaks.is_some() {
        let hist = pgr::libs::kmer::khist::histogram(&table, pgr::libs::kmer::khist::HIST_MAX);
        if let Some(f) = khist_text {
            let mut w = pgr::writer(f)?;
            pgr::libs::kmer::khist::write_khist_text(
                &mut w,
                &hist,
                pgr::libs::kmer::khist::HIST_MAX,
            )?;
            w.flush()?;
        }
        if let Some(f) = peaks {
            let mut w = pgr::writer(f)?;
            let unique = table.keys.len() as u64;
            let peaks_out = pgr::libs::kmer::khist::call_peaks(&hist);
            pgr::libs::kmer::khist::write_peaks_text(&mut w, &peaks_out, k, unique, &hist)?;
            w.flush()?;
        }
    }
    log::info!(
        "==> Wrote histogram of {} distinct {}-mers to {}",
        table.keys.len(),
        k,
        outfile
    );
    Ok(())
}
