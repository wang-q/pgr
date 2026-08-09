use clap::{Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for gsize.
pub fn make_subcommand() -> Command {
    Command::new("gsize")
        .about("Estimates coverage and genome size from k-mer frequencies")
        .after_help(
            r###"
Estimates the k-mer coverage peak and genome size from a count table.
`peak_coverage` is the frequency carried by the most distinct k-mers (the
main mode); `genome_size` is total k-mer instances / peak coverage — the
cheap haploid estimate that precedes GenomeScope-style model fitting.

Give either a sequence file (table built on the fly) or --table to reuse an
existing `.pkt` table.

* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Estimate from reads:
   pgr kmer gsize reads.fq.gz -k 21
2. Estimate from an existing table:
   pgr kmer gsize -t reads.pkt -o stats.tsv
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
        .arg(
            Arg::new("model")
                .long("model")
                .action(clap::ArgAction::SetTrue)
                .help("Fit the GenomeScope model (kmercov/het/genome size)"),
        )
        .arg(
            Arg::new("ploidy")
                .long("ploidy")
                .short('p')
                .num_args(1)
                .default_value("1")
                .value_parser(clap::value_parser!(usize))
                .help("Ploidy for the model (1 or 2; default 1)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the gsize command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let table_path = args.get_one::<String>("table").map(String::as_str);
    let k = super::resolve_k(args.get_one::<usize>("kmer"), table_path)?;

    let table = if let Some(t) = table_path {
        pgr::libs::kmer::count::load(std::path::Path::new(t), k)?
    } else {
        let infile = args.get_one::<String>("infile").unwrap();
        let seqs = super::read_seqs(infile)?;
        pgr::libs::kmer::count::build_table(&seqs, k)?
    };

    if args.get_flag("model") {
        let p = *args.get_one::<usize>("ploidy").unwrap();
        anyhow::ensure!((1..=2).contains(&p), "ploidy must be 1 or 2, got {p}");
        let hist = pgr::libs::kmer::hist::from_table(&table);
        let model = pgr::libs::kmer::genomescope::fit(&hist.hist, k, p);
        // `-o` is the output directory holding summary.txt and model.txt in
        // the GenomeScope formats consumed by anchr's 2_fastk; the summary
        // goes to stdout.
        let outdir = if outfile == "stdout" { "." } else { outfile };
        pgr::libs::kmer::genomescope::write_outputs(std::path::Path::new(outdir), &model)?;
        let mut w = std::io::stdout();
        writeln!(w, "k\t{k}")?;
        writeln!(w, "kmercov\t{:.1}", model.kmercov)?;
        writeln!(w, "bias\t{:.3}", model.bias)?;
        writeln!(w, "d\t{:.4}", model.d)?;
        writeln!(w, "genome_size\t{:.1}", model.length)?;
        writeln!(w, "het_fraction\t{:.4}", model.het)?;
        writeln!(w, "converged\t{}", model.converged)?;
        log::info!(
            "==> GenomeScope fit: kmercov {:.1}, het {:.3}, genome size {:.0} bp -> {}",
            model.kmercov,
            model.het,
            model.length,
            outdir
        );
    } else {
        let mut w = pgr::writer(outfile)?;
        writeln!(w, "k\t{k}")?;
        let est = pgr::libs::kmer::hist::estimate(&table);
        writeln!(w, "peak_coverage\t{}", est.peak_cov)?;
        writeln!(w, "total_distinct\t{}", est.total_distinct)?;
        writeln!(w, "total_kmers\t{}", est.total_kmers)?;
        writeln!(w, "genome_size\t{:.1}", est.genome_size)?;
        log::info!(
            "==> peak coverage {}, estimated genome size {:.0} bp",
            est.peak_cov,
            est.genome_size
        );
    }
    Ok(())
}
