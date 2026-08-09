use clap::{Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for gc.
pub fn make_subcommand() -> Command {
    Command::new("gc")
        .about("Builds a GC-content vs k-mer coverage matrix (.kgc)")
        .after_help(
            r###"
Builds the two-dimensional GC-content × k-mer coverage matrix from a k-mer
table and writes it in the KatGC `.kgc` format (GCP/KF/Count rows, 2x2
neighbor average, values clamped to the peak). Rows are GC counts, columns
are count bins, so the matrix shows how k-mer coverage varies with GC
content (a typical quality diagnostic for read sets).

Give either a sequence file (table built on the fly) or --table to reuse an
existing `.pkt` table. Counts above the count cap are folded into the top
bin; the output x-range is the peak coverage times --xrel (default 2.1),
unless --xmax pins it absolutely (which also sets the count cap, as in
KatGC).

* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Matrix from reads:
   pgr kmer gc reads.fq.gz -k 21 -o reads.kgc
2. Matrix from an existing table:
   pgr kmer gc -t reads.pkt -x 1.9 -o reads.kgc
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
            Arg::new("xmax")
                .long("xmax")
                .short('X')
                .num_args(1)
                .value_parser(clap::value_parser!(usize))
                .help("Absolute x max (also caps the count axis; default: auto)"),
        )
        .arg(
            Arg::new("xrel")
                .long("xrel")
                .short('x')
                .num_args(1)
                .default_value("2.1")
                .value_parser(clap::value_parser!(f64))
                .help("Max x as a multiple of the peak coverage (default: 2.1)"),
        )
        .arg(
            Arg::new("tex")
                .long("tex")
                .action(clap::ArgAction::SetTrue)
                .help("Render a LaTeX heatmap instead of the .kgc matrix"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
}

/// Execute the gc command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = args.get_one::<String>("outfile").unwrap();
    let table_path = args.get_one::<String>("table").map(String::as_str);
    let k = super::resolve_k(args.get_one::<usize>("kmer"), table_path)?;
    let xmax_arg = args.get_one::<usize>("xmax").copied();
    let xrel = *args.get_one::<f64>("xrel").unwrap();

    let table = if let Some(t) = table_path {
        pgr::libs::kmer::count::load(std::path::Path::new(t), k)?
    } else {
        let infile = args.get_one::<String>("infile").unwrap();
        let seqs = super::read_seqs(infile)?;
        pgr::libs::kmer::count::build_table(&seqs, k)?
    };

    let hmax = xmax_arg.unwrap_or(1000);
    let plot = pgr::libs::kmer::gc::gc_matrix(&table, hmax);
    let peak = pgr::libs::kmer::gc::find_peak(&plot, hmax)?;
    let xmax = xmax_arg.unwrap_or_else(|| pgr::libs::kmer::gc::x_limit(peak, xrel, hmax));

    let mut w = pgr::writer(outfile)?;
    if args.get_flag("tex") {
        let hm = pgr::libs::kmer::gc::heatmap(&plot, xmax, peak.zmax);
        let mut context = tera::Context::new();
        context.insert("table", &hm.table);
        context.insert("xlabel", "k-mer coverage");
        context.insert("ylabel", "GC content");
        context.insert("width", &hm.width);
        context.insert("height", &hm.height);
        context.insert("xticks", &hm.xticks);
        context.insert("xtick_labels", &hm.xtick_labels);
        context.insert("ygroups", &hm.ygroups);
        context.insert("yticks", &hm.yticks);
        context.insert("label_len", &hm.label_len);
        pgr::libs::plot::histogram::render_hh_tex(&context, &mut w)?;
    } else {
        pgr::libs::kmer::gc::write_kgc(&mut w, &plot, xmax, peak.zmax)?;
    }
    w.flush()?;
    log::info!(
        "==> Wrote GC x coverage matrix ({} GC rows x {} count bins, peak {}) to {}",
        k,
        xmax,
        peak.xmax,
        outfile
    );
    Ok(())
}
