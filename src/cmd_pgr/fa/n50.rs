use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use pgr::libs::fasta::stat::{calc_n50_stats, transpose};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("n50")
        .about("Calculates N50 and other statistics")
        .after_help(
            r###"
This command calculates various assembly statistics from FASTA files.

Statistics:
* N50/N90: Length where contigs of this length or longer include 50%/90% of the total
* S: Sum of all sequence lengths
* A: Average sequence length
* E: E-size, the expected contig length at which a random base occurs
* C: Count of sequences

Notes:
* N50 is calculated by default, use `-N 0` to skip
* Multiple N-statistics: `-N 50 -N 90`
* Use --genome-size to calculate statistics based on estimated genome size
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Basic N50 calculation:
   pgr fa n50 input.fa

2. Calculate N50, N90 and other statistics:
   pgr fa n50 input.fa -N 50 -N 90 -S -A -E -C

3. Calculate based on genome size:
   pgr fa n50 input.fa -g 3000000

4. Transpose output for better readability:
   pgr fa n50 input.fa -N 50 -N 90 -S --transpose

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("FASTA"))
        .arg(
            Arg::new("no_header")
                .long("no-header")
                .short('H')
                .action(ArgAction::SetTrue)
                .help("Do not display headers"),
        )
        .arg(
            Arg::new("nx")
                .long("nx")
                .short('N')
                .num_args(1)
                .default_value("50")
                .action(ArgAction::Append)
                .value_parser(value_parser!(usize))
                .help("Compute Nx statistic"),
        )
        .arg(
            Arg::new("sum")
                .long("sum")
                .short('S')
                .action(ArgAction::SetTrue)
                .help("Compute the sum of the sizes of all records"),
        )
        .arg(
            Arg::new("average")
                .long("average")
                .short('A')
                .action(ArgAction::SetTrue)
                .help("Compute the average length of all records"),
        )
        .arg(
            Arg::new("esize")
                .long("esize")
                .short('E')
                .action(ArgAction::SetTrue)
                .help("Compute the E-size (from GAGE)"),
        )
        .arg(crate::cmd_pgr::args::count_arg("Count records"))
        .arg(
            Arg::new("genome_size")
                .long("genome-size")
                .short('g')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Size of the genome, not the total size of the files"),
        )
        .arg(
            Arg::new("transpose")
                .long("transpose")
                .action(ArgAction::SetTrue)
                .help("Transpose the outputs"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let is_noheader = args.get_flag("no_header");
    let is_sum = args.get_flag("sum");
    let is_average = args.get_flag("average");
    let is_esize = args.get_flag("esize");
    let is_count = args.get_flag("count");
    let is_transpose = args.get_flag("transpose");

    let opt_nx: Vec<_> = args.get_many::<usize>("nx").unwrap().copied().collect();
    let opt_genome = args.get_one::<usize>("genome_size").copied();
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    //----------------------------
    // Process
    //----------------------------
    let mut lens = vec![];

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut fa_in = pgr::libs::fmt::fa::reader(infile)?;

        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;

            let len = record.sequence().len();

            lens.push(len);
        }
    }

    let stats = calc_n50_stats(lens, &opt_nx, opt_genome);

    //----------------------------
    // Output
    //----------------------------
    let mut outputs = vec![];

    // set N == 0 to skip this
    if !(opt_nx.len() == 1 && opt_nx[0] == 0) {
        for (i, nx) in opt_nx.iter().enumerate() {
            let mut row = vec![];
            if !is_noheader {
                row.push(format!("N{}", nx));
            }
            row.push(format!("{}", stats.nx_sizes[i]));
            outputs.push(row);
        }
    }

    if is_sum {
        let mut row = vec![];
        if !is_noheader {
            row.push("S".to_string());
        }
        row.push(format!("{}", stats.total_size));
        outputs.push(row);
    }

    if is_average {
        let mut row = vec![];
        if !is_noheader {
            row.push("A".to_string());
        }
        row.push(format!(
            "{:.2}",
            stats.total_size as f64 / stats.record_cnt as f64
        ));
        outputs.push(row);
    }

    if is_esize {
        let mut row = vec![];
        if !is_noheader {
            row.push("E".to_string());
        }
        row.push(format!("{:.2}", stats.e_size));
        outputs.push(row);
    }

    if is_count {
        let mut row = vec![];
        if !is_noheader {
            row.push("C".to_string());
        }
        row.push(format!("{}", stats.record_cnt));
        outputs.push(row);
    }

    if is_transpose {
        outputs = transpose(outputs);
    }

    for row in outputs {
        writer.write_fmt(format_args!("{}\n", row.join("\t")))?;
    }

    Ok(())
}
