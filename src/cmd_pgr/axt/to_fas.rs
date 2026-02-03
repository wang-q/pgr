use clap::*;
use intspan::Range;
use pgr::libs::axt::AxtReader;
use pgr::libs::fas::FasEntry;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("to-fas")
        .about("Convert AXT format files to block FA format")
        .after_help(
            r###"
This subcommand converts AXT files into block FA format for further analysis.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- A chromosome sizes file (chr.sizes) for the query genome is required to correctly handle
  coordinates on the negative strand.
- The output file defaults to standard output (stdout). Use the -o option to specify an output file.

Examples:
1. Convert from a file and output to stdout:
   pgr axt to-fas tests/fasr/RM11_1a.chr.sizes tests/fasr/example.axt

2. Read from stdin and output to a file:
   cat tests/fasr/example.axt | pgr axt to-fas tests/fasr/RM11_1a.chr.sizes stdin -o output.fas

3. Specify target and query names:
   pgr axt to-fas tests/fasr/RM11_1a.chr.sizes tests/fasr/example.axt --tname S288c --qname RM11_1a

"###,
        )
        .arg(
            Arg::new("chr.sizes")
                .required(true)
                .index(1)
                .num_args(1)
                .help("Chromosome sizes file for the query genome"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(2)
                .help("Input AXT file(s) to process"),
        )
        .arg(
            Arg::new("tname")
                .long("tname")
                .num_args(1)
                .default_value("target")
                .help("Target name"),
        )
        .arg(
            Arg::new("qname")
                .long("qname")
                .num_args(1)
                .default_value("query")
                .help("Query name"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let sizes = intspan::read_sizes(args.get_one::<String>("chr.sizes").unwrap());

    let tname = args.get_one::<String>("tname").unwrap();
    let qname = args.get_one::<String>("qname").unwrap();

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        let axt_iter = AxtReader::new(reader);

        for axt_res in axt_iter {
            let axt = axt_res?;

            //----------------------------
            // Output
            //----------------------------
            // Target Entry
            let t_start = (axt.t_start + 1) as i32;
            let t_end = axt.t_end as i32;
            let mut t_range = Range::from(&axt.t_name, t_start, t_end);
            *t_range.name_mut() = tname.to_string();
            *t_range.strand_mut() = "+".to_string();

            let t_entry = FasEntry::from(&t_range, axt.t_sym.as_bytes());
            writer.write_all(t_entry.to_string().as_bytes())?;

            // Query Entry
            let q_len = *sizes.get(&axt.q_name).ok_or_else(|| {
                anyhow::anyhow!(".sizes file doesn't contain the needed chr: {}", axt.q_name)
            })?;

            let (q_start, q_end) = if axt.q_strand == '-' {
                let q_s_1 = (axt.q_start + 1) as i32;
                let q_e_1 = axt.q_end as i32;

                let fwd_start = q_len - q_e_1 + 1;
                let fwd_end = q_len - q_s_1 + 1;

                (fwd_start, fwd_end)
            } else {
                ((axt.q_start + 1) as i32, axt.q_end as i32)
            };

            let mut q_range = Range::from(&axt.q_name, q_start, q_end);
            *q_range.name_mut() = qname.to_string();
            *q_range.strand_mut() = axt.q_strand.to_string();

            let q_entry = FasEntry::from(&q_range, axt.q_sym.as_bytes());
            writer.write_all(q_entry.to_string().as_bytes())?;

            // Add a newline to separate blocks
            writer.write_all(b"\n")?;
        }
    }

    Ok(())
}
