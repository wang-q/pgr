use anyhow::Context;
use clap::{ArgMatches, Command};
use intspan::Range;
use pgr::libs::fmt::axt::AxtReader;
use pgr::libs::fmt::fas::FasEntry;
use std::io::Write;

/// Build the clap subcommand for to-fas.
pub fn make_subcommand() -> Command {
    Command::new("to-fas")
        .about("Converts AXT format files to block FA format")
        .after_help(
            r###"
This subcommand converts AXT files into block FA format for further analysis.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* A chromosome sizes file (chr.sizes) for the query genome is required to correctly handle coordinates on the negative strand
* The output file defaults to standard output (stdout). Use the -o option to specify an output file

Examples:
1. Convert from a file and output to stdout:
   pgr axt to-fas tests/fasr/RM11_1a.chr.sizes tests/fasr/example.axt

2. Read from stdin and output to a file:
   cat tests/fasr/example.axt | pgr axt to-fas tests/fasr/RM11_1a.chr.sizes stdin -o output.fas

3. Specify target and query names:
   pgr axt to-fas tests/fasr/RM11_1a.chr.sizes tests/fasr/example.axt --t-name S288c --q-name RM11_1a

"###,
        )
        .arg(crate::cmd_pgr::args::chain_q_sizes_arg().index(1))
        .arg(crate::cmd_pgr::args::infiles_arg_at("AXT", 2))
        .arg(crate::cmd_pgr::args::t_name_arg(Some("target")))
        .arg(crate::cmd_pgr::args::q_name_arg(Some("query")))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the to-fas command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let q_sizes_path = args.get_one::<String>("q_sizes").unwrap();
    let sizes = pgr::read_sizes::<i32>(q_sizes_path)
        .with_context(|| format!("Failed to read sizes from {}", q_sizes_path))?;

    let tname = args.get_one::<String>("t_name").unwrap();
    let qname = args.get_one::<String>("q_name").unwrap();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
        let axt_iter = AxtReader::new(reader);

        for axt_res in axt_iter {
            let axt = axt_res?;

            // Target Entry
            let t_start = axt
                .t_start
                .checked_add(1)
                .and_then(|v| i32::try_from(v).ok())
                .ok_or_else(|| anyhow::anyhow!("t_start {} exceeds i32 range", axt.t_start))?;
            let t_end = i32::try_from(axt.t_end)
                .map_err(|_| anyhow::anyhow!("t_end {} exceeds i32 range", axt.t_end))?;
            let mut t_range = Range::from(&axt.t_name, t_start, t_end);
            *t_range.name_mut() = tname.to_string();
            *t_range.strand_mut() = "+".to_string();

            let t_entry = FasEntry::from(&t_range, axt.t_sym.as_bytes());
            writer.write_all(t_entry.to_string().as_bytes())?;

            // Query Entry
            let q_len = *sizes.get(&axt.q_name).ok_or_else(|| {
                anyhow::anyhow!(".sizes file doesn't contain the needed chr: {}", axt.q_name)
            })?;

            let (q_start, q_end) = pgr::libs::fmt::axt::axt_query_to_forward_coords(
                axt.q_start,
                axt.q_end,
                axt.q_strand,
                q_len,
            )?;

            let mut q_range = Range::from(&axt.q_name, q_start, q_end);
            *q_range.name_mut() = qname.to_string();
            *q_range.strand_mut() = axt.q_strand.to_string();

            let q_entry = FasEntry::from(&q_range, axt.q_sym.as_bytes());
            writer.write_all(q_entry.to_string().as_bytes())?;

            // Add a newline to separate blocks
            writer.write_all(b"\n")?;
        }
    }

    writer.flush()?;
    Ok(())
}
