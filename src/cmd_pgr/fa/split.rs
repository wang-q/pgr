use anyhow::Context;
use clap::{builder::PossibleValue, value_parser, Arg, ArgAction, ArgMatches, Command};
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};

/// Build the clap subcommand for split.
pub fn make_subcommand() -> Command {
    Command::new("split")
        .about("Splits FASTA file(s) into several files")
        .after_help(
            r###"
Split FASTA files into multiple smaller files based on different modes:

1. name: Create separate files for each sequence
   * Uses sequence names as filenames (sanitized)
   * Special characters /\\(): are replaced with _

2. about: Split by approximate size
   * -c SIZE: Split into files of about SIZE bytes each
   * -e: Ensure even number of sequences per file
   * --max-part NUM: Maximum number of output files (default: 999)

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Output files are named as xxx.fa
* For 'name' mode, filenames are sanitized
* For 'about' mode, files are zero-padded numbered

Examples:
1. Split by sequence names:
   pgr fa split name input.fa -o output_dir

2. Split into ~1MB files:
   pgr fa split about input.fa -c 1000000 -o output_dir

3. Split with even sequences:
   pgr fa split about input.fa -c 1000000 -e -o output_dir


"###,
        )
        .arg(
            Arg::new("split_mode")
                .required(true)
                .index(1)
                .action(ArgAction::Set)
                .value_parser([PossibleValue::new("name"), PossibleValue::new("about")])
                .help("Split mode: 'name' or 'about'"),
        )
        .arg(crate::cmd_pgr::args::infiles_arg_at("FASTA", 2))
        .arg(crate::cmd_pgr::args::chunk_size_arg(
            None,
            "Approximate size in bytes (for 'about' mode)",
        ))
        .arg(
            Arg::new("even")
                .long("even")
                .short('e')
                .action(ArgAction::SetTrue)
                .help("Record number in one file should be EVEN"),
        )
        .arg(
            Arg::new("max_part")
                .long("max-part")
                .num_args(1)
                .default_value("999")
                .value_parser(value_parser!(usize))
                .help("Maximum number of output files"),
        )
        .arg(crate::cmd_pgr::args::outdir_arg())
}

/// Execute the split command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mode = args.get_one::<String>("split_mode").unwrap();

    let outdir = args.get_one::<String>("outdir").unwrap();
    if outdir != "stdout" {
        std::fs::create_dir_all(outdir)?;
    }

    let mut fh_of: BTreeMap<String, BufWriter<std::fs::File>> = BTreeMap::new();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();

    // Operating
    if mode == "name" {
        for infile in args.get_many::<String>("infiles").unwrap() {
            let mut fa_in = pgr::libs::fmt::fa::reader(infile)
                .with_context(|| format!("Failed to open reader for {}", infile))?;

            for result in fa_in.records() {
                // obtain record or fail with error
                let record = result?;

                let name = String::from_utf8(record.name().into())
                    .map_err(|e| anyhow::anyhow!("invalid utf8 in record name: {}", e))?;
                let seq = record.sequence();
                let seq_str = std::str::from_utf8(seq.as_ref())
                    .map_err(|e| anyhow::anyhow!("invalid utf8 in sequence: {}", e))?;

                let filename = pgr::libs::io::sanitize_filename(&name);
                write_record_to_fh(outdir, &mut fh_of, &filename, &name, seq_str, &mut out)?;
            }
        }
    } else if mode == "about" {
        let opt_count = if args.contains_id("chunk_size") {
            *args.get_one::<usize>("chunk_size").unwrap()
        } else {
            usize::MAX
        };
        let is_even = args.get_flag("even");
        let opt_max_part = *args.get_one::<usize>("max_part").unwrap();
        anyhow::ensure!(opt_max_part > 0, "--max-part must be positive");

        let mut chunker =
            pgr::libs::fasta::chunk::SizeChunker::new(opt_count, is_even, opt_max_part);
        let part_width = (opt_max_part.checked_ilog10().unwrap_or(0) + 1) as usize;

        'outer: for infile in args.get_many::<String>("infiles").unwrap() {
            let mut fa_in = pgr::libs::fmt::fa::reader(infile)
                .with_context(|| format!("Failed to open reader for {}", infile))?;

            for result in fa_in.records() {
                if chunker.max_files_exceeded() {
                    break 'outer;
                }

                // obtain record or fail with error
                let record = result?;

                let name = String::from_utf8(record.name().into())
                    .map_err(|e| anyhow::anyhow!("invalid utf8 in record name: {}", e))?;

                let seq = record.sequence();
                let seq_str = std::str::from_utf8(seq.as_ref())
                    .map_err(|e| anyhow::anyhow!("invalid utf8 in sequence: {}", e))?;

                let filename = format!("{:0width$}", chunker.file_index(), width = part_width);
                write_record_to_fh(outdir, &mut fh_of, &filename, &name, seq_str, &mut out)?;
                chunker.advance(seq.len());
            } // record
        } // file
    }

    // Explicitly flush all file handles to catch errors on close (e.g. disk full)
    for fh in fh_of.values_mut() {
        fh.flush()?;
    }
    out.flush()?;

    Ok(())
}

fn gen_fh(
    outdir: &str,
    fh_of: &mut BTreeMap<String, BufWriter<std::fs::File>>,
    filename: &str,
) -> anyhow::Result<()> {
    if !fh_of.contains_key(filename) {
        let path = std::path::Path::new(outdir).join(format!("{}.fa", filename));
        let file = std::fs::OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path)?;
        fh_of.insert(filename.to_owned(), BufWriter::new(file));
    }
    Ok(())
}

// Write a record either to stdout or to a file handle in fh_of.
fn write_record_to_fh(
    outdir: &str,
    fh_of: &mut BTreeMap<String, BufWriter<std::fs::File>>,
    filename: &str,
    name: &str,
    seq_str: &str,
    stdout_lock: &mut impl Write,
) -> anyhow::Result<()> {
    if outdir == "stdout" {
        write!(stdout_lock, ">{}\n{}\n", name, seq_str)?;
    } else {
        gen_fh(outdir, fh_of, filename)?;
        let fh = fh_of
            .get_mut(filename)
            .ok_or_else(|| anyhow::anyhow!("file handle not found for: {}", filename))?;
        write!(fh, ">{}\n{}\n", name, seq_str)?;
    }
    Ok(())
}
