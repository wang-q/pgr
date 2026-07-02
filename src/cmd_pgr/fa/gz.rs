use anyhow::Context;
use clap::*;
use noodles_bgzf as bgzf;
use std::io::{Read, Write};

use pgr::libs::fmt::fa;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("gz")
        .about("Compressing a file using the BGZF format")
        .after_help(
            r###"
This command compresses a file using BGZF (Blocked Gzip Format).

Notes:
* Parallel compression with multiple threads
* Creates index file (.gzi) for random access
* Supports stdin as input (use 'stdin' as filename)
* Preserves original file
* Cannot compress already gzipped files
* Default thread count is 1
* Index creation is automatic

Output files:
* <infile>.gz: Compressed file
* <infile>.gz.gzi: Index file

Examples:
1. Compress a file with default settings, and the outfile is input.fa.gz:
   pgr fa gz input.fa

2. Multi-threaded compression:
   pgr fa gz input.fa -p 4

3. Set compression level (0-9, default -1):
   pgr fa gz input.fa --compress-level 9

4. Create index for existing .gz file (reindex):
   pgr fa gz input.fa.gz --reindex

5. From stdin with custom output:
   cat input.fa | pgr fa gz stdin -o output.fa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input FASTA file to compress"),
        )
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(
            Arg::new("compress_level")
                .long("compress-level")
                .value_parser(value_parser!(i32))
                .num_args(1)
                .default_value("-1")
                .help("Compression level (0-9, or -1 for default)"),
        )
        .arg(
            Arg::new("reindex")
                .long("reindex")
                .action(ArgAction::SetTrue)
                .help("Create BGZF index (.gzi) for an existing .gz file"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_optional())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();

    if args.get_flag("reindex") {
        if !std::path::Path::new(infile).exists() {
            return Err(anyhow::anyhow!("Input file not found: {}", infile));
        }
        if !pgr::is_bgzf(infile) {
            return Err(anyhow::anyhow!(
                "Input file is not a valid BGZF file: {}",
                infile
            ));
        }
        fa::build_gzi_index(infile)?;
        return Ok(());
    }

    let opt_parallel: std::num::NonZeroUsize =
        (*args.get_one::<usize>("parallel").unwrap()).try_into()?;
    let compress_level = *args.get_one::<i32>("compress_level").unwrap();

    let outfile = if args.contains_id("outfile") {
        crate::cmd_pgr::args::get_outfile(args).to_string()
    } else {
        format!("{}.gz", infile)
    };

    //----------------------------
    // Input
    //----------------------------
    let mut reader: Box<dyn std::io::BufRead> = if infile == "stdin" {
        // Use 64KB buffer (BGZF block size) to optimize read performance
        Box::new(std::io::BufReader::with_capacity(
            64 * 1024,
            std::io::stdin(),
        ))
    } else {
        let path = std::path::Path::new(infile);
        let file = std::fs::File::open(path)?;

        // Use 64KB buffer (BGZF block size) to optimize read performance
        Box::new(std::io::BufReader::with_capacity(64 * 1024, file))
    };

    let inner_writer = Box::new(std::io::BufWriter::new(
        std::fs::File::create(&outfile).with_context(|| format!("Failed to create: {outfile}"))?,
    ));

    let mut builder =
        bgzf::io::multithreaded_writer::Builder::default().set_worker_count(opt_parallel);

    if (0..=9).contains(&compress_level) {
        use noodles_bgzf::io::writer::CompressionLevel;
        builder =
            builder.set_compression_level(CompressionLevel::new(compress_level as u8).unwrap());
    }

    let mut writer = builder.build_from_writer(inner_writer);

    //----------------------------
    // Output
    //----------------------------
    // Manually read/write in 64KB chunks (BGZF block size) instead of std::io::copy.
    // std::io::copy uses a smaller default buffer (usually 8KB), which causes frequent small writes.
    // For MultithreadedWriter, this increases channel/lock overhead significantly,
    // negating the benefits of parallelism on small/medium files.
    let mut buf = vec![0u8; 64 * 1024];
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        writer.write_all(&buf[..n])?;
    }
    writer.finish()?;

    // Generate GZI index
    fa::build_gzi_index(&outfile)?;

    Ok(())
}
