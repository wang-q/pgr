use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};

/// Build the clap subcommand for split.
pub fn make_subcommand() -> Command {
    Command::new("split")
        .about("Splits block FA files into per-alignment or per-chromosome FA files")
        .after_help(
            r###"
Splits block FA files into per-alignment or per-chromosome FA files.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* By default, each alignment block is written to a separate file
* Use `--chr` to split files by chromosome
* Use `--simple` to simplify headers by keeping only species names

Examples:
1. Split block FA files into per-alignment files:
   pgr fas split tests/fas/example.fas -o output_dir

2. Split block FA files into per-chromosome files:
   pgr fas split tests/fas/example.fas -o output_dir --chr

3. Simplify headers in output files:
   pgr fas split tests/fas/example.fas -o output_dir --simple

4. Use a custom suffix for output files:
   pgr fas split tests/fas/example.fas -o output_dir --suffix .fa

5. Output to stdout:
   pgr fas split tests/fas/example.fas

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::suffix_arg(".fas"))
        .arg(
            Arg::new("chr")
                .long("chr")
                .action(ArgAction::SetTrue)
                .help("Split files by chromosomes"),
        )
        .arg(
            Arg::new("simple")
                .long("simple")
                .action(ArgAction::SetTrue)
                .help("Simplify headers by keeping only species names"),
        )
        .arg(crate::cmd_pgr::args::outdir_arg())
}

/// Execute the split command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outdir = args.get_one::<String>("outdir").unwrap();
    if outdir != "stdout" {
        std::fs::create_dir_all(outdir)?;
    }

    let opt_suffix = args.get_one::<String>("suffix").unwrap();
    let is_chr = args.get_flag("chr");
    let is_simple = args.get_flag("simple");

    let mut file_of: BTreeMap<String, BufWriter<std::fs::File>> = BTreeMap::new();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for (block_idx, block_result) in
            pgr::libs::fmt::fas::iter_fas_blocks(&mut reader).enumerate()
        {
            let block = block_result
                .with_context(|| format!("read block {} from {}", block_idx, infile))?;
            let filename = match pgr::libs::fmt::fas::split_block_key(&block, is_chr) {
                Some(k) => pgr::libs::io::sanitize_filename(&k),
                None => continue,
            };
            let block_str = pgr::libs::fmt::fas::format_split_block(&block, is_simple)?;

            if outdir == "stdout" {
                writeln!(out, "{}", block_str)?;
            } else {
                let file = if let Some(fh) = file_of.get_mut(&filename) {
                    fh
                } else {
                    let path =
                        std::path::Path::new(outdir).join(format!("{}{}", filename, opt_suffix));
                    let file = std::fs::OpenOptions::new()
                        .create(true)
                        .write(true)
                        .truncate(true)
                        .open(path)?;
                    file_of.entry(filename).or_insert(BufWriter::new(file))
                };
                writeln!(file, "{}", block_str)?;
            }
        }
    }

    // Explicitly flush all file handles to catch errors on close (e.g. disk full)
    for fh in file_of.values_mut() {
        fh.flush()?;
    }
    out.flush()?;

    Ok(())
}
