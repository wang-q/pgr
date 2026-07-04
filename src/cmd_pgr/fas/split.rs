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

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            let first = match block.entries.first() {
                Some(e) => e,
                None => continue,
            };
            let filename = if is_chr {
                format!("{}.{}", first.range().name(), first.range().chr())
            } else {
                first.range().to_string()
            }
            .replace(['(', ')', ':'], "_")
            .replace("__", "_");

            for entry in &block.entries {
                let range = entry.range().clone();
                let seq = std::str::from_utf8(entry.seq())?;

                if outdir == "stdout" {
                    let header = if is_simple {
                        range.name().to_string()
                    } else {
                        range.to_string()
                    };
                    print!(">{}\n{}\n", header, seq);
                } else {
                    if !file_of.contains_key(&filename) {
                        let path = std::path::Path::new(outdir).join(filename.clone() + opt_suffix);
                        let file = std::fs::OpenOptions::new()
                            .create(true)
                            .write(true)
                            .truncate(true)
                            .open(path)?;
                        file_of.insert(filename.clone(), BufWriter::new(file));
                    }
                    let file = file_of
                        .get_mut(&filename)
                        .ok_or_else(|| anyhow::anyhow!("file not found: {}", filename))?;
                    write!(file, ">{}\n{}\n", range, seq)?;
                }
            }

            // end of a block
            if outdir == "stdout" {
                println!();
            } else {
                let file = file_of
                    .get_mut(&filename)
                    .ok_or_else(|| anyhow::anyhow!("file not found: {}", filename))?;
                writeln!(file)?;
            }
        }
    }

    // Explicitly flush all file handles to catch errors on close (e.g. disk full)
    for fh in file_of.values_mut() {
        fh.flush()?;
    }

    Ok(())
}
