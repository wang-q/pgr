use anyhow::Result;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::chain::ChainReader;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::fs::{self, File};
use std::hash::{Hash, Hasher};
use std::io::{BufReader, BufWriter};
use std::path::Path;

pub fn make_subcommand() -> Command {
    Command::new("split")
        .about("Split chains up by target or query sequence")
        .arg(Arg::new("out_dir").required(true).help("Output directory"))
        .arg(
            Arg::new("chains")
                .required(true)
                .num_args(1..)
                .help("Input chain file(s)"),
        )
        .arg(
            Arg::new("q")
                .short('q')
                .action(clap::ArgAction::SetTrue)
                .help("Split on query (default is on target)"),
        )
        .arg(
            Arg::new("lump")
                .long("lump")
                .value_parser(clap::value_parser!(usize))
                .help("Lump together so have only N split files"),
        )
}

fn lump_name(name: &str, lump: usize) -> String {
    // Look for integer part of name
    let mut s = name;
    while let Some(idx) = s.find(|c: char| c.is_ascii_digit()) {
        s = &s[idx..];
        let end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
        let digits = &s[..end];
        if let Ok(val) = digits.parse::<usize>() {
            return format!("{:03}", val % lump);
        }
        s = &s[end..];
    }

    // If no digits found, hash it
    let mut hasher = DefaultHasher::new();
    name.hash(&mut hasher);
    let hash = hasher.finish();
    format!("{:03}", (hash as usize) % lump)
}

pub fn execute(args: &ArgMatches) -> Result<()> {
    let out_dir = args.get_one::<String>("out_dir").unwrap();
    let chain_files: Vec<_> = args.get_many::<String>("chains").unwrap().collect();
    let split_on_q = args.get_flag("q");
    let lump = args.get_one::<usize>("lump").copied();

    fs::create_dir_all(out_dir)?;

    // Cache open file handles
    let mut file_cache: HashMap<String, BufWriter<File>> = HashMap::new();

    for file_path in chain_files {
        let f = File::open(file_path)?;
        let mut reader = ChainReader::new(BufReader::new(f));

        while let Some(res) = reader.next() {
            let chain = res?;

            let raw_name = if split_on_q {
                &chain.header.q_name
            } else {
                &chain.header.t_name
            };

            let name = if let Some(l) = lump {
                lump_name(raw_name, l)
            } else {
                raw_name.clone()
            };

            let writer = file_cache.entry(name.clone()).or_insert_with(|| {
                let path = Path::new(out_dir).join(format!("{}.chain", name));
                // We use create (truncate) if it's the first time we see this name in this run?
                // Wait, if we process multiple input files, we might append to the same output file.
                // The C code uses `hashFindVal` to check if it's open.
                // If not open, it does `mustOpen(path, "a")`.
                // But it also deletes the file if it exists? No, it says:
                // `safef(cmd,sizeof(cmd), "cat %s | sort -u > %s", tpath, path);` NO, that is for metadata!

                // Re-reading C code:
                // if ((f = hashFindVal(hash, name)) == NULL) {
                //    ...
                //    f = mustOpen(path, "a");
                //    hashAdd(hash, name, f);
                // }
                // So it always appends.
                // BUT, if we run the command twice, we don't want to append to previous run's output?
                // Usually these tools assume clean output directory or overwrite.
                // The C code `mustOpen(path, "a")` implies appending.
                // However, `makeDir(outDir)` might fail if exists? `makeDir` usually is `mkdir -p`.

                // Let's stick to "append" logic for now, but usually for a new run we might want to truncate.
                // Since we maintain a cache in memory, "first time in this process" = create/truncate?
                // If we process file1 then file2, and both have chr1.
                // Processing file1: chr1 -> open(create).
                // Processing file2: chr1 -> it is already in cache -> append.

                // But what if file1 and file2 both have chr1, but we closed the file handle (e.g. LRU)?
                // Then we re-open. If we used "create" every time we open, we would overwrite file1's chr1 data.
                // So we MUST use "append".
                // AND we should probably delete the file if it exists BEFORE we start writing to it for the first time?
                // The C code doesn't seem to delete files. It just appends.
                // So if you run it twice, you get duplicate data. This is typical for some UCSC tools.

                // To be safe and more user friendly, maybe we should try to be smart?
                // But for exact port, "append" is the way.
                // Wait, Rust `File::create` truncates. `File::options().append(true).open()` appends.
                // If I use `File::create` inside `or_insert_with`, it will truncate ONLY when we first see this name in this process.
                // Which is EXACTLY what we want for a single run!
                // Unless... we have so many files that we drop the handle and re-open it?
                // If we don't implement LRU and just keep all open (simple HashMap), then `File::create` is correct.
                // It truncates on first open (in this process), and keeps it open for subsequent writes.

                // So, using `File::create` is correct for the simple HashMap approach.
                // It ensures that for this run, we start fresh for each output file.

                let file = File::create(path).expect("Failed to create output file");
                BufWriter::new(file)
            });

            chain.write(writer)?;
        }
    }

    Ok(())
}
