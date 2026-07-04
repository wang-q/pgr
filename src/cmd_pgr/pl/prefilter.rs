use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use cmd_lib::run_cmd;
use rayon::prelude::*;
use std::io::{Read, Write};
use std::path::PathBuf;

/// Build the clap subcommand for prefilter.
pub fn make_subcommand() -> Command {
    Command::new("prefilter")
        .about("Prefilters genome/metagenome assembly by amino acid minimizers")
        .after_help(
            r###"
This command filters genome/metagenome sequences by comparing them against protein references using amino acid minimizers.
It processes input files in chunks for efficient memory usage and supports parallel processing.

Process:
1. Splits input genome file into chunks
2. Translates sequences in six frames
3. Calculates amino acid minimizers
4. Compares with reference protein sequences

Parameters:
* --chunk N: Process N bytes at a time (memory control)
* --min-len N: Minimum peptide length to consider (filters short ORFs)
* --kmer/-k N: K-mer size for minimizers
* --window/-w N: Window size for minimizers
* --parallel/-p N: Number of threads

Notes:
* Input file must be FASTA or BGZF-compressed FASTA (.gz)
* Reference file must be Protein FASTA
* Automatic index creation (.loc) if missing
* Cannot read from stdin (requires random access)
* Output format matches `pgr dist seq`

Examples:
1. Basic usage:
   pgr pl prefilter assembly.fa refs.pep.fa

2. Specify chunk size and minimum peptide length:
   pgr pl prefilter assembly.fa refs.pep.fa --chunk 50000 --min-len 20

3. Use custom k-mer and window sizes:
   pgr pl prefilter assembly.fa refs.pep.fa -k 7 -w 2 --parallel 8

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input file containing genome/metagenome assembly",
        ))
        .arg(
            Arg::new("reference")
                .required(true)
                .index(2)
                .help("Reference file containing sequences for filtering"),
        )
        .arg(crate::cmd_pgr::args::chunk_size_arg(
            Some("100000"),
            "Size of each chunk in bytes",
        ))
        .arg(crate::cmd_pgr::args::min_len_arg_with_default(
            "15",
            "Minimum length of the amino acid sequence to consider",
        ))
        .arg(crate::cmd_pgr::args::kmer_arg_with_default("7"))
        .arg(crate::cmd_pgr::args::window_arg())
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the prefilter command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let match_file = args.get_one::<String>("reference").unwrap();

    let opt_chunk = *args.get_one::<usize>("chunk_size").unwrap();
    let opt_len = *args.get_one::<usize>("min_len").unwrap();
    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_window = *args.get_one::<usize>("window").unwrap();

    // Set the number of threads for rayon
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt_parallel)
        .build_global()?;

    let is_bgzf = pgr::is_bgzf(infile);

    // Open files
    let loc_file = format!("{}.loc", infile);
    if !std::path::Path::new(&loc_file).is_file() {
        pgr::libs::loc::create_loc(infile, &loc_file, is_bgzf)?;
    }

    // Split .loc file into chunks
    let chunks = pgr::libs::loc::split_loc_file(&loc_file, opt_chunk)?;

    let pgr = pgr::libs::io::current_exe_string()?;

    // Each parallel task writes its sub-process stdout to a private temp file to
    // avoid interleaved output across rayon workers. Temp files are kept alive
    // until their path is consumed by the serial cat phase below.
    let results: Vec<anyhow::Result<PathBuf>> = chunks
        .par_iter()
        .map(|(_, offset, size)| {
            // Init reader for this chunk
            let mut reader = pgr::libs::loc::open_input(infile, is_bgzf)?;

            let chunk = pgr::libs::loc::read_offset(&mut reader, *offset, *size)?;

            let mut temp_input = tempfile::NamedTempFile::new()?;
            temp_input.write_all(&chunk)?;
            let temp_input_path = temp_input
                .path()
                .to_str()
                .ok_or_else(|| anyhow::anyhow!("temp file path is not valid UTF-8"))?
                .to_string();

            // Persist the input temp file so its path remains valid after the
            // handle is dropped; we clean it up manually after the pipeline.
            let temp_input_persist_path = temp_input.into_temp_path().keep()?;
            let temp_input_keep = temp_input_persist_path.clone();

            let temp_output = tempfile::Builder::new().tempfile_in(std::env::temp_dir())?;
            let temp_output_path = temp_output
                .path()
                .to_str()
                .ok_or_else(|| anyhow::anyhow!("temp file path is not valid UTF-8"))?
                .to_string();
            // Persist output temp file too; we read & delete it in the serial phase.
            let temp_output_persist_path = temp_output.into_temp_path().keep()?;

            run_cmd!(
                ${pgr} fa six-frame ${temp_input_path} --min-len ${opt_len} |
                    ${pgr} dist seq stdin ${match_file} -k ${opt_kmer} -w ${opt_window} > ${temp_output_path}
            )?;

            // Best-effort cleanup of the input temp file; failure is non-fatal.
            let _ = std::fs::remove_file(&temp_input_keep);

            Ok(temp_output_persist_path)
        })
        .collect();

    // Serial phase: stream each chunk's output in order to the writer.
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    for result in results {
        let path = result?;
        let mut f = std::fs::File::open(&path)
            .with_context(|| format!("Failed to open file {}", path.display()))?;
        let mut buf = Vec::new();
        f.read_to_end(&mut buf)?;
        writer.write_all(&buf)?;
        let _ = std::fs::remove_file(&path);
    }
    writer.flush()?;

    Ok(())
}
