use pgr::libs::loc::Input;
use clap::*;
use cmd_lib::*;
use rayon::prelude::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("prefilter")
        .about("Prefilter genome/metagenome assembly by amino acid minimizers")
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
* --len N: Minimum peptide length to consider (filters short ORFs)
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
   pgr pl prefilter assembly.fa refs.pep.fa --chunk 50000 --len 20

3. Use custom k-mer and window sizes:
   pgr pl prefilter assembly.fa refs.pep.fa -k 7 -w 2 --parallel 8

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input file containing genome/metagenome assembly"),
        )
        .arg(
            Arg::new("match")
                .required(true)
                .index(2)
                .help("Match file containing reference sequences for filtering"),
        )
        .arg(
            Arg::new("chunk")
                .long("chunk")
                .short('c')
                .num_args(1)
                .default_value("100000")
                .value_parser(value_parser!(usize))
                .help("Size of each chunk in bytes"),
        )
        .arg(
            Arg::new("len")
                .long("len")
                .num_args(1)
                .default_value("15")
                .value_parser(value_parser!(usize))
                .help("Minimum length of the amino acid sequence to consider"),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("7")
                .value_parser(value_parser!(usize))
                .help("K-mer size"),
        )
        .arg(
            Arg::new("window")
                .long("window")
                .short('w')
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Window size for minimizers"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Number of threads for parallel processing"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let match_file = args.get_one::<String>("match").unwrap();

    let opt_chunk = *args.get_one::<usize>("chunk").unwrap();
    let opt_len = *args.get_one::<usize>("len").unwrap();
    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_window = *args.get_one::<usize>("window").unwrap();

    // Set the number of threads for rayon
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt_parallel)
        .build_global()?;

    let is_bgzf = {
        let path = std::path::Path::new(infile);
        path.extension() == Some(std::ffi::OsStr::new("gz"))
    };

    //----------------------------
    // Open files
    //----------------------------
    let loc_file = format!("{}.loc", infile);
    if !std::path::Path::new(&loc_file).is_file() {
        pgr::libs::loc::create_loc(infile, &loc_file, is_bgzf)?;
    }

    // Split .loc file into chunks
    let chunks = split_loc_file(&loc_file, opt_chunk)?;

    let pgr = std::env::current_exe()?.display().to_string();

    chunks.par_iter().for_each_init(
        || {
            // Init reader for each chunk
            if is_bgzf {
                Input::Bgzf(
                    noodles_bgzf::io::indexed_reader::Builder::default()
                        .build_from_path(infile)
                        .unwrap(),
                )
            } else {
                Input::File(std::fs::File::open(std::path::Path::new(infile)).unwrap())
            }
        },
        |reader, (_, offset, size)| {
            let chunk = pgr::libs::loc::read_offset(reader, *offset, *size).unwrap();

            let mut temp_file = tempfile::NamedTempFile::new().unwrap();
            temp_file.write_all(&chunk).unwrap();
            let temp_path = temp_file.path().to_str().unwrap().to_string();

            run_cmd!(
                ${pgr} fa six-frame ${temp_path} --len ${opt_len} |
                    ${pgr} dist seq stdin ${match_file} -k ${opt_kmer} -w ${opt_window}
            )
            .unwrap();
        },
    );

    Ok(())
}

// Split .loc file into chunks
fn split_loc_file(loc_file: &str, chunk_size: usize) -> anyhow::Result<Vec<(String, u64, usize)>> {
    // Load .loc file
    let loc_of: indexmap::IndexMap<String, (u64, usize)> = pgr::libs::loc::load_loc(loc_file)?;

    let mut chunks: Vec<(String, u64, usize)> = Vec::new();
    let mut cur_size = 0;
    let mut cur_start_offset = 0;
    let mut cur_first_seq = String::new();

    // Iterate over each sequence in the .loc file
    for (seq_id, &(offset, size)) in &loc_of {
        // If the current chunk size exceeds the specified size,
        //   record the current chunk and start a new one
        if cur_size + size > chunk_size && !cur_first_seq.is_empty() {
            chunks.push((cur_first_seq.clone(), cur_start_offset, cur_size));
            cur_size = 0;
            cur_start_offset = offset;
            cur_first_seq = seq_id.clone();
        }

        // If it's the first sequence, record the start offset and sequence name
        if cur_size == 0 {
            cur_start_offset = offset;
            cur_first_seq = seq_id.clone();
        }

        // Update the current chunk size
        cur_size += size;
    }

    // Add the last chunk
    if !cur_first_seq.is_empty() {
        chunks.push((cur_first_seq, cur_start_offset, cur_size));
    }

    Ok(chunks)
}
