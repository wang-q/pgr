use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("window")
        .about("Splits sequences into overlapping windows")
        .after_help(
            r###"
This command splits sequences in a FASTA file into overlapping windows.

Header format:
    >seq_name:start-end

Notes:
* Coordinates are 1-based, inclusive.
* Windows containing only Ns are skipped.
* Output sequences are unwrapped (single line).

Coverage & Overlap:
* Theoretical Coverage = Window Length / Step Size.
* Examples:
  - -l 200 -s 100: 2x coverage (50% overlap).
  - -l 200 -s 200: 1x coverage (no overlap).
  - -l 200 -s 10:  20x coverage (95% overlap).

Splitting & Shuffling:
* --chunk N: Splits output into files with N records each (e.g., output.001.fa).
* --shuffle: Randomizes output records.
  - With --chunk: Buffers N records, shuffles, writes to file, clears buffer (Low memory).
  - Without --chunk: Buffers ALL records, shuffles, writes to single file (High memory).
* --chunk cannot be used with stdout.

Examples:
1. Split into 200bp windows with 100bp step:
   pgr fa window input.fa -l 200 -s 100

2. Split large file into chunks of 1M records with shuffling:
   pgr fa window input.fa --chunk 1000000 --shuffle -o split.fa

3. Use default settings (200bp window, 100bp step):
   pgr fa window input.fa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input FASTA file to process"),
        )
        .arg(
            Arg::new("len")
                .long("len")
                .short('l')
                .value_parser(value_parser!(usize))
                .default_value("200")
                .help("Window length"),
        )
        .arg(
            Arg::new("step")
                .long("step")
                .short('s')
                .value_parser(value_parser!(usize))
                .default_value("100")
                .help("Step size"),
        )
        .arg(
            Arg::new("shuffle")
                .long("shuffle")
                .action(ArgAction::SetTrue)
                .help("Shuffle the output records (uses more memory)"),
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .value_parser(value_parser!(u64))
                .default_value("42")
                .help("Random seed for shuffling"),
        )
        .arg(
            Arg::new("chunk")
                .long("chunk")
                .value_parser(value_parser!(usize))
                .help("Split output into chunks of N records"),
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
    let infile = args.get_one::<String>("infile").unwrap();
    let len = *args.get_one::<usize>("len").unwrap();
    let step = *args.get_one::<usize>("step").unwrap();
    let shuffle = args.get_flag("shuffle");
    let seed = *args.get_one::<u64>("seed").unwrap();
    let chunk_size = args.get_one::<usize>("chunk").copied();
    let outfile = args.get_one::<String>("outfile").unwrap();

    // Check for conflict: chunk + stdout is not allowed because we can't split stdout into files
    if chunk_size.is_some() && outfile == "stdout" {
        return Err(anyhow::anyhow!("Cannot use --chunk with stdout output"));
    }

    let reader = pgr::reader(infile);
    let mut fa_in = noodles_fasta::io::Reader::new(reader);

    // Helper to create writer for a specific part
    let create_writer = |part: usize| -> Box<dyn std::io::Write> {
        if outfile == "stdout" {
            pgr::writer("stdout")
        } else {
            let path = std::path::Path::new(outfile);
            let file_stem = path.file_stem().unwrap().to_str().unwrap();
            let extension = path.extension().unwrap_or_default().to_str().unwrap();
            
            // Simple file splitting: input.fa -> input.001.fa
            // No special handling for .gz (output is usually uncompressed text for downstream tools or user specifies output name)
            let (stem, ext) = (file_stem, extension.to_string());

            let ext_str = if ext.is_empty() { String::new() } else { format!(".{}", ext) };
            let new_filename = format!("{}.{:03}{}", stem, part, ext_str);
            let new_path = path.with_file_name(new_filename);
            pgr::writer(new_path.to_str().unwrap())
        }
    };

    let mut current_part = 1;
    let mut record_count = 0;
    
    // Logic:
    // 1. Shuffle ON:
    //    - Chunked: Buffer one chunk -> Shuffle -> Write to new file -> Clear buffer.
    //    - No chunk: Buffer ALL -> Shuffle -> Write to single file.
    // 2. Shuffle OFF:
    //    - Chunked: Stream to file -> Switch file when limit reached.
    //    - No chunk: Stream to single file.

    let mut fa_out: Option<noodles_fasta::io::Writer<Box<dyn std::io::Write>>> = None;
    
    // Initialize global writer if not chunking.
    if chunk_size.is_none() {
        let writer = pgr::writer(outfile);
        fa_out = Some(noodles_fasta::io::writer::Builder::default()
            .set_line_base_count(usize::MAX)
            .build_from_writer(writer));
    } else if !shuffle {
        // If chunking without shuffle, init first writer
        let writer = create_writer(current_part);
        fa_out = Some(noodles_fasta::io::writer::Builder::default()
            .set_line_base_count(usize::MAX)
            .build_from_writer(writer));
    }

    // Reuse a single buffer to avoid reallocation if not needed, but for shuffle we accumulate.
    // For non-shuffle chunking, we don't need a large buffer.
    let mut records_buffer: Vec<noodles_fasta::Record> = Vec::new();

    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();
        let seq_len = seq.len();

        for start in (0..seq_len).step_by(step) {
            let end = std::cmp::min(start + len, seq_len);
            if start >= end { continue; }

            let start_pos = noodles_core::Position::new(start + 1).unwrap();
            let end_pos = noodles_core::Position::new(end).unwrap();
            let slice = seq.slice(start_pos..=end_pos).unwrap();

            if slice.as_ref().iter().all(|&b| pgr::libs::nt::is_n(b)) { continue; }

            let new_name = format!("{}:{}-{}", name, start + 1, end);
            let definition = noodles_fasta::record::Definition::new(new_name, None);
            let new_record = noodles_fasta::Record::new(definition, slice);
            
            if shuffle {
                records_buffer.push(new_record);
                
                // If chunk limit reached, flush buffer
                if let Some(limit) = chunk_size {
                    if records_buffer.len() >= limit {
                        use rand::seq::SliceRandom;
                        use rand::SeedableRng;
                        // Use deterministic seed derived from base seed + chunk index
                        let chunk_seed = seed + (current_part as u64); 
                        let mut rng = rand::rngs::StdRng::seed_from_u64(chunk_seed);
                        records_buffer.shuffle(&mut rng);
                        
                        // Write to current part file
                        let writer = create_writer(current_part);
                        let mut chunk_out = noodles_fasta::io::writer::Builder::default()
                            .set_line_base_count(usize::MAX)
                            .build_from_writer(writer);
                            
                        for r in &records_buffer {
                            chunk_out.write_record(r)?;
                        }
                        
                        records_buffer.clear();
                        current_part += 1;
                    }
                }
            } else {
                // No shuffle
                if let Some(limit) = chunk_size {
                    if record_count >= limit {
                        current_part += 1;
                        record_count = 0;
                        let writer = create_writer(current_part);
                        fa_out = Some(noodles_fasta::io::writer::Builder::default()
                            .set_line_base_count(usize::MAX)
                            .build_from_writer(writer));
                    }
                }
                
                if let Some(ref mut writer) = fa_out {
                    writer.write_record(&new_record)?;
                    record_count += 1;
                }
            }
        }
    }

    // Flush remaining buffer (Shuffle case)
    if shuffle && !records_buffer.is_empty() {
        use rand::seq::SliceRandom;
        use rand::SeedableRng;
        let chunk_seed = seed + (current_part as u64);
        let mut rng = rand::rngs::StdRng::seed_from_u64(chunk_seed);
        records_buffer.shuffle(&mut rng);

        // Flush remaining records.
        // If chunking, this goes to a new chunk file.
        // If not chunking, this goes to the single global file.
        
        let mut final_out = if chunk_size.is_some() {
             let writer = create_writer(current_part);
             noodles_fasta::io::writer::Builder::default()
                .set_line_base_count(usize::MAX)
                .build_from_writer(writer)
        } else {
             if let Some(writer) = fa_out.take() {
                 writer
             } else {
                 // Fallback (should not be reached if logic holds)
                 let writer = pgr::writer(outfile);
                 noodles_fasta::io::writer::Builder::default()
                    .set_line_base_count(usize::MAX)
                    .build_from_writer(writer)
             }
        };
        
        for record in records_buffer {
            final_out.write_record(&record)?;
        }
    }

    Ok(())
}
