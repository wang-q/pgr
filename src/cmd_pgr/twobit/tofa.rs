use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::twobit::TwoBitFile;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

pub fn make_subcommand() -> Command {
    Command::new("tofa")
        .about("Convert 2bit to FASTA")
        .after_help(
            r###"
Examples:
  # Convert entire 2bit file to FASTA
  pgr twobit tofa input.2bit -o output.fa

  # Extract single sequence
  pgr twobit tofa input.2bit --seq chr1 -o chr1.fa
  pgr twobit tofa input.2bit --seq chr1 --start 0 --end 100 -o chr1_head.fa

  # Extract sequences from list
  pgr twobit tofa input.2bit --seqList list.txt -o out.fa

  # No masking (all uppercase)
  pgr twobit tofa input.2bit --no-mask -o out.fa
"###,
        )
        .arg(
            Arg::new("input")
                .help("Input 2bit file")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output FASTA file")
                .default_value("stdout"),
        )
        .arg(
            Arg::new("seq")
                .long("seq")
                .value_name("NAME")
                .help("Restrict to this sequence"),
        )
        .arg(
            Arg::new("start")
                .long("start")
                .value_name("INT")
                .value_parser(clap::value_parser!(usize))
                .help("Start position (0-based)"),
        )
        .arg(
            Arg::new("end")
                .long("end")
                .value_name("INT")
                .value_parser(clap::value_parser!(usize))
                .help("End position (non-inclusive)"),
        )
        .arg(
            Arg::new("seq_list")
                .long("seqList")
                .value_name("FILE")
                .help("File containing list of sequence names (one per line)"),
        )
        .arg(
            Arg::new("no_mask")
                .long("no-mask")
                .action(ArgAction::SetTrue)
                .help("Convert sequence to all upper case"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("input").unwrap();
    let output_path = args.get_one::<String>("output").unwrap();
    let opt_seq = args.get_one::<String>("seq");
    let opt_start = args.get_one::<usize>("start").copied();
    let opt_end = args.get_one::<usize>("end").copied();
    let opt_seq_list = args.get_one::<String>("seq_list");
    let no_mask = args.get_flag("no_mask");

    let mut tb = TwoBitFile::open(input_path)?;
    let mut writer = intspan::writer(output_path);

    // Determine targets: Vec<(name, start, end)>
    let mut targets = Vec::new();

    if let Some(seq) = opt_seq {
        targets.push((seq.clone(), opt_start, opt_end));
    } else if let Some(list_path) = opt_seq_list {
        let file = File::open(list_path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            
            // Parse name[:start-end]
            if let Some(colon_idx) = line.find(':') {
                let name = line[..colon_idx].to_string();
                let range_part = &line[colon_idx+1..];
                if let Some(dash_idx) = range_part.find('-') {
                    let start_str = &range_part[..dash_idx];
                    let end_str = &range_part[dash_idx+1..];
                    let start = start_str.parse::<usize>().ok();
                    let end = end_str.parse::<usize>().ok();
                    targets.push((name, start, end));
                } else {
                    // Invalid range format, treat as name? Or error?
                    // UCSC allows just name, or name:start-end. 
                    // If no dash, maybe just name? But colon implies range.
                    // Let's assume if colon exists, we try to parse range.
                    targets.push((name, None, None)); 
                }
            } else {
                targets.push((line.to_string(), None, None));
            }
        }
    } else {
        // All sequences
        let names = tb.get_sequence_names();
        for name in names {
            targets.push((name, None, None));
        }
    }

    for (name, start, end) in targets {
        // If start/end provided, we might want to adjust them or validate?
        // read_sequence handles None by using default.
        
        let seq = tb.read_sequence(&name, start, end, no_mask)?;
        
        // Write FASTA
        // Header: >name (if range, maybe add range info? UCSC usually just >name:start-end)
        // UCSC twoBitToFa behavior:
        // If -seq is used: >name
        // If -seqList is used: >name:start-end (if range specified)
        // If whole file: >name
        
        // Let's stick to simple >name unless it's a subslice.
        // If start/end are specified, we should probably indicate it in header
        // to match UCSC or just be helpful.
        // But UCSC `twoBitToFa` output header depends on input.
        // If I do `twoBitToFa hg38.2bit -seq=chr1 -start=0 -end=10 out.fa`
        // Output header is `>chr1:0-10`.
        
        let header = if start.is_some() || end.is_some() {
            // We need to know the actual start/end used.
            // read_sequence doesn't return them.
            // But we can guess.
            // If start is None, it's 0.
            // If end is None, it's seq len. 
            // We don't have seq len easily here without querying index.
            // But `read_sequence` does it internally.
            
            // For now, let's construct header based on request.
            let s = start.unwrap_or(0);
            if let Some(e) = end {
                format!("{}:{}-{}", name, s, e)
            } else {
                // If end is None, we don't know the end without querying size.
                // We can query size from `tb`.
                // But `tb` interface for size is `sequence_offsets`? No, that gives offset.
                // We need to read the record header to get size.
                // Or just use the returned sequence length.
                let e = s + seq.len();
                format!("{}:{}-{}", name, s, e)
            }
        } else {
            name
        };

        writeln!(writer, ">{}", header)?;
        // Wrap lines? UCSC does 50 chars.
        // Let's just write line by line or wrap at 60/80.
        // `fas` usually 60 or 80.
        // Let's use 60.
        
        let mut idx = 0;
        let len = seq.len();
        while idx < len {
            let next_idx = (idx + 60).min(len);
            writeln!(writer, "{}", &seq[idx..next_idx])?;
            idx = next_idx;
        }
    }

    Ok(())
}
