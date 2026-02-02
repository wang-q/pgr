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
            Arg::new("bed")
                .long("bed")
                .value_name("FILE")
                .help("Grab sequences specified by input.bed"),
        )
        .arg(
            Arg::new("bed_pos")
                .long("bedPos")
                .action(ArgAction::SetTrue)
                .help("With -bed, use chrom:start-end as the fasta ID"),
        )
        .arg(
            Arg::new("no_mask")
                .long("no-mask")
                .action(ArgAction::SetTrue)
                .help("Convert sequence to all upper case"),
        )
}

struct Target {
    seq_name: String,
    start: Option<usize>,
    end: Option<usize>,
    header_name: String,
    is_rev_comp: bool,
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("input").unwrap();
    let output_path = args.get_one::<String>("output").unwrap();
    let opt_seq = args.get_one::<String>("seq");
    let opt_start = args.get_one::<usize>("start").copied();
    let opt_end = args.get_one::<usize>("end").copied();
    let opt_seq_list = args.get_one::<String>("seq_list");
    let opt_bed = args.get_one::<String>("bed");
    let bed_pos = args.get_flag("bed_pos");
    let no_mask = args.get_flag("no_mask");

    let mut tb = TwoBitFile::open(input_path)?;
    let mut writer = intspan::writer(output_path);

    // Determine targets
    let mut targets: Vec<Target> = Vec::new();

    if let Some(seq) = opt_seq {
        let header_name = if opt_start.is_some() || opt_end.is_some() {
             // Will be refined later based on actual end
             seq.clone() 
        } else {
             seq.clone()
        };
        
        targets.push(Target {
            seq_name: seq.clone(),
            start: opt_start,
            end: opt_end,
            header_name,
            is_rev_comp: false,
        });
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
                    // For seqList, UCSC behavior is usually to use the line as header?
                    // Or name:start-end.
                    // Let's use the full line spec as header name if range is present.
                    targets.push(Target {
                        seq_name: name,
                        start,
                        end,
                        header_name: line.to_string(),
                        is_rev_comp: false,
                    });
                } else {
                    targets.push(Target {
                        seq_name: name,
                        start: None,
                        end: None,
                        header_name: line.to_string(),
                        is_rev_comp: false,
                    });
                }
            } else {
                targets.push(Target {
                    seq_name: line.to_string(),
                    start: None,
                    end: None,
                    header_name: line.to_string(),
                    is_rev_comp: false,
                });
            }
        }
    } else if let Some(bed_path) = opt_bed {
        let file = File::open(bed_path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 3 {
                continue; // Skip invalid lines
            }
            
            let chrom = fields[0].to_string();
            let start = fields[1].parse::<usize>().ok();
            let end = fields[2].parse::<usize>().ok();
            let name = if fields.len() > 3 { fields[3].to_string() } else { chrom.clone() };
            let strand = if fields.len() > 5 { fields[5] } else { "+" };
            
            let is_rev_comp = strand == "-";
            
            let header_name = if bed_pos {
                if let (Some(s), Some(e)) = (start, end) {
                    format!("{}:{}-{}", chrom, s, e)
                } else {
                    chrom.clone()
                }
            } else {
                name
            };

            targets.push(Target {
                seq_name: chrom,
                start,
                end,
                header_name,
                is_rev_comp,
            });
        }
    } else {
        // All sequences
        let names = tb.get_sequence_names();
        for name in names {
            targets.push(Target {
                seq_name: name.clone(),
                start: None,
                end: None,
                header_name: name,
                is_rev_comp: false,
            });
        }
    }

    for target in targets {
        let mut seq = tb.read_sequence(&target.seq_name, target.start, target.end, no_mask)?;
        
        if target.is_rev_comp {
            seq = reverse_complement(&seq);
        }

        // Refine header if needed (for --seq case where we constructed it lazily)
        let header = if opt_seq.is_some() && (target.start.is_some() || target.end.is_some()) {
            let s = target.start.unwrap_or(0);
            let e = s + seq.len(); // Approximate end based on read length
             format!("{}:{}-{}", target.seq_name, s, e)
        } else {
            target.header_name
        };

        writeln!(writer, ">{}", header)?;
        
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

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T', 'a' => 't',
            'C' => 'G', 'c' => 'g',
            'G' => 'C', 'g' => 'c',
            'T' => 'A', 't' => 'a',
            'U' => 'A', 'u' => 'a',
            'N' => 'N', 'n' => 'n',
            _ => c,
        })
        .collect()
}
