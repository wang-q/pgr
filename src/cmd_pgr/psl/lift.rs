use clap::{Arg, ArgMatches, Command};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::str::FromStr;

use pgr::libs::psl::Psl;

pub fn make_subcommand() -> Command {
    Command::new("lift")
        .about("Lift PSL coordinates from query fragments (e.g., chr1:100-200) to genomic coordinates")
        .arg(
            Arg::new("in_psl")
                .index(1)
                .help("Input PSL file (default: stdin)"),
        )
        .arg(
            Arg::new("out_psl")
                .index(2)
                .help("Output PSL file (default: stdout)"),
        )
        .arg(
            Arg::new("sizes")
                .long("sizes")
                .short('s')
                .num_args(1)
                .help("File containing chromosome sizes (name <tab> size). Required for correct negative strand lifting and q_size updates."),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = args.get_one::<String>("in_psl").map(|s| s.as_str());
    let output = args.get_one::<String>("out_psl").map(|s| s.as_str());
    let sizes_file = args.get_one::<String>("sizes").map(|s| s.as_str());

    // Load sizes if provided
    let mut sizes_map = HashMap::new();
    if let Some(path) = sizes_file {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                if let Ok(size) = parts[1].parse::<u32>() {
                    sizes_map.insert(parts[0].to_string(), size);
                }
            }
        }
    }

    let reader: Box<dyn BufRead> = match input {
        Some(path) => Box::new(BufReader::new(File::open(path)?)),
        None => Box::new(BufReader::new(io::stdin())),
    };

    let writer: Box<dyn Write> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut writer = writer;

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            // Write comments/headers as is?
            // PSL usually doesn't have comments except headers.
            // But pgr tools might preserve them.
            // If it's a header line (doesn't parse as PSL), just write it?
            // Psl::from_str checks for column count.
            // Let's try to parse, if fails, print as is.
            if Psl::from_str(&line).is_err() {
                 writeln!(writer, "{}", line)?;
                 continue;
            }
        }

        let mut psl: Psl = match line.parse() {
            Ok(p) => p,
            Err(_) => {
                // If parsing fails (e.g. header), write as is
                writeln!(writer, "{}", line)?;
                continue;
            }
        };

        // Try to parse q_name as chr:start-end
        // Look for last colon to handle names with colons?
        // Usually format is Name:Start-End
        if let Some(colon_idx) = psl.q_name.rfind(':') {
            let range_part = &psl.q_name[colon_idx + 1..];
            let name_part = &psl.q_name[..colon_idx];

            if let Some(hyphen_idx) = range_part.find('-') {
                let start_str = &range_part[..hyphen_idx];
                // let end_str = &range_part[hyphen_idx + 1..]; // We don't strictly need end for offset

                if let Ok(start) = start_str.parse::<u32>() {
                    // start is 1-based
                    let offset: u32 = start.saturating_sub(1);
                    
                    let old_q_size = psl.q_size;
                    let is_neg = psl.strand.starts_with('-');

                    // Update name
                    psl.q_name = name_part.to_string();

                    // Update size if available
                    if let Some(&real_size) = sizes_map.get(&psl.q_name) {
                        psl.q_size = real_size;
                    }

                    // Calculate coordinate shift
                    // We perform lifting on (+) strand coordinates.
                    
                    if is_neg {
                        // In PSL, if strand is '-', qStart is from end of qSize.
                        
                        if let Some(real_size) = sizes_map.get(&psl.q_name) {
                            // We have sizes, we can do it correctly
                            // new_qEnd_plus = old_q_size - qStart + offset
                            // new_qStart = new_q_size - new_qEnd_plus
                            //            = new_q_size - (old_q_size - qStart + offset)
                            
                            psl.q_start = *real_size as i32 - (old_q_size as i32 - psl.q_start + offset as i32);
                            psl.q_end = *real_size as i32 - (old_q_size as i32 - psl.q_end + offset as i32);
                            
                            // Update blocks
                            for q_start in &mut psl.q_starts {
                                // q_start is start of block.
                                // block end is q_start + block_size.
                                // In '-' strand file, q_start is relative to RC start.
                                // qStart_plus_block = old_q_size - (qStart + size)
                                // new_qStart_plus_block = qStart_plus_block + offset
                                //                       = old_q_size - qStart - size + offset
                                // new_qEnd_plus_block = new_qStart_plus_block + size
                                //                     = old_q_size - qStart + offset
                                
                                // Back to '-' coords:
                                // new_qStart_file = new_q_size - new_qEnd_plus_block
                                //                 = new_q_size - (old_q_size - qStart + offset)
                                
                                *q_start = *real_size - (old_q_size - *q_start + offset as u32);
                            }
                        } else {
                            // No sizes map. We cannot compute correct '-' strand coordinates.
                            // Fallback: Treat as '+' strand but warn? 
                            // Or just apply offset blindly?
                            // If we apply offset blindly: qStart += offset.
                            // This means on RC(Genome), the match is shifted by offset.
                            // But RC(Genome) starts far away.
                            // This would be wrong.
                            
                            // If user doesn't provide sizes, maybe we should assume they want '+' strand output?
                            // But the alignment is '-' strand.
                            // Let's output to stderr warning and skip modification or do best effort?
                            // Best effort: Just add offset? No that's garbage.
                            // Best effort: Convert to '+' strand?
                            //   qStart_plus = old_q_size - qEnd + offset.
                            //   psl.strand = "+".to_string() + ...
                            //   And update qStarts.
                            // This preserves the genomic location correctly, but changes strand representation.
                            // This is probably the most useful behavior if sizes are missing.
                            
                            // Let's try to convert to '+' strand if sizes are missing.
                            // But we need to update t_starts/block order?
                            // No, t_starts are always increasing.
                            // If we switch query strand to '+', the blocks are still in same order on Target.
                            // But on Query(+), are they increasing?
                            // If Query(-) aligned to Target(+), then on Query(-), blocks are increasing (because t is increasing).
                            // On Query(+), blocks would be DECREASING.
                            // PSL requires blocks to be ordered by Target.
                            // So if we switch to Query(+), we still list blocks in Target order.
                            // So `qStarts` will be non-increasing?
                            // "The qStarts array is always listed in increasing order of tStarts."
                            // If the match is inverted, `qStarts` will indeed be decreasing?
                            // Wait, does PSL allow decreasing `qStarts`?
                            // BLAT PSL: "qStarts: ... list of starting positions ...".
                            // If strand is `-`, `qStarts` are on RC. They are increasing because Target is increasing and we walk along the alignment.
                            // If we report on `+` strand for an inverted alignment:
                            // The 5' end of Target matches 3' end of Query(+).
                            // As we move 5'->3' on Target, we move 3'->5' on Query(+).
                            // So `qStarts` (on +) would decrease.
                            // Does PSL allow this?
                            // "qStarts ... increasing order of tStarts".
                            // It doesn't explicitly forbid decreasing values, but usually PSL implies collinear blocks unless it's a "chain" or something.
                            // But `pgr` handles standard PSL.
                            // If `qStarts` are decreasing, it's not a standard linear alignment?
                            // Actually, if it's a single block, it's fine.
                            // If multiple blocks, and it's an inversion, then yes, qStarts decrease.
                            
                            // So, converting to '+' strand is safe-ish.
                            
                            // BUT, I'll stick to the "require sizes" for now to be safe, or just do the naive thing for '+' and warn for '-'.
                            // Actually, the user prompt said "Analysis...". I am implementing.
                            // I will add a warning if `-` strand and no sizes.
                            
                            eprintln!("Warning: '-' strand record found for {} but no sizes provided. Skipping lift for this record.", psl.q_name);
                        }
                    } else {
                        // Strand is '+'. Easy.
                        psl.q_start += offset as i32;
                        psl.q_end += offset as i32;
                        for q_start in &mut psl.q_starts {
                            *q_start += offset;
                        }
                    }
                }
            }
        }

        psl.write_to(&mut writer)?;
    }

    Ok(())
}
