use clap::*;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::fs::File;
use pgr::libs::ms::{parse_header, read_next_sample, perturb_positions, SimpleRng, system_seed};
use pgr::libs::ms2dna::{build_anc_seq, map_positions as map_pos, build_mut_seq, write_fasta};

pub fn make_subcommand() -> Command {
    Command::new("ms2dna")
        .about("Convert ms output haplotypes (0/1) to DNA sequences (FASTA)")
        .arg(
            Arg::new("gc")
                .long("gc")
                .short('g')
                .num_args(1)
                .default_value("0.5")
                .value_parser(value_parser!(f64))
                .help("GC content ratio in ancestral sequence (0..1)"),
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .short('s')
                .num_args(1)
                .value_parser(value_parser!(u64))
                .help("Random seed; default uses system time and PID"),
        )
        .arg(
            Arg::new("no_perturb")
                .long("no-perturb")
                .action(clap::ArgAction::SetTrue)
                .help("Disable positions micro-perturbation"),
        )
        .arg(
            Arg::new("verbose")
                .long("verbose")
                .short('v')
                .action(clap::ArgAction::SetTrue)
                .help("Print runtime information (paths and inputs)"),
        )
        .arg(
            Arg::new("doc")
                .long("doc")
                .action(clap::ArgAction::SetTrue)
                .help("Print the full documentation (markdown)"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
        .arg(
            Arg::new("files")
                .num_args(0..)
                .help("Input files with ms output; reads stdin when omitted"),
        )
        .after_help(
            r###"
Examples:
  # Pipe ms output to pgr ms2dna
  ms 10 1 -t 5 -r 0 1000 | pgr ms2dna --gc 0.5 > out.fa

  # Read from file and write to output
  pgr ms2dna input.ms -o out.fa --seed 12345

  # Disable position perturbation (keep original ms positions)
  pgr ms2dna input.ms --no-perturb

Output Format:
  FASTA format with single-line sequences.
  Headers: >[Lx_][Px_]Sx
    Lx: Batch/Replicate index (if multiple)
    Px: Population index (if multiple)
    Sx: Sample index
"###,
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    if args.get_flag("doc") {
        println!("{}", include_str!("../../docs/ms2dna.md"));
        return Ok(());
    }
    let outfile = args.get_one::<String>("outfile").unwrap();
    let gc = *args.get_one::<f64>("gc").unwrap_or(&0.5);
    let seed = args.get_one::<u64>("seed").copied();
    let no_perturb = args.get_flag("no_perturb");
    let verbose = args.get_flag("verbose");

    //----------------------------
    // Paths
    //----------------------------
    let curdir = std::env::current_dir()?;
    let pgr = std::env::current_exe()?.display().to_string();

    if verbose {
        println!("==> Paths");
        println!("    \"pgr\"     = {}", pgr);
        println!("    \"curdir\"  = {:?}", curdir);
    }

    if verbose {
        println!("==> Inputs");
    }
    let files: Vec<String> = args
        .get_many::<String>("files")
        .map(|vals| vals.map(|s| s.to_string()).collect())
        .unwrap_or_default();
    let abs_files: Vec<String> = files
        .iter()
        .map(|f| intspan::absolute_path(f).unwrap().display().to_string())
        .collect();
    if verbose {
        if abs_files.is_empty() {
            println!("    [stdin]");
        } else {
            println!("    files = {:?}", abs_files);
        }
    }

    let seed_final = seed.unwrap_or(system_seed());
    if verbose {
        println!("==> Seed");
        println!("    using = {}", seed_final);
    }

    let abs_outfile = if outfile == "stdout" {
        outfile.to_string()
    } else {
        intspan::absolute_path(outfile)?.display().to_string()
    };

    // Writer
    let mut writer: Box<dyn Write> = if abs_outfile == "stdout" {
        Box::new(std::io::stdout())
    } else {
        Box::new(intspan::writer(&abs_outfile))
    };

    // Process inputs (stdin or files)
    if abs_files.is_empty() {
        let stdin = std::io::stdin();
        convert_stream(BufReader::new(stdin.lock()), gc, Some(seed_final), &mut writer, no_perturb)?;
    } else {
        for path in abs_files {
            let fp = File::open(&path)?;
            convert_stream(BufReader::new(fp), gc, Some(seed_final), &mut writer, no_perturb)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_convert_stream_warning_and_output() {
        // Header: nsam=2, howmany=1, nsite=2
        // Sample: segsites=3 (> nsite) to trigger warning; haplotypes length=3
        let input = "\
ms 2 1 -r 0 2
//
segsites: 3
positions: 0.1000 0.5000 0.8000
010
001
";
        let mut out = Vec::new();
        let reader = BufReader::new(input.as_bytes());
        convert_stream(reader, 0.5, Some(123), &mut out, true).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert!(lines[0].starts_with("#WARNING: number of segregating sites"));
        assert!(lines[1].starts_with("#Hint: input may come from macs"));
        // Next should be headers and sequences for two samples
        assert!(lines[2].starts_with('>'));
        assert_eq!(lines[3].len(), 2);
        assert!(lines[4].starts_with('>'));
        assert_eq!(lines[5].len(), 2);
    }
}
fn convert_stream<R: BufRead>(
    mut reader: R,
    gc: f64,
    seed: Option<u64>,
    writer: &mut dyn Write,
    no_perturb: bool,
) -> anyhow::Result<()> {
    let mut header = String::new();
    reader.read_line(&mut header)?;
    if header.trim().is_empty() {
        return Ok(());
    }
    let hdr = parse_header(&header)?;
    let nsam = hdr.nsam;
    let howmany = hdr.howmany;
    let nsite = hdr.nsite;
    let npop = hdr.npop;
    let sample_sizes = hdr.sample_sizes;
    if nsite == 0 {
        anyhow::bail!("ERROR [ms2dna]: please use ms with the -r switch (nsite missing).");
    }

    let mut sample_counter = 0usize;
    let seed_final = seed.unwrap_or(system_seed());
    let mut rng = SimpleRng::new(seed_final);
    while let Some(sample) = read_next_sample(&mut reader, nsam)? {
            let segsites = sample.segsites;
            let mut positions = sample.positions;
            let haplotypes = sample.haplotypes;
            // Build sequences
            let seq_anc = build_anc_seq(gc, nsite, &mut rng);
            if segsites > 0 && !no_perturb {
                perturb_positions(&mut positions, &mut rng);
            }
            if segsites > nsite {
                writeln!(
                    writer,
                    "#WARNING: number of segregating sites ({}) > number of mutable sites ({})",
                    segsites, nsite
                )?;
                writeln!(writer, "#Hint: input may come from macs; ensure positions/nsite are compatible")?;
            }
            let map = map_pos(&positions, nsite, &mut rng);
            let seq_mut = build_mut_seq(&seq_anc, &map, gc, &mut rng, nsite);
            sample_counter += 1;
            write_fasta(
                writer,
                nsam,
                nsite,
                &map,
                &seq_anc,
                &seq_mut,
                &haplotypes,
                howmany,
                npop,
                sample_sizes.as_deref(),
                sample_counter,
            )?;
    }
    Ok(())
}
