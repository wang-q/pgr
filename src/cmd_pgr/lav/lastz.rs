use clap::*;
use rayon::prelude::*;
use std::io::Write;
use tempfile::NamedTempFile;

const MATRIX_DEFAULT: &str = r#"
   A    C    G    T
  91 -114  -31 -123
-114  100 -125  -31
 -31 -125  100 -114
-123  -31 -114   91
"#;

const MATRIX_DISTANT: &str = r#"
   A    C    G    T
  91  -90  -25 -100
 -90  100 -100  -25
 -25 -100  100  -90
-100  -25  -90   91
"#;

const MATRIX_SIMILAR: &str = r#"
   A    C    G    T
 100 -300 -150 -300
-300  100 -300 -150
-150 -300  100 -300
-300 -150 -300  100
"#;

#[allow(dead_code)]
const MATRIX_SIMILAR2: &str = r#"
   A    C    G    T
  90 -330 -236 -356
-330  100 -318 -236
-236 -318  100 -330
-356 -236 -330   90
"#;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("lastz")
        .about("Wrapper for lastz alignment (Cactus style)")
        .after_help(
            r###"
This command wraps lastz to perform alignments suitable for the Cactus RepeatMasking workflow.

It handles:
*   Parallel execution for multiple target files.
*   Adding required modifiers: [multiple][nameparse=darkspace].
*   Setting output format to LAV.
*   Setting query depth threshold: --querydepth=keep,nowarn:(period+3).

Presets from UCSC:
    set01: Hg17vsPanTro1 (Human vs Chimp)
           C=0 E=30 K=3000 L=2200 O=400 Y=3400 Q=similar
    set02: Hg19vsPanTro2 (Human vs Primate, more sensitive)
           C=0 E=150 H=2000 K=4500 L=2200 M=254 O=600 T=2 Y=15000 Q=similar2
    set03: Hg17vsMm5 (Human vs Mouse)
           C=0 E=30 K=3000 L=2200 O=400 Q=default
    set04: Hg17vsRheMac2 (Human vs Macaque)
           C=0 E=30 H=2000 K=3000 L=2200 O=400 Q=default
    set05: Hg17vsBosTau2 (Human vs Cow)
           C=0 E=30 H=2000 K=3000 L=2200 M=50 O=400 Q=default
    set06: Hg17vsDanRer3 (Human vs Zebrafish)
           C=0 E=30 H=2000 K=2200 L=6000 O=400 Y=3400 Q=distant
    set07: Hg17vsMonDom1 (Human vs Opossum)
           C=0 E=30 H=2000 K=2200 L=10000 O=400 Y=3400 Q=distant

Examples:
    # Single target with set01
    pgr lav lastz query.fa target.fa --preset set01 -o lastz_out

    # Multiple targets with set03
    pgr lav lastz query.fa target_chr1.fa target_chr2.fa --preset set03 -o lastz_out

    # Show parameters and matrix for set01
    pgr lav lastz --preset set01 --show-preset

"###,
        )
        .arg(
            Arg::new("query")
                .required(true)
                .index(1)
                .help("Query FASTA file (fragments)"),
        )
        .arg(
            Arg::new("target")
                .required(true)
                .index(2)
                .num_args(1..)
                .help("Target FASTA file(s)"),
        )
        .arg(
            Arg::new("period")
                .long("period")
                .default_value("10")
                .value_parser(value_parser!(usize))
                .help("Period for querydepth calculation (depth = period + 3)"),
        )
        .arg(
            Arg::new("preset")
                .long("preset")
                .value_parser([
                    "set01",
                    "set02",
                    "set03",
                    "set04",
                    "set05",
                    "set06",
                    "set07",
                ])
                .help("Use predefined parameter sets (set01..set07)"),
        )
        .arg(
            Arg::new("show_preset")
                .long("show-preset")
                .action(ArgAction::SetTrue)
                .help("Display the configuration (parameters & matrix) for the selected preset and exit"),
        )
        .arg(
            Arg::new("lastz_args")
                .long("lastz-args")
                .help("Additional arguments passed directly to lastz (overrides preset)"),
        )
        .arg(
            Arg::new("output")
                .long("output")
                .short('o')
                .default_value("lastz_out")
                .help("Output directory"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .value_parser(value_parser!(usize))
                .default_value("4")
                .help("Number of parallel threads"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let preset = matches.get_one::<String>("preset");
    
    // Check if show-preset is requested
    if matches.get_flag("show_preset") {
        if let Some(preset_name) = preset {
            match preset_name.as_str() {
                "set01" => {
                    println!("Preset: set01");
                    println!("Description: Hg17vsPanTro1 (Human vs Chimp)");
                    println!("Parameters: C=0 E=30 K=3000 L=2200 O=400 Y=3400 Q=similar");
                    println!("\nMatrix Content:\n{}", MATRIX_SIMILAR);
                }
                "set02" => {
                    println!("Preset: set02");
                    println!("Description: Hg19vsPanTro2 (Human vs Primate, more sensitive)");
                    println!("Parameters: C=0 E=150 H=2000 K=4500 L=2200 M=254 O=600 T=2 Y=15000 Q=similar2");
                    println!("\nMatrix Content:\n{}", MATRIX_SIMILAR2);
                }
                "set03" => {
                    println!("Preset: set03");
                    println!("Description: Hg17vsMm5 (Human vs Mouse)");
                    println!("Parameters: C=0 E=30 K=3000 L=2200 O=400 Q=default");
                    println!("\nMatrix Content:\n{}", MATRIX_DEFAULT);
                }
                "set04" => {
                    println!("Preset: set04");
                    println!("Description: Hg17vsRheMac2 (Human vs Macaque)");
                    println!("Parameters: C=0 E=30 H=2000 K=3000 L=2200 O=400 Q=default");
                    println!("\nMatrix Content:\n{}", MATRIX_DEFAULT);
                }
                "set05" => {
                    println!("Preset: set05");
                    println!("Description: Hg17vsBosTau2 (Human vs Cow)");
                    println!("Parameters: C=0 E=30 H=2000 K=3000 L=2200 M=50 O=400 Q=default");
                    println!("\nMatrix Content:\n{}", MATRIX_DEFAULT);
                }
                "set06" => {
                    println!("Preset: set06");
                    println!("Description: Hg17vsDanRer3 (Human vs Zebrafish)");
                    println!("Parameters: C=0 E=30 H=2000 K=2200 L=6000 O=400 Y=3400 Q=distant");
                    println!("\nMatrix Content:\n{}", MATRIX_DISTANT);
                }
                "set07" => {
                    println!("Preset: set07");
                    println!("Description: Hg17vsMonDom1 (Human vs Opossum)");
                    println!("Parameters: C=0 E=30 H=2000 K=2200 L=10000 O=400 Y=3400 Q=distant");
                    println!("\nMatrix Content:\n{}", MATRIX_DISTANT);
                }
                _ => unreachable!(),
            }
            return Ok(());
        } else {
            anyhow::bail!("--show-preset requires --preset to be specified.");
        }
    }

    let query_file = matches.get_one::<String>("query").unwrap();
    let targets: Vec<_> = matches.get_many::<String>("target").unwrap().collect();
    let period = *matches.get_one::<usize>("period").unwrap();
    let lastz_args = matches.get_one::<String>("lastz_args");
    let output_dir = matches.get_one::<String>("output").unwrap();
    let parallel = *matches.get_one::<usize>("parallel").unwrap();

    // Check if lastz is installed
    if which::which("lastz").is_err() {
        anyhow::bail!("lastz not found in PATH. Please install lastz first.");
    }

    std::fs::create_dir_all(output_dir)?;

    // Handle targets
    // Unlike previous version, we process each target separately if they are separate files
    // But wait, lastz.pm logic was:
    // "Lastz will take the first sequence in target fasta files... For less confusions, each fasta files should contain only one sequence."
    // So we should iterate over target files.

    // Calculate depth
    let depth = period + 3;

    // Prepare matrix file if needed
    let mut _temp_matrix_handle: Option<NamedTempFile> = None;
    let mut matrix_arg = String::new();

    if let Some(preset_name) = preset {
        let mut t = NamedTempFile::new()?;
        match preset_name.as_str() {
            "set01" => {
                t.write_all(MATRIX_SIMILAR.as_bytes())?;
            }
            "set02" => {
                t.write_all(MATRIX_SIMILAR2.as_bytes())?;
            }
            "set03" | "set04" | "set05" => {
                t.write_all(MATRIX_DEFAULT.as_bytes())?;
            }
            "set06" | "set07" => {
                t.write_all(MATRIX_DISTANT.as_bytes())?;
            }
            _ => unreachable!(),
        }
        let path = t.path().to_string_lossy().to_string();
        matrix_arg = format!("Q={}", path);
        _temp_matrix_handle = Some(t);
    }

    // Build common args
    let mut common_args = Vec::new();
    common_args.push(format!("--querydepth=keep,nowarn:{}", depth));
    common_args.push("--format=lav".to_string());
    common_args.push("--markend".to_string());

    if let Some(preset_name) = preset {
         match preset_name.as_str() {
            "set01" => {
                common_args.push("C=0".to_string());
                common_args.push("E=30".to_string());
                common_args.push("K=3000".to_string());
                common_args.push("L=2200".to_string());
                common_args.push("O=400".to_string());
                common_args.push("Y=3400".to_string());
            }
            "set02" => {
                common_args.push("C=0".to_string());
                common_args.push("E=150".to_string());
                common_args.push("H=2000".to_string());
                common_args.push("K=4500".to_string());
                common_args.push("L=2200".to_string());
                common_args.push("M=254".to_string());
                common_args.push("O=600".to_string());
                common_args.push("T=2".to_string());
                common_args.push("Y=15000".to_string());
            }
            "set03" => {
                common_args.push("C=0".to_string());
                common_args.push("E=30".to_string());
                common_args.push("K=3000".to_string());
                common_args.push("L=2200".to_string());
                common_args.push("O=400".to_string());
            }
            "set04" => {
                common_args.push("C=0".to_string());
                common_args.push("E=30".to_string());
                common_args.push("H=2000".to_string());
                common_args.push("K=3000".to_string());
                common_args.push("L=2200".to_string());
                common_args.push("O=400".to_string());
            }
            "set05" => {
                common_args.push("C=0".to_string());
                common_args.push("E=30".to_string());
                common_args.push("H=2000".to_string());
                common_args.push("K=3000".to_string());
                common_args.push("L=2200".to_string());
                common_args.push("M=50".to_string());
                common_args.push("O=400".to_string());
            }
            "set06" => {
                common_args.push("C=0".to_string());
                common_args.push("E=30".to_string());
                common_args.push("H=2000".to_string());
                common_args.push("K=2200".to_string());
                common_args.push("L=6000".to_string());
                common_args.push("O=400".to_string());
                common_args.push("Y=3400".to_string());
            }
            "set07" => {
                common_args.push("C=0".to_string());
                common_args.push("E=30".to_string());
                common_args.push("H=2000".to_string());
                common_args.push("K=2200".to_string());
                common_args.push("L=10000".to_string());
                common_args.push("O=400".to_string());
                common_args.push("Y=3400".to_string());
            }
            _ => unreachable!(),
        }
        common_args.push(matrix_arg.clone());
    }

    if let Some(args) = lastz_args {
        for arg in args.split_whitespace() {
            common_args.push(arg.to_string());
        }
    }

    // Parallel execution
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(parallel)
        .build()
        .unwrap();

    pool.install(|| {
        targets.par_iter().for_each(|target_file| {
             let t_base = std::path::Path::new(target_file)
                .file_stem()
                .unwrap()
                .to_string_lossy();
            let q_base = std::path::Path::new(query_file)
                .file_stem()
                .unwrap()
                .to_string_lossy();

            // Output filename: [target]vs[query].lav
            // Note: If multiple targets have same basename, this will overwrite. 
            // Assuming unique basenames for now as per lastz.pm logic.
            let out_name = format!("[{}]vs[{}].lav", t_base, q_base);
            let out_path = std::path::Path::new(output_dir).join(out_name);

            // Note: lastz expects "file[mod]" as a single argument
            // lastz.pm didn't strictly use [multiple] for single files, but it's safer.
            // However, if we follow lastz.pm exactly:
            // "Lastz will take the first sequence in target fasta files and all sequences in query fasta files."
            // So for target we might want to be careful.
            // Let's stick to simple "file" for now unless user explicitly wants modifiers.
            // Wait, Cactus mode REQUIRED modifiers. 
            // But now we are in "Hybrid" mode.
            // Let's apply [multiple][nameparse=darkspace] to be safe and consistent with Cactus/general best practices.
            
            let target_arg = format!("{}[multiple][nameparse=darkspace]", target_file);
            let query_arg = format!("{}[nameparse=darkspace]", query_file);

            let mut cmd = std::process::Command::new("lastz");
            cmd.arg(&target_arg).arg(&query_arg);
            
            for arg in &common_args {
                cmd.arg(arg);
            }

            // Redirect stdout to file
            let file = std::fs::File::create(&out_path).expect("Failed to create output file");
            cmd.stdout(file);

            let status = cmd.status().expect("Failed to execute lastz");
            
            if !status.success() {
                eprintln!("Warning: lastz failed for {} vs {}", t_base, q_base);
            } else {
                println!("Finished: {} vs {}", t_base, q_base);
            }
        });
    });

    Ok(())
}
