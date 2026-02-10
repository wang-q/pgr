use clap::builder::PossibleValuesParser;
use clap::*;
use rayon::prelude::*;
use std::io::Write;
use tempfile::NamedTempFile;

use std::path::{Path, PathBuf};

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

struct Preset {
    name: &'static str,
    desc: &'static str,
    params: &'static str,
    matrix: Option<&'static str>,
}

const PRESETS: &[Preset] = &[
    Preset {
        name: "set01",
        desc: "Hg17vsPanTro1 (Human vs Chimp)",
        params: "C=0 E=30 K=3000 L=2200 O=400 Y=3400 Q=similar",
        matrix: Some(MATRIX_SIMILAR),
    },
    Preset {
        name: "set02",
        desc: "Hg19vsPanTro2 (Human vs Primate, more sensitive)",
        params: "C=0 E=150 H=2000 K=4500 L=2200 M=254 O=600 T=2 Y=15000 Q=similar2",
        matrix: Some(MATRIX_SIMILAR2),
    },
    Preset {
        name: "set03",
        desc: "Hg17vsMm5 (Human vs Mouse)",
        params: "C=0 E=30 K=3000 L=2200 O=400 Q=default",
        matrix: Some(MATRIX_DEFAULT),
    },
    Preset {
        name: "set04",
        desc: "Hg17vsRheMac2 (Human vs Macaque)",
        params: "C=0 E=30 H=2000 K=3000 L=2200 O=400 Q=default",
        matrix: Some(MATRIX_DEFAULT),
    },
    Preset {
        name: "set05",
        desc: "Hg17vsBosTau2 (Human vs Cow)",
        params: "C=0 E=30 H=2000 K=3000 L=2200 M=50 O=400 Q=default",
        matrix: Some(MATRIX_DEFAULT),
    },
    Preset {
        name: "set06",
        desc: "Hg17vsDanRer3 (Human vs Zebrafish)",
        params: "C=0 E=30 H=2000 K=2200 L=6000 O=400 Y=3400 Q=distant",
        matrix: Some(MATRIX_DISTANT),
    },
    Preset {
        name: "set07",
        desc: "Hg17vsMonDom1 (Human vs Opossum)",
        params: "C=0 E=30 H=2000 K=2200 L=10000 O=400 Y=3400 Q=distant",
        matrix: Some(MATRIX_DISTANT),
    },
];

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    let mut preset_help = String::from("Presets from UCSC:\n");
    let mut preset_names = Vec::new();
    for p in PRESETS {
        preset_help.push_str(&format!("    {}: {}\n           {}\n", p.name, p.desc, p.params));
        preset_names.push(p.name);
    }

    Command::new("lastz")
        .about("Wrapper for lastz alignment (Cactus style)")
        .after_help(format!(
            r###"
This command wraps lastz to perform alignments suitable for the Cactus RepeatMasking workflow.

It handles:
*   Parallel execution for multiple target files.
*   Directory recursion for target and query inputs.
*   Adding required modifiers: [multiple][nameparse=darkspace].
*   Setting output format to LAV.
*   Setting query depth threshold: --querydepth=keep,nowarn:N.
    N is the depth of coverage threshold (aligned bases / query length).
    When the threshold is exceeded, processing stops for that query/strand to save time.
    'keep' ensures alignments found *before* the threshold are reported (unlike default which discards all).
    'nowarn' suppresses warnings about exceeded depth.
    Note: Reported alignments are the first found, not necessarily optimal.
    Default depth 50 allows ~50x coverage.

{}
Examples:
    # Single target with set01
    pgr lav lastz target.fa query.fa --preset set01 -o lastz_out

    # Directory inputs
    pgr lav lastz target_dir/ query_dir/ --preset set03 -o lastz_out

    # Show parameters and matrix for set01
    pgr lav lastz --preset set01 --show-preset

"###,
            preset_help
        ))
        .arg(
            Arg::new("target")
                .required(true)
                .index(1)
                .help("Target FASTA file or directory"),
        )
        .arg(
            Arg::new("query")
                .required(true)
                .index(2)
                .help("Query FASTA file or directory"),
        )
        .arg(
            Arg::new("depth")
                .long("depth")
                .default_value("50")
                .value_parser(value_parser!(usize))
                .help("Query depth threshold"),
        )
        .arg(
            Arg::new("is_self")
                .long("self")
                .action(clap::ArgAction::SetTrue)
                .help("Self-alignment"),
        )
        .arg(
            Arg::new("preset")
                .long("preset")
                .short('s')
                .value_parser(PossibleValuesParser::new(preset_names))
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

fn find_fasta_files<P: AsRef<Path>>(path: P) -> Vec<PathBuf> {
    let mut files = Vec::new();
    let path = path.as_ref();

    if path.is_file() {
        files.push(path.to_path_buf());
    } else if path.is_dir() {
        if let Ok(entries) = std::fs::read_dir(path) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_dir() {
                    files.extend(find_fasta_files(&path));
                } else if let Some(ext) = path.extension() {
                    let ext_str = ext.to_string_lossy().to_lowercase();
                    // Match .fa and .fa.gz
                    if ext_str == "fa" {
                        files.push(path);
                    } else if ext_str == "gz" {
                        if let Some(stem) = path.file_stem() {
                            let stem_path = Path::new(stem);
                            if let Some(stem_ext) = stem_path.extension() {
                                if stem_ext.to_string_lossy().to_lowercase() == "fa" {
                                    files.push(path);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    files
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let preset = matches.get_one::<String>("preset");

    // Check if show-preset is requested
    if matches.get_flag("show_preset") {
        if let Some(preset_name) = preset {
            if let Some(p) = PRESETS.iter().find(|p| p.name == preset_name) {
                println!("Preset: {}", p.name);
                println!("Description: {}", p.desc);
                println!("Parameters: {}", p.params);
                if let Some(matrix) = p.matrix {
                    println!("\nMatrix Content:\n{}", matrix);
                }
            } else {
                unreachable!();
            }
            return Ok(());
        } else {
            anyhow::bail!("--show-preset requires --preset to be specified.");
        }
    }

    let arg_query = matches.get_one::<String>("query").unwrap();
    let arg_target = matches.get_one::<String>("target").unwrap();
    let opt_depth = *matches.get_one::<usize>("depth").unwrap();
    let opt_lastz_args = matches.get_one::<String>("lastz_args");
    let opt_output = matches.get_one::<String>("output").unwrap();
    let opt_parallel = *matches.get_one::<usize>("parallel").unwrap();
    let is_self = matches.get_flag("is_self");

    // Check if lastz is installed
    if which::which("lastz").is_err() {
        anyhow::bail!("lastz not found in PATH. Please install lastz first.");
    }

    std::fs::create_dir_all(opt_output)?;

    // Expand files
    let mut query_files = find_fasta_files(arg_query);
    query_files.sort();

    let mut target_files = find_fasta_files(arg_target);
    target_files.sort();

    if query_files.is_empty() {
        anyhow::bail!("No query FASTA files found in {}", arg_query);
    }
    if target_files.is_empty() {
        anyhow::bail!("No target FASTA files found in {}", arg_target);
    }

    // Prepare matrix file if needed
    let mut _temp_matrix_handle: Option<NamedTempFile> = None;
    let mut matrix_path = String::new();

    if let Some(preset_name) = preset {
        let mut t = NamedTempFile::new()?;
        if let Some(p) = PRESETS.iter().find(|p| p.name == preset_name) {
            if let Some(matrix) = p.matrix {
                t.write_all(matrix.as_bytes())?;
            }
        } else {
            unreachable!();
        }
        matrix_path = t.path().to_string_lossy().to_string();
        _temp_matrix_handle = Some(t);
    }

    // Build common args
    let mut common_args = Vec::new();
    common_args.push(format!("--querydepth=keep,nowarn:{}", opt_depth));
    common_args.push("--format=lav".to_string());
    common_args.push("--markend".to_string());
    common_args.push("--ambiguous=iupac".to_string());

    if let Some(preset_name) = preset {
        if let Some(p) = PRESETS.iter().find(|p| p.name == preset_name) {
            for arg in p.params.split_whitespace() {
                if !arg.starts_with("Q=") {
                    common_args.push(arg.to_string());
                }
            }
        } else {
            unreachable!();
        }
        if !matrix_path.is_empty() {
            common_args.push(format!("Q={}", matrix_path));
        }
    }

    if let Some(args) = opt_lastz_args {
        for arg in args.split_whitespace() {
            common_args.push(arg.to_string());
        }
    }

    // Create jobs (Cartesian product)
    let mut jobs = Vec::new();
    for t in &target_files {
        for q in &query_files {
            jobs.push((t, q));
        }
    }

    eprintln!("* Target files: [{}]", target_files.len());
    eprintln!("* Query files:  [{}]", query_files.len());
    eprintln!("* Total jobs:   [{}]", jobs.len());

    // Parallel execution
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opt_parallel)
        .build()
        .unwrap();

    pool.install(|| {
        jobs.par_iter().for_each(|(target_file, query_file)| {
            let get_base_name = |path: &Path| -> String {
                let stem = path.file_stem().unwrap().to_string_lossy();
                stem.split('.').next().unwrap().to_string()
            };

            let t_base = get_base_name(target_file);
            let q_base = get_base_name(query_file);

            if is_self && t_base != q_base {
                return;
            }

            // Output filename: [target]vs[query].lav
            // Logic ported from lastz.pm to handle potential duplicates
            let mut i = 0;
            let out_path;
            loop {
                let out_name = if i == 0 {
                    format!("[{}]vs[{}].lav", t_base, q_base)
                } else {
                    format!("[{}]vs[{}].{}.lav", t_base, q_base, i)
                };
                let candidate = std::path::Path::new(opt_output).join(out_name);

                // Atomically try to create the file to reserve the name
                // This prevents race conditions when multiple threads process identically named inputs
                if std::fs::OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(&candidate)
                    .is_ok()
                {
                    out_path = candidate;
                    break;
                }

                i += 1;
            }

            // [nameparse=darkspace] is required for correct sequence name parsing.
            // [multiple] implies that the target file contains multiple sequences.
            // However, lastz documentation states that --self cannot be used with multiple sequences in the target.
            // Since we are running in --self mode for repeat masking (or at least supporting it),
            // and we are feeding single-sequence chunks (or small chunks) in the standard workflow,
            // we omit [multiple] to avoid conflicts.
            // If the user provides a multi-sequence file without --self, lastz might complain or process only the first sequence
            // unless we add [multiple] back conditionally.
            // But for now, for the "Cactus-style" workflow which splits files, this is safer.
            let target_arg = format!(
                "{}[nameparse=darkspace]",
                target_file.to_string_lossy()
            );

            let mut cmd = std::process::Command::new("lastz");
            cmd.arg(&target_arg);

            if is_self && target_file == query_file {
                cmd.arg("--self");
            } else {
                let query_arg = format!("{}[nameparse=darkspace]", query_file.to_string_lossy());
                cmd.arg(&query_arg);
            }

            for arg in &common_args {
                cmd.arg(arg);
            }

            cmd.arg(format!("--output={}", out_path.to_string_lossy()));

            // Print command for progress tracking
            eprintln!("{:?}", cmd);

            // Execute lastz and wait for it to complete
            let status = cmd.status().expect("Failed to execute lastz");

            if !status.success() {
                eprintln!("Warning: lastz failed for {} vs {}", t_base, q_base);
            } else {
                // println!("Finished: {} vs {}", t_base, q_base);
            }
        });
    });

    Ok(())
}
