//! Lastz aligner presets and scoring matrices ported from UCSC.

use rayon::prelude::*;
use std::path::PathBuf;

/// Default scoring matrix for lastz (Human vs Mouse / Macaque / Cow).
pub const MATRIX_DEFAULT: &str = "   A    C    G    T
A  91 -114  -31 -123
C -114  100 -125  -31
G  -31 -125  100 -114
T -123  -31 -114   91
";

/// Distant-species scoring matrix (Human vs Zebrafish / Opossum).
pub const MATRIX_DISTANT: &str = "   A    C    G    T
A  91  -90  -25 -100
C -90  100 -100  -25
G -25 -100  100  -90
T -100  -25  -90   91
";

/// Close-species scoring matrix (Human vs Chimp).
pub const MATRIX_SIMILAR: &str = "   A    C    G    T
A  100 -300 -150 -300
C -300  100 -300 -150
G -150 -300  100 -300
T -300 -150 -300  100
";

/// Close-species scoring matrix variant (Human vs Primate, more sensitive).
pub const MATRIX_SIMILAR2: &str = "   A    C    G    T
A  90 -330 -236 -356
C -330  100 -318 -236
G -236 -318  100 -330
T -356 -236 -330   90
";

/// A predefined lastz parameter set with optional scoring matrix.
#[derive(Debug)]
pub struct Preset {
    pub name: &'static str,
    pub desc: &'static str,
    pub params: &'static str,
    pub matrix: Option<&'static str>,
}

/// UCSC-derived lastz presets for common pairwise vertebrate alignments.
pub const PRESETS: &[Preset] = &[
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

/// Look up a preset by name.
pub fn find_preset(name: &str) -> Option<&'static Preset> {
    PRESETS.iter().find(|p| p.name == name)
}

/// Collect all preset names (for clap PossibleValuesParser).
pub fn preset_names() -> Vec<&'static str> {
    PRESETS.iter().map(|p| p.name).collect()
}

/// Build the preset help string used in `--help` output.
pub fn preset_help() -> String {
    let mut help = String::from("Presets from UCSC:\n");
    for p in PRESETS {
        help.push_str(&format!(
            "    {}: {}\n           {}\n",
            p.name, p.desc, p.params
        ));
    }
    help
}

/// Options controlling a batch lastz run.
pub struct RunLastzOptions {
    /// Query depth threshold (informational; the actual flag lives in `common_args`).
    pub depth: usize,
    /// Self-alignment mode: skip pairs with different basenames and use `--self`
    /// when target and query paths are identical.
    pub is_self: bool,
    /// Arguments passed to every lastz invocation (depth, format, preset, user args).
    pub common_args: Vec<String>,
    /// Output directory (created if missing).
    pub output_dir: String,
    /// Number of worker threads.
    pub parallel: usize,
}

/// Run lastz for the cartesian product of `target_files` and `query_files`.
///
/// For each (target, query) pair, lastz is invoked with `opts.common_args` and
/// the result is written to `[t_base]vs[q_base].lav` in `opts.output_dir`.
/// Filename collisions are resolved by appending an incrementing counter
/// (`[t]vs[q].1.lav`, ...). When `opts.is_self` is set, pairs with different
/// file basenames are skipped, and identical target/query paths use lastz's
/// `--self` flag instead of passing a separate query argument.
pub fn run_lastz(
    target_files: Vec<PathBuf>,
    query_files: Vec<PathBuf>,
    opts: RunLastzOptions,
) -> anyhow::Result<()> {
    std::fs::create_dir_all(&opts.output_dir)?;

    let n_targets = target_files.len();
    let n_queries = query_files.len();
    let mut jobs: Vec<(PathBuf, PathBuf)> = Vec::with_capacity(n_targets * n_queries);
    for t in &target_files {
        for q in &query_files {
            jobs.push((t.clone(), q.clone()));
        }
    }

    log::info!("* Target files: [{}]", n_targets);
    log::info!("* Query files:  [{}]", n_queries);
    log::info!("* Total jobs:   [{}]", jobs.len());

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.parallel)
        .build()?;

    let common_args = opts.common_args;
    let output_dir = opts.output_dir;
    let is_self = opts.is_self;

    pool.install(move || {
        jobs.par_iter().for_each(|(target_file, query_file)| {
            let t_base =
                crate::libs::io::get_basename(&target_file.to_string_lossy()).unwrap_or_default();
            let q_base =
                crate::libs::io::get_basename(&query_file.to_string_lossy()).unwrap_or_default();

            if is_self && t_base != q_base {
                return;
            }

            // Output filename: [target]vs[query].lav.
            // Logic ported from lastz.pm to handle potential duplicates:
            // atomically reserve the name via create_new to prevent race
            // conditions when multiple threads process identically named inputs.
            let mut i = 0;
            let out_path;
            loop {
                let out_name = if i == 0 {
                    format!("[{}]vs[{}].lav", t_base, q_base)
                } else {
                    format!("[{}]vs[{}].{}.lav", t_base, q_base, i)
                };
                let candidate = std::path::Path::new(&output_dir).join(out_name);

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
            let target_arg = format!("{}[nameparse=darkspace]", target_file.to_string_lossy());

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

            log::info!("{:?}", cmd);

            match cmd.status() {
                Ok(status) if status.success() => {}
                Ok(status) => {
                    log::warn!(
                        "lastz failed (exit {:?}) for {} vs {}",
                        status,
                        t_base,
                        q_base
                    );
                }
                Err(err) => {
                    log::error!(
                        "failed to spawn lastz for {} vs {}: {}",
                        t_base,
                        q_base,
                        err
                    );
                }
            }
        });
    });

    Ok(())
}
