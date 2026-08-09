use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for spectra.
pub fn make_subcommand() -> Command {
    Command::new("spectra")
        .about("Plots k-mer spectra (GenomeScope-style) as LaTeX")
        .after_help(
            r###"
Plots the k-mer coverage spectra in the four GenomeScope views (linear,
transformed, log, transformed-log) into one standalone LaTeX document:
observed histogram, fitted full model, error region, and k-mer peak
markers, with the model summary in the title.

The histogram is a `.hist` file (`pgr kmer hist`); the model is a
`model.txt` written by `pgr kmer gsize --model`.

* Compile with tectonic to get a PDF

Examples:
1. Plot from a histogram and model:
   pgr kmer hist reads.fq.gz -k 21 -o reads.hist
   pgr kmer gsize reads.fq.gz -k 21 --model -o gs_out
   pgr plot spectra reads.hist gs_out/model.txt -o spectra.tex
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input .hist histogram file",
        ))
        .arg(
            Arg::new("model")
                .required(true)
                .num_args(1)
                .help("Model.txt written by pgr kmer gsize --model"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the spectra command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let hist_file = args.get_one::<String>("infile").unwrap();
    let model_file = args.get_one::<String>("model").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let hist = pgr::libs::kmer::hist::load(std::path::Path::new(hist_file))
        .with_context(|| format!("failed to load histogram {hist_file}"))?;
    let (params, se_kmercov, p) = parse_model(model_file)?;
    let k = hist.k;
    let pairs: Vec<(f64, f64)> = hist
        .hist
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(i, &c)| ((i + 1) as f64, c as f64))
        .collect();

    let err = pgr::libs::plot::spectra::compute_error_rate(&pairs, k, &params, se_kmercov);
    let summary = pgr::libs::plot::spectra::SpectraSummary {
        k,
        p,
        params,
        se_kmercov,
        len: params.length,
        unique: 1.0 - params.d,
        error_rate: err,
    };

    let mut w = pgr::writer(outfile)?;
    pgr::libs::plot::spectra::render_spectra(&mut w, &pairs, &summary)?;
    w.flush()?;
    log::info!(
        "==> Wrote k-mer spectra (k={k}, kmercov={:.2}) to {}",
        params.kmercov,
        outfile
    );
    Ok(())
}

/// Parse `d/kmercov/bias/length[/r1]` rows with Estimate and Std. Error.
fn parse_model(
    path: &str,
) -> anyhow::Result<(pgr::libs::kmer::genomescope::ModelParams, f64, usize)> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read model file {path}"))?;
    let mut d = 0.0f64;
    let mut kmercov = 0.0f64;
    let mut bias = 0.0f64;
    let mut length = 0.0f64;
    let mut r1 = 0.0f64;
    let mut se_kmercov = 0.0f64;
    let mut has_r1 = false;
    for line in text.lines() {
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        let est = cols[1].parse::<f64>().unwrap_or(f64::NAN);
        let se = cols
            .get(2)
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
        match cols[0] {
            "d" => d = est,
            "kmercov" => {
                kmercov = est;
                se_kmercov = se;
            }
            "bias" => bias = est,
            "length" => length = est,
            "r1" => {
                r1 = est;
                has_r1 = true;
            }
            _ => {}
        }
    }
    anyhow::ensure!(kmercov > 0.0, "model file must contain a kmercov estimate");
    let params = pgr::libs::kmer::genomescope::ModelParams {
        d,
        kmercov,
        bias,
        length,
        r1,
    };
    let p = if has_r1 { 2 } else { 1 };
    Ok((params, se_kmercov, p))
}
