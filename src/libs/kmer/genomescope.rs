//! GenomeScope-style negative-binomial mixture fitting (genescopefk.R port).

use crate::libs::lm::{fdjac2, invert_matrix, lmdif, qrfac};

/// Number of optimization rounds (R `NUM_ROUNDS`).
pub const NUM_ROUNDS: usize = 4;
/// Coverage steps trimmed off between rounds (R `START_SHIFT`).
pub const START_SHIFT: usize = 5;
/// Typical sequencing-error cutoff (R `TYPICAL_ERROR`).
pub const TYPICAL_ERROR: usize = 15;
/// Max iterations for the LM fit (R `MAX_ITERATIONS`).
pub const MAX_ITERATIONS: usize = 200;

/// Natural log of the gamma function (Lanczos approximation, GSL-style).
#[allow(clippy::excessive_precision)]
fn lgamma(x: f64) -> f64 {
    const C: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI.ln() - (std::f64::consts::PI * x).sin().ln() - lgamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let t = x + 7.5;
        let mut a = C[0];
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Negative binomial PMF, R `dnbinom(x, size, mu)` semantics.
fn dnbinom(x: f64, size: f64, mu: f64) -> f64 {
    if mu <= 0.0 {
        return 0.0;
    }
    if !size.is_finite() {
        // Poisson limit as size -> Inf (R dnbinom handles size = Inf).
        return (-mu + x * mu.ln() - lgamma(x + 1.0)).exp();
    }
    if size <= 0.0 {
        return 0.0;
    }
    // log P = lgamma(x+size) - lgamma(size) - lgamma(x+1)
    //        + size*ln(size/(size+mu)) + x*ln(mu/(size+mu))
    let log_p = lgamma(x + size) - lgamma(size) - lgamma(x + 1.0)
        + size * (size / (size + mu)).ln()
        + x * (mu / (size + mu)).ln();
    log_p.exp()
}

/// Model parameters shared by the fitted mixtures.
#[derive(Debug, Clone, Copy)]
pub struct ModelParams {
    /// Repetitiveness (fraction of duplicated k-mers).
    pub d: f64,
    /// K-mer coverage of the homozygous peak.
    pub kmercov: f64,
    /// Negative-binomial overdispersion (R `bias`).
    pub bias: f64,
    /// Genome k-mer positions (R `length`).
    pub length: f64,
    /// Nucleotide heterozygosity (p=2, `r1`).
    pub r1: f64,
}

/// Fitted model with point estimates and 2-SE ranges (R `nls_peak` output).
#[derive(Debug, Clone)]
pub struct FittedModel {
    /// Ploidy.
    pub p: usize,
    /// Topology.
    pub top: usize,
    /// Point estimates.
    pub params: ModelParams,
    /// Deviance (residual sum of squares of the transformed fit).
    pub deviance: f64,
    /// Standard errors of the fitted parameters (R summary() order).
    pub se: Vec<f64>,
    /// 2-SE ranges for d, kmercov, bias, length.
    pub d_range: [f64; 2],
    pub kcov_range: [f64; 2],
    pub bias_range: [f64; 2],
    pub length_range: [f64; 2],
    /// Heterozygosity range (p=2 only; p=1 is [0, 0]).
    pub het_range: [f64; 2],
}

/// Model mixture prediction `predict_p_top(...)`: the negative-binomial
/// mixture without the `x * length` scale factor (the formula layer applies
/// `x^transform_exp * length * predict(...)`, matching R).
fn predict(p: usize, top: usize, params: &ModelParams, k: usize, x: f64) -> f64 {
    match (p, top) {
        (1, 1) | (1, 0) => {
            // r0 = 1, t0 = s0 = 1
            let alpha1 = 1.0 - params.d;
            let alpha2 = params.d;
            alpha1 * dnbinom(x, params.kmercov / params.bias, params.kmercov)
                + alpha2 * dnbinom(x, 2.0 * params.kmercov / params.bias, 2.0 * params.kmercov)
        }
        (2, 1) | (2, 0) => {
            let r0 = 1.0 - params.r1;
            if r0 < 0.0 || params.d > 1.0 {
                return 0.0;
            }
            let t0 = r0.powi(k as i32);
            let s0 = t0;
            let s1 = 1.0 - t0;
            let alpha1 = (1.0 - params.d) * (2.0 * s1) + params.d * (2.0 * s0 * s1 + 2.0 * s1 * s1);
            let alpha2 = (1.0 - params.d) * s0 + params.d * (s1 * s1);
            let alpha3 = params.d * (2.0 * s0 * s1);
            let alpha4 = params.d * (s0 * s0);
            alpha1 * dnbinom(x, params.kmercov / params.bias, params.kmercov)
                + alpha2 * dnbinom(x, 2.0 * params.kmercov / params.bias, 2.0 * params.kmercov)
                + alpha3 * dnbinom(x, 3.0 * params.kmercov / params.bias, 3.0 * params.kmercov)
                + alpha4 * dnbinom(x, 4.0 * params.kmercov / params.bias, 4.0 * params.kmercov)
        }
        _ => 0.0,
    }
}
pub struct GenomeScopeResult {
    /// Ploidy.
    pub p: usize,
    /// K-mer length.
    pub k: usize,
    /// Fitted k-mer coverage.
    pub kmercov: f64,
    /// Negative-binomial overdispersion.
    pub bias: f64,
    /// Standard errors of the fitted parameters (R parameter order).
    pub se: Vec<f64>,
    /// Repetitiveness.
    pub d: f64,
    /// Genome k-mer positions.
    pub length: f64,
    /// Nucleotide heterozygosity (p=2).
    pub het: f64,
    /// Genome haploid length (min, max over 2 SE).
    pub genome_haploid: [f64; 2],
    /// Repeat length range.
    pub repeat_len: [f64; 2],
    /// Unique length range.
    pub unique_len: [f64; 2],
    /// Model fit (percent kmers modeled, all/full).
    pub model_fit: [f64; 2],
    /// Read error rate range.
    pub error_rate: [f64; 2],
    /// Whether any model converged.
    pub converged: bool,
}

/// Fit the GenomeScope model to a k-mer histogram (counts indexed by
/// coverage, `hist[c-1]` = number of k-mers with count `c`).
///
/// Mirrors `genescopefk.R`: trim sequencing errors across `NUM_ROUNDS`
/// rounds, try `p` kmercov candidates (and 6 heterozygosity starts for
/// p=2), keep the lowest-deviance model, then score it.
pub fn fit(hist: &[u64], k: usize, p: usize) -> GenomeScopeResult {
    // Sparse (coverage, count) pairs, matching R's histogram file where
    // zero-count rows are absent; drop a leading zero and the trailing
    // position like R (`minkmerx` / `get rid of the last position`).
    let pairs: Vec<(f64, f64)> = hist
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(i, &c)| ((i + 1) as f64, c as f64))
        .collect();
    let minkmerx = if hist.first() == Some(&0) { 1 } else { 0 };
    let prof = &pairs[minkmerx..pairs.len().saturating_sub(1)];
    if prof.len() < 3 {
        return empty_result(p, k);
    }

    // Transformed profile and the error/peak anchors (transform_exp = 1).
    let kmer_trans: Vec<f64> = prof.iter().map(|&(x, c)| x * c).collect();
    let typical = TYPICAL_ERROR.min(kmer_trans.len());
    let min_err = kmer_trans[..typical]
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let start0 = kmer_trans[..typical]
        .iter()
        .rposition(|&v| v == min_err)
        .unwrap_or(0);
    let peak_val = kmer_trans[start0..]
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let start_max = start0
        + kmer_trans[start0..]
            .iter()
            .position(|&v| v == peak_val)
            .unwrap_or(0);

    let mut best: Option<(FittedModel, f64)> = None;
    let mut start = start0;
    let mut round = 0;
    let end = prof.len();
    while round < NUM_ROUNDS && start + 5 < end {
        let x: Vec<f64> = prof[start..end].iter().map(|&(xv, _)| xv).collect();
        let y: Vec<f64> = prof[start..end].iter().map(|&(_, c)| c).collect();
        if x.len() < 5 {
            break;
        }
        let model = estimate_peak(&x, &y, k, p, prof[start_max].0);
        if let Some(m) = &model {
            let dev = m.deviance;
            let keep = match &best {
                None => true,
                Some((b, _)) => {
                    let pdiff = (dev - b.deviance).abs() / dev.max(b.deviance);
                    if pdiff < 0.20 {
                        // Similar score: prefer higher heterozygosity.
                        m.params.r1 > b.params.r1
                    } else {
                        dev < b.deviance
                    }
                }
            };
            if keep {
                best = Some((m.clone(), dev));
            }
        }
        start += START_SHIFT;
        round += 1;
    }

    match best {
        Some((model, _)) => {
            let score = score_model(hist, &model, k);
            let total = model.params.length;
            let genome_haploid = [model.length_range[0].max(0.0), model.length_range[1]];
            GenomeScopeResult {
                p,
                k,
                kmercov: model.params.kmercov,
                bias: model.params.bias,
                se: model.se.clone(),
                d: model.params.d,
                length: model.params.length,
                het: model.params.r1,
                genome_haploid,
                repeat_len: [model.d_range[0] * total, model.d_range[1] * total],
                unique_len: [
                    (1.0 - model.d_range[1]) * total,
                    (1.0 - model.d_range[0]) * total,
                ],
                model_fit: [score.all_score, score.full_score],
                error_rate: [0.0, 0.0],
                converged: true,
            }
        }
        None => {
            let mut r = empty_result(p, k);
            r.converged = false;
            r
        }
    }
}

fn empty_result(p: usize, k: usize) -> GenomeScopeResult {
    GenomeScopeResult {
        p,
        k,
        kmercov: 0.0,
        bias: 0.0,
        se: vec![],
        d: 0.0,
        length: 0.0,
        het: 0.0,
        genome_haploid: [0.0, 0.0],
        repeat_len: [0.0, 0.0],
        unique_len: [0.0, 0.0],
        model_fit: [0.0, 0.0],
        error_rate: [0.0, 0.0],
        converged: false,
    }
}

/// Try the `p` kmercov candidates (R `estimate_Genome_peakp`) and return the
/// best model.
fn estimate_peak(
    x: &[f64],
    y: &[f64],
    k: usize,
    p: usize,
    est_kmercov: f64,
) -> Option<FittedModel> {
    let numof_kmers: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
    let mut best: Option<FittedModel> = None;
    for i in 1..=p {
        let kcov = est_kmercov / i as f64;
        let est_length = numof_kmers / kcov;
        if let Some(m) = nls_peak(x, y, k, p, kcov, est_length) {
            let keep = match &best {
                None => true,
                Some(b) => m.deviance < b.deviance,
            };
            if keep {
                best = Some(m);
            }
        }
    }
    best
}

/// Single LM fit with the R parameter initialization sets.
fn nls_peak(
    x: &[f64],
    y: &[f64],
    k: usize,
    p: usize,
    est_kmercov: f64,
    est_length: f64,
) -> Option<FittedModel> {
    // Truncate to the first 2000 points like R.
    let n = x.len().min(2000);
    let x = &x[..n];
    let y = &y[..n];
    let bias_init = 0.5;
    let mut best: Option<(ModelParams, f64, Vec<f64>, Vec<f64>)> = None;
    let mut best_dev = f64::INFINITY;
    // estLength = sum(x*count)/kmercov is already genome-scale; try starts
    // around it (the transformed profile carries an x factor, not x^3, so
    // no kmercov^2 correction is needed).
    let length_inits = [
        est_length / p as f64,
        est_length / p as f64 * 1.25,
        est_length / p as f64 * 0.8,
    ];
    let r_starts: Vec<f64> = if p == 1 {
        vec![0.0]
    } else {
        vec![0.001, 0.001 * 1.0, 0.001, 0.01, 0.01, 0.01]
    };
    // d starts span the plausible range; R uses a single 0.10 start, but
    // the extra starts help escape the d/length local optimum that the
    // numeric LM can otherwise settle into.
    let d_starts = [0.001, 0.10, 0.50];
    let maxfev = 100 * (spec_params_len(p) + 1);
    let residual = |params: &[f64]| -> Vec<f64> {
        let mp = params_to_model(params, p);
        x.iter()
            .zip(y)
            .map(|(&xi, &yi)| xi * yi - xi * mp.length * predict(p, 1, &mp, k, xi))
            .collect()
    };
    for length_init in length_inits {
        for d_init in d_starts {
            for r1 in &r_starts {
                let init: Vec<f64> = if p == 1 {
                    vec![d_init, est_kmercov, bias_init, length_init]
                } else {
                    vec![d_init, *r1, est_kmercov, bias_init, length_init]
                };
                let lower = vec![0.0; init.len()];
                let mut upper = vec![f64::INFINITY; init.len()];
                upper[0] = 1.0; // d <= 1
                if p == 2 {
                    upper[1] = 1.0; // r1 <= 1
                }
                let res = lmdif(
                    &residual, &init, &lower, &upper, 1.49e-8, // ftol = sqrt(eps)
                    1.49e-8, // ptol = sqrt(eps)
                    0.0,     // gtol (nls.lm default)
                    maxfev, 0.0, // epsfcn: machine precision
                    0.1, // factor (genescopefk control)
                );
                let fit_x = res.x;
                let fvec = res.fvec;
                let info = res.info;
                if (1..=8).contains(&info) {
                    let mp = params_to_model(&fit_x, p);
                    let dev: f64 = fvec.iter().map(|v| v * v).sum();
                    if dev < best_dev {
                        best_dev = dev;
                        best = Some((mp, dev, fit_x.clone(), fvec));
                    }
                }
            }
        }
    }
    let (params, deviance, fit_x, fvec) = best?;
    // Standard errors from the final Jacobian: hessian = P^T R^T R P,
    // SE = sqrt(diag(hessian^-1) * resvar), mirroring R's summary().
    let np = fit_x.len();
    let m = x.len();
    let fj = fdjac2(&residual, &fit_x, &fvec, 0.0, m);
    let mut fjac = fj;
    let mut ipvt = vec![0usize; np];
    let mut rdiag = vec![0.0; np];
    let mut acnorm = vec![0.0; np];
    qrfac(&mut fjac, m, np, &mut ipvt, &mut rdiag, &mut acnorm);
    // Build R (upper triangular, diagonal from rdiag).
    let mut r = vec![0.0; np * np];
    for j in 0..np {
        for i in 0..=j {
            r[j * np + i] = if i == j { rdiag[j] } else { fjac[j * m + i] };
        }
    }
    // R^T R.
    let mut rtr = vec![0.0; np * np];
    for a in 0..np {
        for b in 0..np {
            let mut s = 0.0;
            for i in 0..np {
                s += r[i * np + a] * r[i * np + b];
            }
            rtr[a * np + b] = s;
        }
    }
    // Permute to the original parameter order: hessian = P^T (R^T R) P.
    let mut hess = vec![0.0; np * np];
    for j1 in 0..np {
        for j2 in 0..np {
            let l1 = ipvt[j1];
            let l2 = ipvt[j2];
            hess[l1 * np + l2] = rtr[j1 * np + j2];
        }
    }
    let resvar = deviance / ((m - np).max(1) as f64);
    let se = match invert_matrix(&hess, np) {
        Some(inv) => (0..np)
            .map(|j| (inv[j * np + j].max(0.0) * resvar).sqrt())
            .collect(),
        None => vec![0.0; np],
    };
    Some(FittedModel {
        p,
        top: 1,
        params,
        deviance,
        se: se.clone(),
        d_range: [params.d - 2.0 * se[0], params.d + 2.0 * se[0]],
        kcov_range: [
            params.kmercov - 2.0 * se[if p == 1 { 1 } else { 2 }],
            params.kmercov + 2.0 * se[if p == 1 { 1 } else { 2 }],
        ],
        bias_range: [
            params.bias - 2.0 * se[if p == 1 { 2 } else { 3 }],
            params.bias + 2.0 * se[if p == 1 { 2 } else { 3 }],
        ],
        length_range: [
            params.length - 2.0 * se[if p == 1 { 3 } else { 4 }],
            params.length + 2.0 * se[if p == 1 { 3 } else { 4 }],
        ],
        het_range: [params.r1, params.r1],
    })
}

/// Number of fitted parameters for a ploidy.
fn spec_params_len(p: usize) -> usize {
    if p == 1 {
        4
    } else {
        5
    }
}

/// Map the fitted parameter vector to [`ModelParams`] (R parameter order:
/// d, r..., kmercov, bias, length).
fn params_to_model(params: &[f64], p: usize) -> ModelParams {
    if p == 1 {
        ModelParams {
            d: params[0],
            kmercov: params[1],
            bias: params[2],
            length: params[3],
            r1: 0.0,
        }
    } else {
        ModelParams {
            d: params[0],
            kmercov: params[2],
            bias: params[3],
            length: params[4],
            r1: params[1],
        }
    }
}

/// Model score (R `score_model`): residual sum of squares and percent of
/// kmers modeled, excluding the sequencing-error region.
struct ModelScore {
    all_score: f64,
    full_score: f64,
}

fn score_model(hist: &[u64], model: &FittedModel, k: usize) -> ModelScore {
    // Sparse (x, count) pairs, matching R's histogram file where zero-count
    // rows are absent.
    let pairs: Vec<(f64, f64)> = hist
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(i, &c)| ((i + 1) as f64, c as f64))
        .collect();
    let end = pairs.len();
    let y_transform: Vec<f64> = pairs.iter().map(|&(x, c)| x * c).collect();
    let pred: Vec<f64> = pairs
        .iter()
        .map(|&(x, _)| x * model.params.length * predict(model.p, model.top, &model.params, k, x))
        .collect();
    let se_k = model.se[if model.p == 1 { 1 } else { 2 }];
    let kcovfloor = ((model.params.kmercov - 2.0 * se_k).floor().max(1.0)) as usize;
    let err_cut = pairs
        .iter()
        .position(|&(x, _)| x > kcovfloor as f64)
        .unwrap_or(end);
    let recip = |i: usize| 1.0 / pairs[i].0;

    // Truncate error residuals once they drop below 1 (R `first_zero`).
    let mut first_zero = err_cut;
    for i in 0..err_cut {
        let e = recip(i) * (y_transform[i] - pred[i]);
        if e < 1.0 {
            first_zero = i + 1;
            break;
        }
    }
    let mut all_abs = 0.0;
    let mut all_sum = 0.0;
    for i in first_zero.saturating_sub(1)..end {
        let e = y_transform[i] - pred[i];
        all_abs += e.abs();
        all_sum += y_transform[i];
    }
    let full_end = ((2 * model.p + 1) * kcovfloor).min(end);
    let mut full_abs = 0.0;
    let mut full_sum = 0.0;
    for i in first_zero.saturating_sub(1)..full_end {
        let e = y_transform[i] - pred[i];
        full_abs += e.abs();
        full_sum += y_transform[i];
    }
    ModelScore {
        all_score: 1.0 - all_abs / all_sum.max(1e-300),
        full_score: 1.0 - full_abs / full_sum.max(1e-300),
    }
}

/// Write `summary.txt` and `model.txt` in the GenomeScope formats consumed
/// by anchr's `2_fastk` (`grep ^kmercov model.txt`).
pub fn write_outputs(outdir: &std::path::Path, result: &GenomeScopeResult) -> anyhow::Result<()> {
    std::fs::create_dir_all(outdir)?;
    let summary = outdir.join("summary.txt");
    let model = outdir.join("model.txt");
    let bp = |v: f64| format!("{:.0} bp", v);
    let pct = |v: f64| format!("{:.4}%", v * 100.0);
    if result.converged {
        let mut s = String::new();
        s.push_str("GenomeScope version 2.0 (pgr native)\n");
        s.push_str(&format!("p = {}\n", result.p));
        s.push_str(&format!("k = {}\n", result.k));
        s.push_str(&format!("{:<30}{:<18}{:<18}\n", "property", "min", "max"));
        if result.p == 1 {
            s.push_str(&format!(
                "{:<30}{:<18}{:<18}\n",
                "Homozygous (a)",
                pct(1.0),
                pct(1.0)
            ));
        } else {
            s.push_str(&format!(
                "{:<30}{:<18}{:<18}\n",
                "Homozygous (aa)",
                pct(1.0 - result.het),
                pct(1.0 - result.het)
            ));
            s.push_str(&format!(
                "{:<30}{:<18}{:<18}\n",
                "Heterozygous (ab)",
                pct(result.het),
                pct(result.het)
            ));
        }
        s.push_str(&format!(
            "{:<30}{:<18}{:<18}\n",
            "Genome Haploid Length",
            bp(result.genome_haploid[0]),
            bp(result.genome_haploid[1])
        ));
        s.push_str(&format!(
            "{:<30}{:<18}{:<18}\n",
            "Genome Repeat Length",
            bp(result.repeat_len[0]),
            bp(result.repeat_len[1])
        ));
        s.push_str(&format!(
            "{:<30}{:<18}{:<18}\n",
            "Genome Unique Length",
            bp(result.unique_len[0]),
            bp(result.unique_len[1])
        ));
        s.push_str(&format!(
            "{:<30}{:<18}{:<18}\n",
            "Model Fit",
            pct(result.model_fit[0]),
            pct(result.model_fit[1])
        ));
        std::fs::write(&summary, s)?;

        // model.txt mirrors an nls summary: parameters with an Estimate
        // column, so `grep '^kmercov' | cut -f 2` works.
        let mut m = String::new();
        m.push_str("Formula: y_transform ~ x*length*predict\n\nParameters:\n");
        m.push_str(&format!(
            "{:<12}{:<18}{:<18}\n",
            "", "Estimate", "Std. Error"
        ));
        let se = |i: usize| result.se.get(i).copied().unwrap_or(0.0);
        m.push_str(&format!("{:<12}{:<18.6}{:<18.6}\n", "d", result.d, se(0)));
        if result.p == 2 {
            m.push_str(&format!(
                "{:<12}{:<18.6}{:<18.6}\n",
                "r1",
                result.het,
                se(1)
            ));
        }
        m.push_str(&format!(
            "{:<12}{:<18.6}{:<18.6}\n",
            "kmercov",
            result.kmercov,
            se(if result.p == 1 { 1 } else { 2 })
        ));
        m.push_str(&format!(
            "{:<12}{:<18.6}{:<18.6}\n",
            "bias",
            result.bias,
            se(if result.p == 1 { 2 } else { 3 })
        ));
        m.push_str(&format!(
            "{:<12}{:<18.6}{:<18.6}\n",
            "length",
            result.length,
            se(if result.p == 1 { 3 } else { 4 })
        ));
        std::fs::write(&model, m)?;
    } else {
        std::fs::write(
            &summary,
            "GenomeScope version 2.0 (pgr native)\np = ...\nFailed to converge.\n",
        )?;
        std::fs::write(&model, "Failed to converge.")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lgamma_and_dnbinom_are_sane() {
        // lgamma(7) == ln(6!) == ln(720)
        assert!((lgamma(7.0) - (720.0f64).ln()).abs() < 1e-10);
        // dnbinom sums to ~1 over the support for a Poisson-like shape.
        let size = 20.0 / 0.5;
        let total: f64 = (0..200).map(|x| dnbinom(x as f64, size, 20.0)).sum();
        assert!((total - 1.0).abs() < 1e-6, "dnbinom sum {total}");
    }

    #[test]
    fn predicts_peak_at_coverage() {
        let mp = ModelParams {
            d: 0.0,
            kmercov: 30.0,
            bias: 0.5,
            length: 100_000.0,
            r1: 0.0,
        };
        let y = |x: f64| predict(1, 1, &mp, 21, x);
        let mut peak_x = 1.0;
        let mut peak_v = y(1.0);
        for i in 2..=100 {
            let v = y(i as f64);
            if v > peak_v {
                peak_v = v;
                peak_x = i as f64;
            }
        }
        assert!((peak_x - 30.0).abs() <= 2.0, "peak at {peak_x}");
    }

    #[test]
    fn fit_recovers_synthetic_p1_parameters() {
        // Noiseless p=1 spectrum: d=0.1, kmercov=30, bias=0.5, length=100k.
        let mp = ModelParams {
            d: 0.1,
            kmercov: 30.0,
            bias: 0.5,
            length: 100_000.0,
            r1: 0.0,
        };
        let mut hist = vec![0u64; 300];
        for x in 1..=300 {
            hist[x - 1] = (mp.length * predict(1, 1, &mp, 21, x as f64)).round() as u64;
        }
        let result = fit(&hist, 21, 1);
        assert!(result.converged, "model must converge");
        assert!(
            (result.kmercov - 30.0).abs() <= 2.0,
            "kmercov {} far from 30",
            result.kmercov
        );
        assert!(
            (result.length - 100_000.0).abs() / 100_000.0 < 0.1,
            "length {} far from 100k",
            result.length
        );
        assert!((result.d - 0.1).abs() < 0.05, "d {} far from 0.1", result.d);
    }
}
