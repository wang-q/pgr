//! Minpack `lmdif` Levenberg-Marquardt nonlinear least squares.
//!
//! Ported from minpack.lm (lmder.f/lmdif.f family) for use by the
//! GenomeScope fit; the algorithm is generic over any residual function.

/// Result of an `lmdif` run.
#[derive(Debug, Clone)]
pub struct LmResult {
    /// Fitted parameters.
    pub x: Vec<f64>,
    /// Residuals at the fitted parameters.
    pub fvec: Vec<f64>,
    /// Termination code (minpack info): 1/2/3 converged, 4 gtol, 5 maxfev,
    /// 6/7/8 machine precision.
    pub info: i32,
    /// Number of residual evaluations.
    pub nfev: usize,
}

/// Euclidean norm, minpack `enorm`: scaled three-sum accumulation so no
/// overflow/underflow occurs for extreme components (matches `enorm.f`).
fn enorm(x: &[f64]) -> f64 {
    const RDWARF: f64 = 3.834e-20;
    const RGIANT: f64 = 1.304e19;
    let mut s1 = 0.0;
    let mut s2 = 0.0;
    let mut s3 = 0.0;
    let mut x1max = 0.0;
    let mut x3max = 0.0;
    let agiant = RGIANT / x.len() as f64;
    for &v in x {
        let xabs = v.abs();
        if RDWARF < xabs && xabs < agiant {
            s2 += xabs * xabs;
            continue;
        }
        if xabs <= RDWARF {
            // Sum for small components.
            if xabs <= x3max {
                if xabs != 0.0 {
                    s3 += (xabs / x3max).powi(2);
                }
            } else {
                s3 = 1.0 + s3 * (x3max / xabs).powi(2);
                x3max = xabs;
            }
        } else {
            // Sum for large components.
            if xabs <= x1max {
                s1 += (xabs / x1max).powi(2);
            } else {
                s1 = 1.0 + s1 * (x1max / xabs).powi(2);
                x1max = xabs;
            }
        }
    }
    if s1 != 0.0 {
        x1max * (s1 + (s2 / x1max) / x1max).sqrt()
    } else if s2 != 0.0 {
        if s2 >= x3max {
            (s2 * (1.0 + (x3max / s2) * (x3max * s3))).sqrt()
        } else {
            (x3max * ((s2 / x3max) + (x3max * s3))).sqrt()
        }
    } else {
        x3max * s3.sqrt()
    }
}

/// Forward-difference Jacobian, minpack `fdjac2` (epsfcn = 0 uses the
/// machine precision).
pub(crate) fn fdjac2<F: Fn(&[f64]) -> Vec<f64>>(
    f: &F,
    x: &[f64],
    fvec: &[f64],
    epsfcn: f64,
    m: usize,
) -> Vec<f64> {
    let n = x.len();
    let eps = epsfcn.max(f64::EPSILON).sqrt();
    let mut fjac = vec![0.0; m * n];
    for j in 0..n {
        let temp = x[j];
        let mut h = eps * temp.abs();
        if h == 0.0 {
            h = eps;
        }
        let mut xp = x.to_vec();
        xp[j] = temp + h;
        let fp = f(&xp);
        for i in 0..m {
            fjac[j * m + i] = (fp[i] - fvec[i]) / h;
        }
    }
    fjac
}

/// QR factorization with column pivoting, minpack `qrfac` (column-major `a`).
pub(crate) fn qrfac(
    a: &mut [f64],
    m: usize,
    n: usize,
    ipvt: &mut [usize],
    rdiag: &mut [f64],
    acnorm: &mut [f64],
) {
    let epsmch = f64::EPSILON;
    let p05 = 0.5;
    for j in 0..n {
        acnorm[j] = enorm(&a[j * m..(j + 1) * m]);
        rdiag[j] = acnorm[j];
        ipvt[j] = j;
    }
    // Working copy of the column norms (minpack `wa`): swapped and updated
    // during the factorization, while `acnorm` keeps the original column
    // norms (lmdif uses them for the diag scaling and gnorm afterwards).
    let mut wa = acnorm.to_vec();
    let minmn = m.min(n);
    for j in 0..minmn {
        let mut kmax = j;
        for k in j..n {
            if rdiag[k] > rdiag[kmax] {
                kmax = k;
            }
        }
        if kmax != j {
            for i in 0..m {
                a.swap(j * m + i, kmax * m + i);
            }
            rdiag[kmax] = rdiag[j];
            wa[kmax] = wa[j];
            ipvt.swap(j, kmax);
        }
        let mut ajnorm = enorm(&a[j * m + j..j * m + m]);
        if ajnorm != 0.0 {
            if a[j * m + j] < 0.0 {
                ajnorm = -ajnorm;
            }
            for i in j..m {
                a[j * m + i] /= ajnorm;
            }
            a[j * m + j] += 1.0;
            for k in j + 1..n {
                let mut sum = 0.0;
                for i in j..m {
                    sum += a[j * m + i] * a[k * m + i];
                }
                let temp = sum / a[j * m + j];
                for i in j..m {
                    a[k * m + i] -= temp * a[j * m + i];
                }
                if rdiag[k] != 0.0 {
                    let t = a[k * m + j] / rdiag[k];
                    rdiag[k] *= (0.0f64.max(1.0 - t * t)).sqrt();
                    if p05 * (rdiag[k] / wa[k]).powi(2) <= epsmch {
                        rdiag[k] = enorm(&a[k * m + j + 1..k * m + m]);
                        wa[k] = rdiag[k];
                    }
                }
            }
        }
        rdiag[j] = -ajnorm;
    }
}

/// Solve the damped triangular system, minpack `qrsolv`.
#[allow(clippy::too_many_arguments)]
fn qrsolv(
    n: usize,
    r: &mut [f64],
    ipvt: &[usize],
    diag: &[f64],
    qtb: &[f64],
    x: &mut [f64],
    sdiag: &mut [f64],
    wa: &mut [f64],
) {
    for j in 0..n {
        for i in (j + 1)..n {
            // r(i,j) = r(j,i): copy the strict upper triangle into the
            // lower triangle (column-major: upper is r[i*n+j]).
            let tmp = r[i * n + j];
            r[j * n + i] = tmp;
        }
        x[j] = r[j * n + j];
        wa[j] = qtb[j];
    }
    for j in 0..n {
        let l = ipvt[j];
        if diag[l] == 0.0 {
            continue;
        }
        for v in sdiag[j..n].iter_mut() {
            *v = 0.0;
        }
        sdiag[j] = diag[l];
        let mut qtbpj = 0.0;
        for k in j..n {
            if sdiag[k] == 0.0 {
                continue;
            }
            let (sin, cos);
            if r[k * n + k].abs() >= sdiag[k].abs() {
                let tan = sdiag[k] / r[k * n + k];
                cos = 0.5 / (0.25 + 0.25 * tan * tan).sqrt();
                sin = cos * tan;
            } else {
                let cotan = r[k * n + k] / sdiag[k];
                sin = 0.5 / (0.25 + 0.25 * cotan * cotan).sqrt();
                cos = sin * cotan;
            }
            r[k * n + k] = cos * r[k * n + k] + sin * sdiag[k];
            let temp = cos * wa[k] + sin * qtbpj;
            qtbpj = -sin * wa[k] + cos * qtbpj;
            wa[k] = temp;
            for i in (k + 1)..n {
                // r(i,k) = cos*r(i,k) + sin*sdiag(i): lower triangle
                // (row i, col k) is r[k*n+i] in column-major order.
                let temp = cos * r[k * n + i] + sin * sdiag[i];
                sdiag[i] = -sin * r[k * n + i] + cos * sdiag[i];
                r[k * n + i] = temp;
            }
        }
        sdiag[j] = r[j * n + j];
        r[j * n + j] = x[j];
    }
    let mut nsing = n;
    for j in 0..n {
        if sdiag[j] == 0.0 && nsing == n {
            nsing = j;
        }
        if nsing < n {
            wa[j] = 0.0;
        }
    }
    if nsing >= 1 {
        for kk in 0..nsing {
            let j = nsing - kk - 1;
            let mut sum = 0.0;
            for i in (j + 1)..nsing {
                // r(i,j) is at r[j*n+i] (column-major).
                sum += r[j * n + i] * wa[i];
            }
            wa[j] = (wa[j] - sum) / sdiag[j];
        }
    }
    for j in 0..n {
        let l = ipvt[j];
        x[l] = wa[j];
    }
}

/// Determine the Levenberg-Marquardt parameter, minpack `lmpar`.
#[allow(clippy::too_many_arguments)]
fn lmpar(
    n: usize,
    r: &mut [f64],
    ipvt: &[usize],
    diag: &[f64],
    qtb: &[f64],
    delta: f64,
    mut par: f64,
    x: &mut [f64],
    sdiag: &mut [f64],
    wa1: &mut [f64],
    wa2: &mut [f64],
) -> f64 {
    let dwarf = f64::MIN_POSITIVE;
    let mut nsing = n;
    for j in 0..n {
        wa1[j] = qtb[j];
        if r[j * n + j] == 0.0 && nsing == n {
            nsing = j;
        }
        if nsing < n {
            wa1[j] = 0.0;
        }
    }
    if nsing >= 1 {
        for kk in 0..nsing {
            let j = nsing - kk - 1;
            wa1[j] /= r[j * n + j];
            let temp = wa1[j];
            for i in 0..j {
                wa1[i] -= r[j * n + i] * temp;
            }
        }
    }
    for j in 0..n {
        let l = ipvt[j];
        x[l] = wa1[j];
    }
    let mut iter = 0;
    for j in 0..n {
        wa2[j] = diag[j] * x[j];
    }
    let mut dxnorm = enorm(&wa2[..n]);
    let mut fp = dxnorm - delta;
    if fp <= 0.1 * delta {
        return 0.0;
    }
    let mut parl = 0.0;
    if nsing == n {
        for j in 0..n {
            let l = ipvt[j];
            wa1[j] = diag[l] * (wa2[l] / dxnorm);
        }
        for j in 0..n {
            let mut sum = 0.0;
            for i in 0..j {
                sum += r[j * n + i] * wa1[i];
            }
            wa1[j] = (wa1[j] - sum) / r[j * n + j];
        }
        let temp = enorm(&wa1[..n]);
        parl = ((fp / delta) / temp) / temp;
    }
    for j in 0..n {
        let mut sum = 0.0;
        for i in 0..=j {
            sum += r[j * n + i] * qtb[i];
        }
        let l = ipvt[j];
        wa1[j] = sum / diag[l];
    }
    let gnorm = enorm(&wa1[..n]);
    let mut paru = gnorm / delta;
    if paru == 0.0 {
        paru = dwarf / delta.min(0.1);
    }
    par = par.max(parl).min(paru);
    if par == 0.0 {
        par = gnorm / dxnorm;
    }
    loop {
        iter += 1;
        if par == 0.0 {
            par = dwarf.max(0.001 * paru);
        }
        let temp = par.sqrt();
        for j in 0..n {
            wa1[j] = temp * diag[j];
        }
        qrsolv(n, r, ipvt, &wa1[..n], qtb, x, sdiag, wa2);
        for j in 0..n {
            wa2[j] = diag[j] * x[j];
        }
        dxnorm = enorm(&wa2[..n]);
        let temp2 = fp;
        fp = dxnorm - delta;
        if fp.abs() <= 0.1 * delta || (parl == 0.0 && fp <= temp2 && temp2 < 0.0) || iter == 10 {
            break;
        }
        for j in 0..n {
            let l = ipvt[j];
            wa1[j] = diag[l] * (wa2[l] / dxnorm);
        }
        for j in 0..n {
            wa1[j] /= sdiag[j];
            let temp3 = wa1[j];
            for i in (j + 1)..n {
                // r(i,j) is at r[j*n+i] (column-major).
                wa1[i] -= r[j * n + i] * temp3;
            }
        }
        let temp3 = enorm(&wa1[..n]);
        let parc = ((fp / delta) / temp3) / temp3;
        if fp > 0.0 {
            parl = parl.max(par);
        }
        if fp < 0.0 {
            paru = paru.min(par);
        }
        par = parl.max(par + parc);
    }
    if iter == 0 {
        par = 0.0;
    }
    par
}

/// Minpack `lmdif`: Levenberg-Marquardt with forward-difference Jacobian.
///
/// Parameters are projected into `[lower, upper]` at every evaluation
/// (R's nls.lm fcn wrapper clamps), `factor` scales the initial step bound,
/// and `info` follows minpack: 1 (ftol), 2 (xtol), 3 (both), 4 (gtol),
/// 5 (maxfev), 6/7/8 (machine precision).
#[allow(clippy::too_many_arguments)]
pub fn lmdif<F: Fn(&[f64]) -> Vec<f64>>(
    f: &F,
    x0: &[f64],
    lower: &[f64],
    upper: &[f64],
    ftol: f64,
    ptol: f64,
    gtol: f64,
    maxfev: usize,
    epsfcn: f64,
    factor: f64,
) -> LmResult {
    let n = x0.len();
    let mut x = x0.to_vec();
    for j in 0..n {
        x[j] = x[j].clamp(lower[j], upper[j]);
    }
    let mut fvec = f(&x);
    let m = fvec.len();
    let epsmch = f64::EPSILON;
    let (one, p1, p5, p25, p75, p0001) = (1.0, 0.1, 0.5, 0.25, 0.75, 0.0001);

    let mut fjac = vec![0.0; m * n];
    let mut diag = vec![0.0; n];
    let mut qtf = vec![0.0; n];
    let mut wa1 = vec![0.0; n];
    let mut wa2 = vec![0.0; n];
    let mut wa3 = vec![0.0; n];
    let mut wa4 = vec![0.0; m];
    let mut ipvt = vec![0usize; n];
    let mut acnorm = vec![0.0; n];
    let mut rdiag = vec![0.0; n];
    let mut sdiag = vec![0.0; n];

    let mut nfev = 1;
    let mut fnorm = enorm(&fvec);
    let mut par = 0.0;
    let mut delta = 0.0;
    let mut xnorm = 0.0;
    let mut info = 0i32;
    let mut iter = 1;

    'outer: loop {
        let fj = fdjac2(
            &|p: &[f64]| {
                let mut pc = p.to_vec();
                for j in 0..n {
                    pc[j] = pc[j].clamp(lower[j], upper[j]);
                }
                f(&pc)
            },
            &x,
            &fvec,
            epsfcn,
            m,
        );
        fjac.copy_from_slice(&fj);
        qrfac(&mut fjac, m, n, &mut ipvt, &mut rdiag, &mut acnorm);

        if iter == 1 {
            for j in 0..n {
                diag[j] = acnorm[j];
                if diag[j] == 0.0 {
                    diag[j] = 1.0;
                }
            }
            for j in 0..n {
                wa3[j] = diag[j] * x[j];
            }
            xnorm = enorm(&wa3[..n]);
            delta = factor * xnorm;
            if delta == 0.0 {
                delta = factor;
            }
        }

        wa4[..m].copy_from_slice(&fvec);
        for j in 0..n {
            if fjac[j * m + j] != 0.0 {
                let mut sum = 0.0;
                for i in j..m {
                    sum += fjac[j * m + i] * wa4[i];
                }
                let temp = -sum / fjac[j * m + j];
                for i in j..m {
                    wa4[i] += fjac[j * m + i] * temp;
                }
            }
            fjac[j * m + j] = rdiag[j];
            qtf[j] = wa4[j];
        }

        let mut gnorm = 0.0f64;
        if fnorm != 0.0 {
            for j in 0..n {
                let l = ipvt[j];
                if acnorm[l] == 0.0 {
                    continue;
                }
                let mut sum = 0.0;
                for i in 0..=j {
                    sum += fjac[j * m + i] * (qtf[i] / fnorm);
                }
                gnorm = gnorm.max((sum / acnorm[l]).abs());
            }
        }
        if gnorm <= gtol {
            info = 4;
            break;
        }
        for j in 0..n {
            diag[j] = diag[j].max(acnorm[j]);
        }

        let mut rmat = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                rmat[j * n + i] = fjac[j * m + i];
            }
        }

        'inner: loop {
            let new_par = lmpar(
                n, &mut rmat, &ipvt, &diag, &qtf, delta, par, &mut wa1, &mut sdiag, &mut wa2,
                &mut wa3,
            );
            par = new_par;
            for j in 0..n {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            let pnorm = enorm(&wa3[..n]);
            if iter == 1 {
                delta = delta.min(pnorm);
            }
            let mut trial = wa2[..n].to_vec();
            for j in 0..n {
                trial[j] = trial[j].clamp(lower[j], upper[j]);
            }
            let fvec1 = f(&trial);
            nfev += 1;
            let fnorm1 = enorm(&fvec1);
            let actred = if p1 * fnorm1 < fnorm {
                one - (fnorm1 / fnorm).powi(2)
            } else {
                -one
            };
            for j in 0..n {
                wa3[j] = 0.0;
                let l = ipvt[j];
                let temp = wa1[l];
                for i in 0..=j {
                    wa3[i] += rmat[i * n + j] * temp;
                }
            }
            let temp1 = enorm(&wa3[..n]) / fnorm;
            let temp2 = (par.sqrt() * pnorm) / fnorm;
            let prered = temp1 * temp1 + temp2 * temp2 / p5;
            let dirder = -(temp1 * temp1 + temp2 * temp2);
            let ratio = if prered != 0.0 { actred / prered } else { 0.0 };

            if ratio <= p25 {
                let temp = if actred >= 0.0 {
                    p5
                } else {
                    p5 * dirder / (dirder + p5 * actred)
                };
                let temp = if p1 * fnorm1 >= fnorm || temp < p1 {
                    p1
                } else {
                    temp
                };
                delta = temp * delta.min(pnorm / p1);
                par /= temp;
            } else if par == 0.0 || ratio >= p75 {
                delta = pnorm / p5;
                par *= p5;
            }

            if ratio >= p0001 {
                for j in 0..n {
                    x[j] = trial[j];
                    wa2[j] = diag[j] * x[j];
                }
                fvec.copy_from_slice(&fvec1);
                xnorm = enorm(&wa2[..n]);
                fnorm = fnorm1;
                iter += 1;
            }

            if actred.abs() <= ftol && prered <= ftol && p5 * ratio <= one {
                info = 1;
            }
            if delta <= ptol * xnorm {
                info = if info == 1 { 3 } else { 2 };
            }
            if info != 0 {
                break 'outer;
            }
            if nfev >= maxfev {
                info = 5;
                break 'outer;
            }
            if actred.abs() <= epsmch && prered <= epsmch && p5 * ratio <= one {
                info = 6;
                break 'outer;
            }
            if delta <= epsmch * xnorm {
                info = 7;
                break 'outer;
            }
            if gnorm <= epsmch {
                info = 8;
                break 'outer;
            }
            if ratio < p0001 {
                continue 'inner;
            }
            break 'inner;
        }
    }
    LmResult {
        x,
        fvec,
        info,
        nfev,
    }
}

/// Invert a small square matrix by Gauss-Jordan elimination.
pub(crate) fn invert_matrix(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut aug = vec![0.0; n * 2 * n];
    for i in 0..n {
        for j in 0..n {
            aug[i * 2 * n + j] = a[i * n + j];
        }
        aug[i * 2 * n + n + i] = 1.0;
    }
    for col in 0..n {
        let mut piv = col;
        for r in col + 1..n {
            if aug[r * 2 * n + col].abs() > aug[piv * 2 * n + col].abs() {
                piv = r;
            }
        }
        if aug[piv * 2 * n + col].abs() < 1e-300 {
            return None;
        }
        if piv != col {
            for c in 0..2 * n {
                aug.swap(col * 2 * n + c, piv * 2 * n + c);
            }
        }
        let d = aug[col * 2 * n + col];
        for c in 0..2 * n {
            aug[col * 2 * n + c] /= d;
        }
        for r in 0..n {
            if r != col {
                let f = aug[r * 2 * n + col];
                for c in 0..2 * n {
                    aug[r * 2 * n + c] -= f * aug[col * 2 * n + c];
                }
            }
        }
    }
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            inv[i * n + j] = aug[i * 2 * n + n + j];
        }
    }
    Some(inv)
}

#[cfg(test)]
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;

    #[test]
    fn lmdif_fits_linear_model() {
        // y = 2x + 3 with noise-free points; fit (slope, intercept).
        let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| 2.0 * xi + 3.0).collect();
        let residual = |p: &[f64]| -> Vec<f64> {
            x.iter()
                .zip(&y)
                .map(|(&xi, &yi)| yi - (p[0] * xi + p[1]))
                .collect()
        };
        let res = lmdif(
            &residual,
            &[0.0, 0.0],
            &[f64::NEG_INFINITY, f64::NEG_INFINITY],
            &[f64::INFINITY, f64::INFINITY],
            1.49e-8,
            1.49e-8,
            0.0,
            300,
            0.0,
            0.1,
        );
        assert!((1..=8).contains(&res.info), "info {}", res.info);
        assert!((res.x[0] - 2.0).abs() < 1e-6, "slope {}", res.x[0]);
        assert!((res.x[1] - 3.0).abs() < 1e-6, "intercept {}", res.x[1]);
    }

    #[test]
    fn lmdif_respects_bounds() {
        // Start near the solution so the trust region is well scaled.
        let residual = |p: &[f64]| vec![p[0] - 1.0, p[1] - 2.0];
        let res = lmdif(
            &residual,
            &[1.5, 2.5],
            &[0.0, 0.0],
            &[3.0, f64::INFINITY],
            1.49e-8,
            1.49e-8,
            0.0,
            300,
            0.0,
            0.1,
        );
        assert!(res.x[0] <= 3.0 + 1e-6, "x0 must respect upper bound");
        assert!((res.x[1] - 2.0).abs() < 1e-6, "x1 {}", res.x[1]);
    }

    #[test]
    fn lmpar_matches_minpack_lm_reference() {
        // First-iteration inputs from the minpack.lm Fortran run on the
        // GenomeScope problem; expected output par=0.56244753664 and
        // x=(0.095564939, 4.656859496, 0.679095862, -90.452266883).
        let n = 4;
        let r = vec![
            -7955.050625225653,
            0.0,
            0.0,
            0.0, //
            -1614.2909977540087,
            -1218.6343721288254,
            0.0,
            0.0, //
            -68.334174844890711,
            285.81791036053585,
            -307.59028780530474,
            0.0, //
            8.403675929025966,
            8.500842227421385e-6,
            -1.8756556049487534e-4,
            -1.0352689913445627e-4,
        ];
        let ipvt = vec![0usize, 2, 1, 3];
        let diag = vec![
            7955.050625225653,
            425.40947625405056,
            2022.6233357606282,
            8.403675931761132,
        ];
        let qtb = vec![
            -3362.422530966394,
            -212.386257577499,
            -3543.6065262283173,
            -451.94132229970836,
        ];
        let delta = 2613.2047615359565;
        let mut x = vec![0.0; n];
        let mut sdiag = vec![0.0; n];
        let mut wa1 = vec![0.0; n];
        let mut wa2 = vec![0.0; n];
        let par: f64 = lmpar(
            n,
            &mut r.clone(),
            &ipvt,
            &diag,
            &qtb,
            delta,
            0.0,
            &mut x,
            &mut sdiag,
            &mut wa1,
            &mut wa2,
        );
        assert!((par - 0.56244753664351843).abs() < 1e-12, "par {par}");
        let expected = [
            0.0955649391598923,
            4.656859495638244,
            0.6790958621285849,
            -90.45226688323875,
        ];
        for (i, e) in expected.iter().enumerate() {
            assert!((x[i] - e).abs() < 1e-8, "x[{i}] {} vs {e}", x[i]);
        }
    }

    #[test]
    fn qrsolv_matches_minpack_lm_reference() {
        // First lmpar iteration of the GenomeScope run: par=1.5063873743616401e-6.
        let n = 4;
        let mut r = vec![
            -7955.050625225653,
            0.0,
            0.0,
            0.0, //
            -1614.2909977540087,
            -1218.6343721288254,
            0.0,
            0.0, //
            -68.334174844890711,
            285.81791036053585,
            -307.59028780530474,
            0.0, //
            8.403675929025966,
            8.500842227421385e-6,
            -1.8756556049487534e-4,
            -1.0352689913445627e-4,
        ];
        let ipvt = [0usize, 2, 1, 3];
        let orig_diag = [
            7955.050625225653,
            425.40947625405056,
            2022.6233357606282,
            8.403675931761132,
        ];
        let par: f64 = 1.5063873743616401e-6;
        let temp = par.sqrt();
        let diag: Vec<f64> = orig_diag.iter().map(|d| temp * d).collect();
        let qtb = vec![
            -3362.422530966394,
            -212.386257577499,
            -3543.6065262283173,
            -451.94132229970836,
        ];
        let mut x = vec![0.0; n];
        let mut sdiag = vec![0.0; n];
        let mut wa = vec![0.0; n];
        qrsolv(n, &mut r, &ipvt, &diag, &qtb, &mut x, &mut sdiag, &mut wa);
        let expected = [
            0.10233276198364068f64,
            11.520260552142499,
            2.8762302542552676,
            342.93815320894311,
        ];
        for (i, e) in expected.iter().enumerate() {
            assert!((x[i] - e).abs() < 1e-8f64, "x[{i}] {} vs {e}", x[i]);
        }
        let expected_s = [
            -7955.056616917308f64,
            -1218.6385112553457,
            -307.5917711573897,
            0.014587118848163102,
        ];
        for (i, e) in expected_s.iter().enumerate() {
            assert!(
                (sdiag[i] - e).abs() < 1e-8f64,
                "sdiag[{i}] {} vs {e}",
                sdiag[i]
            );
        }
    }

    #[test]
    fn invert_matrix_recovers_inverse() {
        let a = vec![4.0, 7.0, 2.0, 6.0];
        let inv = invert_matrix(&a, 2).unwrap();
        // A * A^-1 = I
        let prod = |i: usize, j: usize| a[i * 2] * inv[j] + a[i * 2 + 1] * inv[2 + j];
        assert!((prod(0, 0) - 1.0).abs() < 1e-10);
        assert!((prod(0, 1) - 0.0).abs() < 1e-10);
        assert!((prod(1, 0) - 0.0).abs() < 1e-10);
        assert!((prod(1, 1) - 1.0).abs() < 1e-10);
    }
}
