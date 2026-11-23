//! Root finding.
//!
//! [`brentq`] is a line-by-line port of SciPy's `brentq` (scipy/optimize/Zeros/brentq.c), including its defaults, so it takes the
//! same steps as `scipy.optimize.root_scalar(f, bracket=[a, b], method="brentq")`. Roots can differ from SciPy's in the last bit,
//! as SciPy's compiled code may fuse some multiply-adds.

/// SciPy's default absolute tolerance for `brentq`.
pub const BRENTQ_XTOL: f64 = 2e-12;

/// SciPy's default (and minimum allowed) relative tolerance for `brentq`, four times the machine epsilon.
pub const BRENTQ_RTOL: f64 = 4.0 * f64::EPSILON;

/// SciPy's default maximum number of iterations for `brentq`.
pub const BRENTQ_MAXITER: usize = 100;

/// The result of a bracketed root search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BrentqResult {
    /// The root was found to within the tolerance.
    Converged(f64),
    /// The iteration limit was reached. Holds the last estimate.
    NotConverged(f64),
    /// The function has the same sign at both ends of the bracket.
    SignError,
}

/// Finds a root of `f` in the bracket `[xa, xb]` with Brent's method.
///
/// This reproduces SciPy's implementation exactly: the same sequence of function evaluations, interpolation and bisection
/// steps, and convergence test. As in SciPy, a NaN function value is an error.
///
/// # Arguments
///
/// * `f` - The function, which may fail with an error message.
/// * `xa`, `xb` - The ends of the bracket.
/// * `xtol`, `rtol` - Absolute and relative tolerances. The root is accepted when the bracket is narrower than
///   `xtol + rtol * |x|`.
/// * `maxiter` - The maximum number of iterations.
///
/// # Returns
///
/// The outcome of the search, or an error if `f` failed or returned NaN.
pub fn brentq<F>(
    mut f: F,
    xa: f64,
    xb: f64,
    xtol: f64,
    rtol: f64,
    maxiter: usize,
) -> Result<BrentqResult, String>
where
    F: FnMut(f64) -> Result<f64, String>,
{
    let mut eval = |x: f64| -> Result<f64, String> {
        let fx = f(x)?;
        if fx.is_nan() {
            return Err(format!(
                "The function value at x={x:.6} is NaN; solver cannot converge."
            ));
        }
        Ok(fx)
    };

    let mut xpre = xa;
    let mut xcur = xb;
    let mut xblk = 0.0;
    let mut fblk = 0.0;
    let mut spre = 0.0;
    let mut scur = 0.0;

    let mut fpre = eval(xpre)?;
    let mut fcur = eval(xcur)?;
    if fpre == 0.0 {
        return Ok(BrentqResult::Converged(xpre));
    }
    if fcur == 0.0 {
        return Ok(BrentqResult::Converged(xcur));
    }
    if fpre.is_sign_negative() == fcur.is_sign_negative() {
        return Ok(BrentqResult::SignError);
    }

    for _ in 0..maxiter {
        if fpre != 0.0 && fcur != 0.0 && fpre.is_sign_negative() != fcur.is_sign_negative() {
            xblk = xpre;
            fblk = fpre;
            spre = xcur - xpre;
            scur = spre;
        }
        if fblk.abs() < fcur.abs() {
            xpre = xcur;
            xcur = xblk;
            xblk = xpre;

            fpre = fcur;
            fcur = fblk;
            fblk = fpre;
        }

        // The tolerance is 2 * delta
        let delta = (xtol + rtol * xcur.abs()) / 2.0;
        let sbis = (xblk - xcur) / 2.0;
        if fcur == 0.0 || sbis.abs() < delta {
            return Ok(BrentqResult::Converged(xcur));
        }

        if spre.abs() > delta && fcur.abs() < fpre.abs() {
            let stry = if xpre == xblk {
                // Interpolate
                -fcur * (xcur - xpre) / (fcur - fpre)
            } else {
                // Extrapolate
                let dpre = (fpre - fcur) / (xpre - xcur);
                let dblk = (fblk - fcur) / (xblk - xcur);
                -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre))
            };
            if 2.0 * stry.abs() < spre.abs().min(3.0 * sbis.abs() - delta) {
                // Good short step
                spre = scur;
                scur = stry;
            } else {
                // Bisect
                spre = sbis;
                scur = sbis;
            }
        } else {
            // Bisect
            spre = sbis;
            scur = sbis;
        }

        xpre = xcur;
        fpre = fcur;
        if scur.abs() > delta {
            xcur += scur;
        } else {
            xcur += if sbis > 0.0 { delta } else { -delta };
        }

        fcur = eval(xcur)?;
    }
    Ok(BrentqResult::NotConverged(xcur))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_simple_roots() {
        let root = |f: &dyn Fn(f64) -> f64, a, b| match brentq(
            |x| Ok(f(x)),
            a,
            b,
            BRENTQ_XTOL,
            BRENTQ_RTOL,
            BRENTQ_MAXITER,
        )
        .unwrap()
        {
            BrentqResult::Converged(x) => x,
            other => panic!("did not converge: {other:?}"),
        };
        assert!((root(&|x| x * x - 2.0, 0.0, 2.0) - 2f64.sqrt()).abs() < 1e-11);
        assert!((root(&|x| x.cos() - x, 0.0, 1.0) - 0.739_085_133_215_160_6).abs() < 1e-11);
    }

    #[test]
    fn endpoint_roots_and_sign_errors() {
        let r = brentq(
            |x| Ok(x),
            0.0,
            1.0,
            BRENTQ_XTOL,
            BRENTQ_RTOL,
            BRENTQ_MAXITER,
        )
        .unwrap();
        assert_eq!(r, BrentqResult::Converged(0.0));
        let r = brentq(
            |x| Ok(x + 5.0),
            0.0,
            1.0,
            BRENTQ_XTOL,
            BRENTQ_RTOL,
            BRENTQ_MAXITER,
        )
        .unwrap();
        assert_eq!(r, BrentqResult::SignError);
        assert!(
            brentq(
                |_| Ok(f64::NAN),
                0.0,
                1.0,
                BRENTQ_XTOL,
                BRENTQ_RTOL,
                BRENTQ_MAXITER
            )
            .is_err()
        );
    }
}
