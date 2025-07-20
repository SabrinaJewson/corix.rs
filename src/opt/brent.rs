//! [`Brent`]: Optimization using Brent’s method.

/// Optimization using Brent’s method.
///
/// # Examples
///
/// Estimating the maximum likelihood parameters for a half-normal distribution:
///
/// ```
/// let data = (0..50).map(|i| i as f64).collect::<Vec<_>>();
///
/// let output = corix::opt::Brent::new()
///     .maximize(|σ| {
///         if σ <= 0.0 {
///             return f64::NEG_INFINITY;
///         }
///         data.iter().map(|&x| half_norm_ln_pdf(σ, x)).sum::<f64>()
///     })
///     .unwrap();
///
/// assert_float_eq::assert_float_relative_eq!(output.x, 28.43413472232329);
/// assert_float_eq::assert_float_relative_eq!(output.f, -203.66908460766538);
/// assert_eq!(output.iters, 11);
///
/// fn half_norm_ln_pdf(σ: f64, x: f64) -> f64 {
///     -0.225_791_352_644_727_44 - σ.ln() - 0.5 * (x / σ).powi(2)
/// }
/// ```
#[derive(Debug)]
#[non_exhaustive]
pub struct Brent {
    /// The bracket in which to search.
    pub bracket: Bracket,

    /// Relative tolerance for convergence. Defaults to 1.48e-8.
    pub rel_tol: f64,

    /// The maximum number of iterations before the algorithm fails to converge.
    /// Defaults to 500.
    pub max_iters: u64,
}

/// A bracket that [`Brent`] will search in.
#[derive(Debug, Clone, Copy)]
pub enum Bracket {
    /// Compute a bracket from the given start- and end-points (uses [`opt::Bracket`]).
    Computed {
        /// The start point to pass into [`opt::Bracket`].
        start: f64,
        /// The end point to pass into [`opt::Bracket`].
        end: f64,
    },
    /// Use these values as the bracket.
    Given([f64; 3]),
}

impl Brent {
    /// Construct a new `Brent` with default parameters.
    #[must_use]
    pub fn new() -> Self {
        Self {
            bracket: Bracket::Computed {
                start: 0.0,
                end: 1.0,
            },
            // https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_optimize.py#L2441
            rel_tol: 1.48e-8,
            max_iters: 500,
        }
    }

    /// Configure the start- and end-points to pass into [`opt::Bracket`]
    /// when computing the bracket to search.
    #[must_use]
    pub fn with_computed_bracket(mut self, start: f64, end: f64) -> Self {
        self.bracket = Bracket::Computed { start, end };
        self
    }

    /// Set the three values to be used as the bracket.
    ///
    /// You will get errors if you fail to satisfy the conditions that
    /// + `a < b < c` or `c < b < a`, and
    /// + `f(b) ≤ f(a)` and `f(b) ≤ f(c)`,
    ///   where one inequality holds strictly.
    #[must_use]
    pub fn with_bracket(mut self, x: [f64; 3]) -> Self {
        self.bracket = Bracket::Given(x);
        self
    }

    crate::util::setter!(with_rel_tol => rel_tol: f64);
    crate::util::setter!(with_max_iters => max_iters: u64);

    /// Maximize the given function.
    pub fn maximize(self, mut f: impl FnMut(f64) -> f64) -> Result<Output, Error> {
        self.try_maximize(|x| Ok(f(x)))
    }

    /// Minimize the given function.
    pub fn minimize(self, mut f: impl FnMut(f64) -> f64) -> Result<Output, Error> {
        self.try_minimize(|x| Ok(f(x)))
    }

    /// Maximize the given fallible function.
    pub fn try_maximize<E>(
        self,
        mut f: impl FnMut(f64) -> Result<f64, E>,
    ) -> Result<Output, Error<E>> {
        let mut output = self.try_minimize(|x| Ok(-f(x)?))?;
        output.f = -output.f;
        Ok(output)
    }

    /// Minimize the given fallible function.
    #[expect(clippy::many_single_char_names)]
    pub fn try_minimize<E>(
        self,
        mut f: impl FnMut(f64) -> Result<f64, E>,
    ) -> Result<Output, Error<E>> {
        // https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_optimize.py#L2439-L2617

        InvalidTol::check_rel(false, self.rel_tol).map_err(Error::InvalidTol)?;

        let ([mut a, mut b, c], [_fa, fb, _fc]) = self.bracket.get(&mut f)?;

        let mut f = |arg: f64| -> Result<f64, Error<E>> {
            let res = f(arg).map_err(Error::Function)?;
            if res.is_nan() {
                return Err(Error::Nan(NanError));
            }
            Ok(res)
        };

        let [mut x, mut w, mut v] = [b; 3];
        let [mut fx, mut fw, mut fv] = [fb; 3];
        if a < c {
            b = c;
        } else {
            [a, b] = [c, a];
        }
        let mut δx = 0.0_f64;
        let mut rat = f64::NAN;

        for iters in 0..self.max_iters {
            let tol1 = self.rel_tol * x.abs() + 1e-11;
            let tol2 = 2.0 * tol1;
            let xmid = f64::midpoint(a, b);
            if (x - xmid).abs() < tol2 - 0.5 * (b - a) {
                return Ok(Output { x, f: fx, iters });
            }
            if δx.abs() <= tol1 {
                δx = if x >= xmid { a - x } else { b - x };
                rat = TWO_SUB_GOLDEN_RATIO * δx;
            } else {
                let tmp1 = (x - w) * (fx - fv);
                let tmp2 = (x - v) * (fx - fw);
                let mut p = (x - v) * tmp2 - (x - w) * tmp1;
                let tmp2 = 2.0 * (tmp2 - tmp1);
                if tmp2 > 0.0 {
                    p = -p;
                }
                let tmp2 = tmp2.abs();
                let δx_temp = δx;
                δx = rat;
                // check parabolic fit
                if p > tmp2 * (a - x)
                    && p < tmp2 * (b - x)
                    && p.abs() < (0.5 * tmp2 * δx_temp).abs()
                {
                    // if parabolic step is useful
                    rat = p / tmp2;
                    let u = x + rat;
                    if (u - a) < tol2 || (b - u) < tol2 {
                        rat = if xmid - x >= 0.0 { tol1 } else { -tol1 };
                    }
                } else {
                    δx = if x >= xmid { a - x } else { b - x };
                    rat = TWO_SUB_GOLDEN_RATIO * δx;
                }
            }

            let u = if rat.abs() < tol1 {
                if rat >= 0.0 { x + tol1 } else { x - tol1 }
            } else {
                x + rat
            };
            let fu = f(u)?;

            if fu > fx {
                if u < x {
                    a = u;
                } else {
                    b = u;
                }
                #[expect(clippy::float_cmp)]
                if fu <= fw || w == x {
                    [v, fv] = [w, fw];
                    [w, fw] = [u, fu];
                } else if fu <= fv || v == x || v == w {
                    [v, fv] = [u, fu];
                }
            } else {
                if u >= x {
                    a = x;
                } else {
                    b = x;
                }
                [v, fv] = [w, fw];
                [w, fw] = [x, fx];
                [x, fx] = [u, fu];
            }
        }

        Err(Error::Convergence(ConvergenceError {
            calls: false,
            val: self.max_iters,
        }))
    }
}

impl Default for Brent {
    fn default() -> Self {
        Self::new()
    }
}

impl Bracket {
    #[expect(clippy::many_single_char_names)]
    fn get<E>(
        self,
        mut f: impl FnMut(f64) -> Result<f64, E>,
    ) -> Result<([f64; 3], [f64; 3]), Error<E>> {
        match self {
            Self::Computed { start, end } => {
                let output = opt::Bracket::new(start, end)
                    .try_minimum(&mut f)
                    .map_err(Error::Bracket)?;
                Ok((output.x, output.f))
            }
            Self::Given(x) => {
                let mut f = |arg: f64| -> Result<f64, Error<E>> {
                    let res = f(arg).map_err(Error::Function)?;
                    if res.is_nan() {
                        return Err(Error::Nan(NanError));
                    }
                    Ok(res)
                };
                let [a, b, c] = x;
                if !(a < b && b < c || c < b && b < a) {
                    return Err(Error::InvalidBracket(InvalidBracket {
                        x,
                        kind: InvalidBracketKind::NotCollinear,
                    }));
                }
                let [fa, fb, fc] = [f(a)?, f(b)?, f(c)?];
                if !(fb < fc && fb <= fa || fb < fa && fb <= fc) {
                    return Err(Error::InvalidBracket(InvalidBracket {
                        x,
                        kind: InvalidBracketKind::NotBracket([fa, fb, fc]),
                    }));
                }
                Ok(if c < a {
                    ([c, b, a], [fc, fb, fa])
                } else {
                    ([a, b, c], [fa, fb, fc])
                })
            }
        }
    }
}

/// The output of [`Brent::minimize`].
#[derive(Debug)]
#[non_exhaustive]
pub struct Output {
    /// The `x`-value for which `f(x)` is minimized.
    pub x: f64,
    /// `f(x)`.
    pub f: f64,
    /// The number of iterations taken.
    pub iters: u64,
}

/// Errors from [`Brent`].
#[derive(Debug)]
pub enum Error<E = Infallible> {
    /// A tolerance parameter was invalid.
    InvalidTol(InvalidTol),

    /// An error caused by an invalid [`Bracket::Given`].
    InvalidBracket(InvalidBracket),

    /// An error finding a [`Bracket::Computed`].
    Bracket(opt::bracket::Error<E>),

    /// Failed to converge after a number of iterations.
    Convergence(ConvergenceError),

    /// The function returned `NaN`.
    Nan(NanError),

    /// An error from the function being optimized.
    Function(E),
}

/// An error caused by an invalid [`Bracket::Given`].
#[derive(Debug)]
#[non_exhaustive]
pub struct InvalidBracket {
    x: [f64; 3],
    kind: InvalidBracketKind,
}

#[derive(Debug)]
enum InvalidBracketKind {
    NotCollinear,
    NotBracket([f64; 3]),
}

impl Display for InvalidBracket {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "bracket points {:?} do not satisfy ", self.x)?;
        match self.kind {
            InvalidBracketKind::NotCollinear => f.write_str("a < b < c or c < b < a"),
            InvalidBracketKind::NotBracket(f_) => {
                write!(f, "f(b) < min(f(a), f(c)) (f(…) = {f_:?}")
            }
        }
    }
}

impl std::error::Error for InvalidBracket {}

// 2 - φ
const TWO_SUB_GOLDEN_RATIO: f64 = 0.381_966_0;

use super::ConvergenceError;
use super::InvalidTol;
use super::NanError;
use crate::opt;
use std::convert::Infallible;
use std::fmt;
use std::fmt::Display;
use std::fmt::Formatter;
