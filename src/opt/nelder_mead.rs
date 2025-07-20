//! Optimization using Nelder–Mead. See [`NelderMead`].
//!
//! # Examples
//!
//! Fitting a GEV distribution using maximum likelihood.
//!
//! ```
//! let data = (0..50).map(|i| i as f64).collect::<Vec<_>>();
//!
//! let mean = data.iter().sum::<f64>() / (data.len() as f64);
//! let var = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (data.len() as f64 - 1.0);
//!
//! let output = corix::opt::NelderMead::new(&[mean, var.sqrt(), 0.0])
//!     .minimize(|&[μ, σ, ξ]| {
//!         if σ <= 0.0 {
//!             return f64::INFINITY;
//!         }
//!         -data.iter().map(|&x| gev_ln_pdf(μ, σ, ξ, x)).sum::<f64>()
//!     })
//!     .unwrap();
//! assert_eq!(output.iters, 139);
//! let [μ, σ, ξ] = output.x;
//! assert_float_eq::assert_float_relative_eq!(μ, 20.530_555_331_726_88);
//! assert_float_eq::assert_float_relative_eq!(σ, 15.192_413_744_842_586);
//! assert_float_eq::assert_float_relative_eq!(ξ, -0.440_690_136_761_227_8);
//! assert_float_eq::assert_float_relative_eq!(output.f, 203.122_397_429_223_53);
//!
//! fn gev_ln_pdf(μ: f64, σ: f64, ξ: f64, x: f64) -> f64 {
//!     let x = (x - μ) / σ;
//!     -σ.ln()
//!         - match ξ {
//!             0.0 => x + (-x).exp(),
//!             ξ if ξ * x < -1.0 => return -f64::INFINITY,
//!             ξ => (1.0 + ξ.recip()) * (ξ * x).ln_1p() + (1.0 + ξ * x).powf(-ξ.recip()),
//!         }
//!
//! }
//! ```

/// Optimization using Nelder–Mead.
#[derive(Debug)]
pub struct NelderMead<V: ?Sized + Vector> {
    n: V::N,
    initial_simplex: V::Simplex,

    /// Absolute tolerance in `x`-values for convergence. Defaults to 1e-4.
    pub x_abs_tol: f64,

    /// Absolute tolerance in `f(x)`-values for convergence. Defaults to 1e-4.
    pub f_abs_tol: f64,

    /// The maximum number of iterations before the algorithm fails to converge.
    /// Defaults to `N × 200`, where `N` is the number of variables.
    pub max_iters: u64,
    /// The maximum number of function calls before the algorithm fails to converge.
    /// Defaults to `u64::MAX`.
    // TODO: implement
    pub max_calls: u64,

    /// Reflection coefficient (ρ or α).
    pub reflection: f64,
    /// Expansion coefficient (χ or γ).
    pub expansion: f64,
    /// Contraction coefficient (ψ or ρ).
    pub contraction: f64,
    /// Shrink coefficient (σ).
    pub shrink: f64,
}

impl<V: ?Sized + Vector> NelderMead<V> {
    /// Construct a `NelderMead` with default parameters.
    #[must_use]
    pub fn new(initial: &V) -> Self {
        let (n_store, mut initial_simplex) = initial.to_simplex();
        let n = V::dimensions(n_store);

        for (i, chunk) in V::simplex_as_slice_mut(&mut initial_simplex)[n..]
            .chunks_exact_mut(n)
            .enumerate()
        {
            chunk[i] = match chunk[i] {
                0.0 => 0.00025,
                _ => 1.05 * chunk[i],
            };
        }

        Self::with_simplex(n_store, initial_simplex)
    }

    /// Construct a `NelderMead` with the given initial simplex.
    #[must_use]
    pub fn with_simplex(n: V::N, initial_simplex: V::Simplex) -> Self {
        let n_usize = V::dimensions(n);
        let n_f64 = n_usize as f64;

        Self {
            n,
            initial_simplex,
            // https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_optimize.py#L703
            x_abs_tol: 1e-4,
            f_abs_tol: 1e-4,
            // https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_optimize.py#L715
            max_iters: n_usize.try_into().unwrap_or(u64::MAX).saturating_mul(200),
            max_calls: u64::MAX,

            reflection: 1.0,
            expansion: 1.0 + 2.0 / n_f64,
            contraction: 0.75 - (2.0 * n_f64).recip(),
            shrink: 1.0 - n_f64.recip(),
        }
    }

    /// Configure the `NelderMead` to use standard coefficient values.
    ///
    /// By default, adaptive values will be used.
    #[must_use]
    pub fn with_standard_coefficients(mut self) -> Self {
        self.reflection = 1.0;
        self.expansion = 2.0;
        self.contraction = 0.5;
        self.shrink = 0.5;
        self
    }

    /// Set [`Self::x_abs_tol`].
    #[must_use]
    pub fn with_x_abs_tol(mut self, x_abs_tol: f64) -> Self {
        self.x_abs_tol = x_abs_tol;
        self
    }

    /// Set [`Self::f_abs_tol`].
    #[must_use]
    pub fn with_f_abs_tol(mut self, f_abs_tol: f64) -> Self {
        self.f_abs_tol = f_abs_tol;
        self
    }

    /// Set [`Self::max_iters`].
    #[must_use]
    pub fn with_max_iters(mut self, max_iters: u64) -> Self {
        self.max_iters = max_iters;
        self
    }

    /// Set [`Self::max_iters`].
    #[must_use]
    pub fn with_max_calls(mut self, max_calls: u64) -> Self {
        self.max_calls = max_calls;
        self
    }

    /// Minimize the result of the given function.
    pub fn minimize(self, mut f: impl FnMut(&V) -> f64) -> Result<Output<V>, Error> {
        self.try_minimize(|x| Ok(f(x)))
    }

    /// Minimize the result of the given fallible function.
    pub fn try_minimize<E>(
        mut self,
        mut f: impl FnMut(&V) -> Result<f64, E>,
    ) -> Result<Output<V>, Error<E>> {
        // https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_optimize.py#L700-L967

        // TODO: set bounds
        // TODO: state machine
        // TODO: big buffer
        InvalidTol::check_abs(false, self.x_abs_tol).map_err(Error::InvalidTol)?;
        InvalidTol::check_abs(true, self.f_abs_tol).map_err(Error::InvalidTol)?;

        let n = V::dimensions(self.n);
        let mut sim = V::simplex_as_slice_mut(&mut self.initial_simplex);
        let mut f = |slice: &[f64]| -> Result<f64, Error<E>> {
            let res = f(V::from_slice(slice)).map_err(Error::Function)?;
            if res.is_nan() {
                return Err(Error::Nan(NanError));
            }
            Ok(res)
        };

        let mut fsim = V::f_simplex_from_fn(self.n, |i| Ok((i, f(&sim[i * n..][..n])?)))?;
        let fsim = V::f_simplex_as_slice_mut(&mut fsim);

        let mut new_sim = V::arbitrary_simplex(self.n);
        let mut new_sim = V::simplex_as_slice_mut(&mut new_sim);

        let mut iters = 0_u64;

        loop {
            // Sort `fsim` and `sim` simulateously based on the values of `fsim`. There’s no super
            // easy way to do this; we just have to use indices.
            // PANICS: We already ensure that the slice contains no `NaN`s.
            fsim.sort_unstable_by(|&(_, a), &(_, b)| a.partial_cmp(&b).unwrap());
            for (new_i, (old_i, _)) in fsim.iter_mut().enumerate() {
                new_sim[new_i * n..][..n].copy_from_slice(&sim[*old_i * n..][..n]);
                *old_i = new_i;
            }
            (sim, new_sim) = (new_sim, sim);

            if sim[n..]
                .iter()
                .enumerate()
                .all(|(i, &param)| (param - sim[i % n]).abs() < self.x_abs_tol)
                && fsim[1..]
                    .iter()
                    .all(|(_, v)| (fsim[0].1 - v).abs() < self.f_abs_tol)
            {
                break;
            }

            if iters == self.max_iters {
                return Err(Error::Convergence(ConvergenceError(self.max_iters)));
            }
            iters += 1;

            // Use `new_sim` as scratch space. This can’t panic since `n ≥ 1`.
            let [xbar, x_reflect] = new_sim.get_disjoint_mut([0..n, n..2 * n]).unwrap();
            for i in 0..n {
                xbar[i] = sim[..n * n].chunks_exact(n).map(|c| c[i]).sum::<f64>() / (n as f64);
                x_reflect[i] = (1.0 + self.reflection) * xbar[i] - self.reflection * sim[n * n + i];
            }
            let f_reflect = f(x_reflect)?;

            let new_last = if f_reflect < fsim[0].1 {
                // expansion
                for i in 0..n {
                    xbar[i] = (1.0 + self.reflection * self.expansion) * xbar[i]
                        - self.reflection * self.expansion * sim[n * n + i];
                }
                let x_expand = xbar;
                let f_expand = f(x_expand)?;

                Some(if f_expand < f_reflect {
                    (x_expand, f_expand)
                } else {
                    (x_reflect, f_reflect)
                })
            } else if f_reflect < fsim[n - 1].1 {
                Some((x_reflect, f_reflect))
            } else if f_reflect < fsim[n].1 {
                // outside contraction
                for i in 0..n {
                    xbar[i] = (1.0 + self.contraction * self.reflection) * xbar[i]
                        - self.contraction * self.reflection * sim[n * n + i];
                }
                let x_contract = xbar;
                let f_contract = f(x_contract)?;

                (f_contract < f_reflect).then_some((x_contract, f_contract))
            } else {
                // inside contraction
                for i in 0..n {
                    xbar[i] =
                        (1.0 - self.contraction) * xbar[i] + self.contraction * sim[n * n + i];
                }
                let x_contract = xbar;
                let f_contract = f(x_contract)?;

                (f_contract < fsim[n].1).then_some((x_contract, f_contract))
            };

            match new_last {
                Some((new_x, new_f)) => {
                    sim[n * n..].copy_from_slice(new_x);
                    fsim[n].1 = new_f;
                }
                None => {
                    // do shrink
                    for j in 1..=n {
                        for i in 0..n {
                            sim[j * n + i] = sim[i] + self.shrink * (sim[j * n + i] - sim[i]);
                        }
                        fsim[j].1 = f(&sim[j * n..][..n])?;
                    }
                }
            }
        }

        Ok(Output {
            x: V::simplex_to_output(self.n, self.initial_simplex),
            f: fsim[0].1,
            iters,
        })
    }
}

pub struct Output<V: ?Sized + Vector> {
    /// The output `x`-vector.
    pub x: V::Output,
    /// `f(x)`.
    pub f: f64,
    /// The number of iterations required to achieve the result.
    /// This is always `< max_iters`.
    pub iters: u64,
}

/// An abstraction over `[f64; N]` and `[f64]`.
/// The latter uses heap allocation.
pub trait Vector {
    /// Storage of the number of dimensions.
    type N: std::fmt::Debug + Copy;

    // - 2N² + 4N 
    /// Retrieve the number of dimensions.
    fn dimensions(n: Self::N) -> usize;

    /// A list with `N(N + 1)` floats.
    type Simplex: std::fmt::Debug;

    /// Construct a `Simplex` by repeating this vector `N + 1` times.
    fn to_simplex(&self) -> (Self::N, Self::Simplex);

    /// Construct a simplex filled with arbitrary values.
    fn arbitrary_simplex(n: Self::N) -> Self::Simplex;

    /// Return the data of the `Simplex` list as a slice.
    fn simplex_as_slice_mut(simplex: &mut Self::Simplex) -> &mut [f64];

    /// Reconstruct a vector from a slice of floats.
    ///
    /// # Panics
    ///
    /// Panics if the slice has the wrong length.
    fn from_slice(vector: &[f64]) -> &Self;

    /// A list with `N + 1` `T`s.
    type FSimplex<T>;

    /// Construct an `FSimplex` with a function.
    fn f_simplex_from_fn<T, E>(
        n: Self::N,
        f: impl FnMut(usize) -> Result<T, E>,
    ) -> Result<Self::FSimplex<T>, E>;

    /// Return the data of the `FSimplex` list as a slice.
    fn f_simplex_as_slice_mut<T>(f_simplex: &mut Self::FSimplex<T>) -> &mut [T];

    /// The output of the computation.
    type Output;

    /// Construct the output of the computation from the first value of the simplex.
    fn simplex_to_output(n: Self::N, simplex: Self::Simplex) -> Self::Output;
}

impl<const N: usize> Vector for [f64; N] {
    type N = ();
    fn dimensions((): Self::N) -> usize {
        N
    }
    type Simplex = ReprCTuple<[[f64; N]; N], [f64; N]>;
    fn to_simplex(&self) -> (Self::N, Self::Simplex) {
        let simplex = ReprCTuple {
            a: [*self; N],
            b: *self,
        };
        ((), simplex)
    }
    fn arbitrary_simplex((): Self::N) -> Self::Simplex {
        ReprCTuple {
            a: [[0.0; N]; N],
            b: [0.0; N],
        }
    }
    fn simplex_as_slice_mut(simplex: &mut Self::Simplex) -> &mut [f64] {
        simplex.as_first_array_flattened_mut().as_flattened_mut()
    }
    fn from_slice(vector: &[f64]) -> &Self {
        vector.try_into().unwrap()
    }
    type FSimplex<T> = ReprCTuple<[T; N], T>;
    fn f_simplex_from_fn<T, E>(
        (): (),
        mut f: impl FnMut(usize) -> Result<T, E>,
    ) -> Result<Self::FSimplex<T>, E> {
        Ok(ReprCTuple {
            a: array_try_from_fn(&mut f)?,
            b: f(N)?,
        })
    }
    fn f_simplex_as_slice_mut<T>(f_simplex: &mut Self::FSimplex<T>) -> &mut [T] {
        f_simplex.as_first_array_flattened_mut()
    }
    type Output = [f64; N];
    fn simplex_to_output((): Self::N, simplex: Self::Simplex) -> Self::Output {
        simplex.a[0]
    }
}

impl Vector for [f64] {
    type N = usize;
    fn dimensions(n: Self::N) -> usize {
        n
    }
    type Simplex = Box<[f64]>;
    fn to_simplex(&self) -> (Self::N, Self::Simplex) {
        (self.len(), self.repeat(self.len() + 1).into_boxed_slice())
    }
    fn arbitrary_simplex(n: Self::N) -> Self::Simplex {
        vec![0.0; n * (n + 1)].into_boxed_slice()
    }
    fn simplex_as_slice_mut(simplex: &mut Self::Simplex) -> &mut [f64] {
        simplex
    }
    fn from_slice(vector: &[f64]) -> &Self {
        vector
    }
    type FSimplex<T> = Box<[T]>;
    fn f_simplex_from_fn<T, E>(
        n: usize,
        f: impl FnMut(usize) -> Result<T, E>,
    ) -> Result<Self::FSimplex<T>, E> {
        (0..=n).map(f).collect()
    }
    fn f_simplex_as_slice_mut<T>(f_simplex: &mut Self::FSimplex<T>) -> &mut [T] {
        f_simplex
    }
    type Output = Vec<f64>;
    fn simplex_to_output(n: Self::N, simplex: Self::Simplex) -> Self::Output {
        let mut vec = Vec::from(simplex);
        vec.truncate(n);
        vec
    }
}

mod hidden {}

#[derive(Debug)]
pub enum Error<E = Infallible> {
    /// A tolerance parameter was invalid.
    InvalidTol(InvalidTol),

    /// Failed to converge after a number of iterations.
    Convergence(ConvergenceError),

    /// An error from the function returning `NaN`.
    Nan(NanError),

    /// An error from the function being optimized.
    Function(E),
}

impl<E> Display for Error<E> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("error performing Nelder–Mead optimization")
    }
}

impl<E: 'static + std::error::Error> std::error::Error for Error<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidTol(e) => Some(e),
            Self::Convergence(e) => Some(e),
            Self::Nan(e) => Some(e),
            Self::Function(e) => Some(e),
        }
    }
}

use super::ConvergenceError;
use super::InvalidTol;
use crate::NanError;
use crate::util::ReprCTuple;
use crate::util::array_try_from_fn;
use std::convert::Infallible;
use std::fmt;
use std::fmt::Display;
use std::fmt::Formatter;
