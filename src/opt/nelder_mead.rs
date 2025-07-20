//! [`NelderMead`]: Optimization using Nelder–Mead.

/// Optimization using Nelder–Mead.
///
/// [`Self::begin_minimize`] is the lowest-level interface for advanced use cases,
/// but we provide convenient wrappers too.
///
/// # Examples
///
/// Fitting a GEV distribution using maximum likelihood.
///
/// ```
/// let data = (0..50).map(|i| i as f64).collect::<Vec<_>>();
///
/// let mean = data.iter().sum::<f64>() / (data.len() as f64);
/// let var = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (data.len() as f64 - 1.0);
///
/// let output = corix::opt::NelderMead::new(&[mean, var.sqrt(), 0.0])
///     .with_lower_bound([f64::NEG_INFINITY, f64::EPSILON, f64::NEG_INFINITY])
///     .maximize(|&[μ, σ, ξ]| data.iter().map(|&x| gev_ln_pdf(μ, σ, ξ, x)).sum::<f64>())
///     .unwrap();
/// let [μ, σ, ξ] = output.x;
/// assert_float_eq::assert_float_relative_eq!(μ, 20.530_555_331_726_88);
/// assert_float_eq::assert_float_relative_eq!(σ, 15.192_413_744_842_586);
/// assert_float_eq::assert_float_relative_eq!(ξ, -0.440_690_136_761_227_8);
/// assert_float_eq::assert_float_relative_eq!(output.f, -203.122_397_429_223_53);
/// assert_eq!(output.iters, 139);
/// assert_eq!(output.calls, 250);
///
/// fn gev_ln_pdf(μ: f64, σ: f64, ξ: f64, x: f64) -> f64 {
///     let x = (x - μ) / σ;
///     -σ.ln()
///         - match ξ {
///             0.0 => x + (-x).exp(),
///             ξ if ξ * x < -1.0 => return f64::NEG_INFINITY,
///             ξ => (1.0 + ξ.recip()) * (ξ * x).ln_1p() + (1.0 + ξ * x).powf(-ξ.recip()),
///         }
///
/// }
/// ```
pub struct NelderMead<V: ?Sized + Vector> {
    n: usize,
    // The buffer itself; of the form [f64; 2n² + 4n + 2]. It contains:
    //
    // - [[f64; n]; n + 1]: The current simplex
    // - [[f64; n]; n + 1]: Alternative storage or scratch space
    // - [(u64, f64); n + 1]: f(current simplex) and indices (used for sorting)
    buf: V::Buf,

    /// The lower bound of the accepted input.
    ///
    /// Everything in the simplex, including the initial guess,
    /// to this bound.
    pub lower_bound: V::Owned,
    /// The upper bound of the accepted input.
    ///
    /// Everything in the simplex, including the initial guess,
    /// will be silently clipped to this bound.
    pub upper_bound: V::Owned,

    /// Absolute tolerance in `x`-values for convergence. Defaults to 1e-4.
    pub x_abs_tol: f64,

    /// Absolute tolerance in `f(x)`-values for convergence. Defaults to 1e-4.
    pub f_abs_tol: f64,

    /// The maximum number of iterations before the algorithm fails to converge.
    /// Defaults to `N × 200`, where `N` is the number of variables.
    pub max_iters: u64,
    /// The maximum number of function calls before the algorithm fails to converge.
    /// Defaults to `u64::MAX`.
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
    /// Construct a `NelderMead` with default parameters, given an initial guess.
    #[must_use]
    pub fn new(initial: &V) -> Self {
        Self::with_simplex(initial, |i, v| {
            let x = &mut v.borrow_mut()[i];
            *x = match *x {
                0.0 => 0.00025,
                x => 1.05 * x,
            };
        })
    }

    /// Construct a `NelderMead` using a function to generate the initial simplex.
    ///
    /// The function is called `n` times, with values from `0..n`.
    #[must_use]
    pub fn with_simplex(initial: &V, mut f: impl FnMut(usize, &mut V)) -> Self {
        let n = initial.borrow().len();
        let mut buf = V::allocate(n);

        let ([initial_sim, _], _) = split_buf_mut::<V>(n, &mut buf);
        crate::util::write_copy_of_slice(&mut initial_sim[..n], initial.borrow());
        for (i, chunk) in initial_sim[n..].chunks_exact_mut(n).enumerate() {
            let chunk = crate::util::write_copy_of_slice(chunk, initial.borrow());
            f(i, V::from_slice_mut(chunk));
        }

        let n_f64 = n as f64;
        Self {
            n,
            buf,

            lower_bound: V::unbounded(false),
            upper_bound: V::unbounded(true),

            // https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_optimize.py#L703
            x_abs_tol: 1e-4,
            f_abs_tol: 1e-4,
            // https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_optimize.py#L715
            max_iters: n.try_into().unwrap_or(u64::MAX).saturating_mul(200),
            max_calls: u64::MAX,

            reflection: 1.0,
            expansion: 1.0 + 2.0 / n_f64,
            contraction: 0.75 - (2.0 * n_f64).recip(),
            shrink: 1.0 - n_f64.recip(),
        }
    }

    /// Get the initial simplex.
    ///
    /// This is not yet clipped to the bounds; that will happen later.
    #[must_use]
    pub fn initial_simplex(&self) -> impl Clone + ExactSizeIterator<Item = &V> {
        let ([initial_sim, _], _) = split_buf::<V>(self.n, &self.buf);
        // SAFETY: Initialized by `with_simplex`.
        let initial_sim = unsafe { crate::util::slice_assume_init(initial_sim) };
        initial_sim.chunks_exact(self.n).map(V::from_slice)
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

    crate::util::setter!(with_lower_bound => lower_bound: V::Owned);
    crate::util::setter!(with_upper_bound => upper_bound: V::Owned);
    crate::util::setter!(with_x_abs_tol => x_abs_tol: f64);
    crate::util::setter!(with_f_abs_tol => f_abs_tol: f64);
    crate::util::setter!(with_max_iters => max_iters: u64);
    crate::util::setter!(with_max_calls => max_calls: u64);

    /// Maximize the result of the given function.
    pub fn maximize(self, mut f: impl FnMut(&V) -> f64) -> Result<Output<V>, Error> {
        self.try_maximize(|x| Ok(f(x)))
    }

    /// Minimize the result of the given function.
    pub fn minimize(self, mut f: impl FnMut(&V) -> f64) -> Result<Output<V>, Error> {
        self.try_minimize(|x| Ok(f(x)))
    }

    /// Maximize the result of the given fallible function.
    pub fn try_maximize<E>(
        self,
        mut f: impl FnMut(&V) -> Result<f64, E>,
    ) -> Result<Output<V>, Error<E>> {
        let mut output = self.try_minimize(|x| Ok(-f(x)?))?;
        output.f = -output.f;
        Ok(output)
    }

    /// Minimize the result of the given fallible function.
    pub fn try_minimize<E>(
        self,
        f: impl FnMut(&V) -> Result<f64, E>,
    ) -> Result<Output<V>, Error<E>> {
        let mut minimizer = self.begin_minimize(f)?;
        while !minimizer.step()? {}
        Ok(minimizer.into_output())
    }

    /// Begin a minimization operation. Returns a [`Minimizer`] that may be used to step through
    /// minimization.
    pub fn begin_minimize<E, F: FnMut(&V) -> Result<f64, E>>(
        mut self,
        mut function: F,
    ) -> Result<Minimizer<V, F>, Error<E>> {
        // https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_optimize.py#L700-L967

        InvalidTol::check_abs(false, self.x_abs_tol).map_err(Error::InvalidTol)?;
        InvalidTol::check_abs(true, self.f_abs_tol).map_err(Error::InvalidTol)?;

        let lower = self.lower_bound.borrow().borrow();
        let upper = self.upper_bound.borrow().borrow();
        for (i, (&lower, &upper)) in lower.iter().zip(upper).enumerate() {
            if !(lower < upper) {
                return Err(Error::InvalidBound(InvalidBound { i, lower, upper }));
            }
        }

        let n = self.n;

        let ([initial_sim, scratch_space], fsim) = split_buf_mut::<V>(n, &mut self.buf);
        // SAFETY: Initialized by `with_simplex`.
        let initial_sim = unsafe { crate::util::slice_assume_init_mut(initial_sim) };

        for (i, initial) in initial_sim.chunks_exact_mut(n).enumerate() {
            for (k, initial) in initial.iter_mut().enumerate() {
                let upper = upper.get(k).copied().unwrap_or(f64::INFINITY);
                let lower = lower.get(k).copied().unwrap_or(f64::NEG_INFINITY);
                if i != 0 && *initial > upper {
                    *initial = 2.0 * upper - *initial;
                }
                *initial = initial.clamp(lower, upper);
            }
        }

        // Wrapper around the function call, keeping track of NaNs.
        let mut f = |slice: &[f64]| -> Result<f64, Error<E>> {
            let res = function(V::from_slice(slice)).map_err(Error::Function)?;
            if res.is_nan() {
                return Err(Error::Nan(NanError));
            }
            Ok(res)
        };

        // Initializing `fsim`.
        let calls = (n + 1) as u64;
        if self.max_calls < calls {
            return Err(Error::Convergence(ConvergenceError {
                calls: true,
                val: self.max_calls,
            }));
        }
        let fsim = crate::util::slice_try_init_with(fsim, |i| {
            Ok(Pair {
                index: i as u64,
                f: f(&initial_sim[i * n..][..n])?,
            })
        })?;

        // Sort `fsim` and `sim` simulateously based on the values of `fsim`. There’s no super
        // easy way to do this; we just have to use indices.
        // PANICS: We already ensure that the slice contains no `NaN`s.
        fsim.sort_unstable_by(|&a, &b| a.f.partial_cmp(&b.f).unwrap());
        for (new_i, Pair { index: old_i, .. }) in fsim.iter_mut().enumerate() {
            #[expect(clippy::cast_possible_truncation, reason = "came from a usize")]
            let old_i_usize = *old_i as usize;
            crate::util::write_copy_of_slice(
                &mut scratch_space[new_i * n..][..n],
                &initial_sim[old_i_usize * n..][..n],
            );
            *old_i = new_i as u64;
        }

        // Premultiply.
        self.expansion *= self.reflection;

        Ok(Minimizer {
            config: self,
            function,
            iters: 0,
            calls,
        })
    }
}

impl<V: ?Sized + Vector + Debug> Debug for NelderMead<V> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut f = f.debug_struct("NelderMead");
        f.field("initial_sim", &DebugIter(self.initial_simplex()));
        config_fields(&mut f, self);
        f.finish()
    }
}

fn config_fields<V: ?Sized + Vector + Debug>(
    f: &mut fmt::DebugStruct<'_, '_>,
    config: &NelderMead<V>,
) {
    f.field("n", &config.n)
        .field("lower_bound", &config.lower_bound.borrow())
        .field("upper_bound", &config.upper_bound.borrow())
        .field("x_abs_tol", &config.x_abs_tol)
        .field("f_abs_tol", &config.f_abs_tol)
        .field("max_iters", &config.max_iters)
        .field("max_calls", &config.max_calls)
        .field("reflection", &config.reflection)
        .field("expansion", &config.expansion)
        .field("contraction", &config.contraction)
        .field("shrink", &config.shrink);
}

/// A state machine that drives [`NelderMead`] minimization.
#[non_exhaustive]
pub struct Minimizer<V: ?Sized + Vector, F> {
    config: NelderMead<V>,
    function: F,
    iters: u64,
    calls: u64,
}

impl<V: ?Sized + Vector, F: FnMut(&V) -> Result<f64, E>, E> Minimizer<V, F> {
    /// Run a step of the Nelder–Mead algorithm.
    ///
    /// Returns `Ok(true)` if minimization completed successfully. Returns `Ok(false)` if more
    /// steps should be run.
    ///
    /// Calling this function after it has returned `Ok(true)` will just continue to return
    /// `Ok(true)`.
    #[expect(clippy::too_many_lines)]
    pub fn step(&mut self) -> Result<bool, Error<E>> {
        let config = &mut self.config;
        let n = config.n;

        // Split up and prepare the buffer.
        let ([sim_a, sim_b], fsim) = split_buf_mut::<V>(n, &mut config.buf);
        // SAFETY: Initialized by `begin_minimize`.
        let sim_a = unsafe { crate::util::slice_assume_init_mut(sim_a) };
        let sim_b = unsafe { crate::util::slice_assume_init_mut(sim_b) };
        let fsim = unsafe { crate::util::slice_assume_init_mut(fsim) };

        let (sim, scratch_space) = match self.iters % 2 != 0 {
            false => (sim_b, sim_a),
            true => (sim_a, sim_b),
        };

        // Wrapper around the function call, keeping track of maximum calls and NaNs.
        let mut f = |slice: &mut [f64]| -> Result<f64, Error<E>> {
            if self.calls == config.max_calls {
                return Err(Error::Convergence(ConvergenceError {
                    calls: true,
                    val: self.calls,
                }));
            }
            self.calls += 1;

            let lower = config.lower_bound.borrow().borrow();
            let upper = config.upper_bound.borrow().borrow();
            slice
                .iter_mut()
                .zip(lower)
                .for_each(|(x, &b)| *x = x.max(b));
            slice
                .iter_mut()
                .zip(upper)
                .for_each(|(x, &b)| *x = x.min(b));

            let res = (self.function)(V::from_slice(slice)).map_err(Error::Function)?;
            if res.is_nan() {
                return Err(Error::Nan(NanError));
            }
            Ok(res)
        };

        // Check for convergence
        if sim[n..]
            .iter()
            .enumerate()
            .all(|(i, &param)| (param - sim[i % n]).abs() < config.x_abs_tol)
            && fsim[1..]
                .iter()
                .all(|Pair { f, .. }| (fsim[0].f - f).abs() < config.f_abs_tol)
        {
            return Ok(true);
        }

        if self.iters == config.max_iters {
            return Err(Error::Convergence(ConvergenceError {
                calls: false,
                val: self.iters,
            }));
        }

        // This can’t panic since `n ≥ 1`.
        let [xbar, x_reflect] = scratch_space.get_disjoint_mut([0..n, n..2 * n]).unwrap();
        for (i, (xbar, x_reflect)) in xbar.iter_mut().zip(&mut *x_reflect).enumerate() {
            *xbar = sim[..n * n].chunks_exact(n).map(|c| c[i]).sum::<f64>() / (n as f64);
            *x_reflect = (1.0 + config.reflection) * *xbar - config.reflection * sim[n * n + i];
        }
        let f_reflect = f(x_reflect)?;

        let new_last = if f_reflect < fsim[0].f {
            // expansion
            for i in 0..n {
                xbar[i] = (1.0 + config.expansion) * xbar[i] - config.expansion * sim[n * n + i];
            }
            let x_expand = xbar;
            let f_expand = f(x_expand)?;
            Some(if f_expand < f_reflect {
                (x_expand, f_expand)
            } else {
                (x_reflect, f_reflect)
            })
        } else if f_reflect < fsim[n - 1].f {
            Some((x_reflect, f_reflect))
        } else if f_reflect < fsim[n].f {
            // outside contraction
            for i in 0..n {
                xbar[i] = (1.0 + config.contraction * config.reflection) * xbar[i]
                    - config.contraction * config.reflection * sim[n * n + i];
            }
            let x_contract = xbar;
            let f_contract = f(x_contract)?;
            (f_contract < f_reflect).then_some((x_contract, f_contract))
        } else {
            // inside contraction
            for i in 0..n {
                xbar[i] =
                    (1.0 - config.contraction) * xbar[i] + config.contraction * sim[n * n + i];
            }
            let x_contract = xbar;
            let f_contract = f(x_contract)?;
            (f_contract < fsim[n].f).then_some((x_contract, f_contract))
        };

        match new_last {
            Some((new_x, new_f)) => {
                sim[n * n..].copy_from_slice(new_x);
                fsim[n].f = new_f;
            }
            None => {
                // do shrink
                for j in 1..=n {
                    for i in 0..n {
                        sim[j * n + i] = sim[i] + config.shrink * (sim[j * n + i] - sim[i]);
                    }
                    fsim[j].f = f(&mut sim[j * n..][..n])?;
                }
            }
        }

        // Sort `fsim` and `sim` simulateously based on the values of `fsim`. There’s no super
        // easy way to do this; we just have to use indices.
        // PANICS: We already ensure that the slice contains no `NaN`s.
        fsim.sort_unstable_by(|&a, &b| a.f.partial_cmp(&b.f).unwrap());
        for (new_i, Pair { index: old_i, .. }) in fsim.iter_mut().enumerate() {
            #[expect(clippy::cast_possible_truncation, reason = "came from a usize")]
            let old_i_usize = *old_i as usize;
            scratch_space[new_i * n..][..n].copy_from_slice(&sim[old_i_usize * n..][..n]);
            *old_i = new_i as u64;
        }
        // Incrementing `self.iterations` also swaps `sim` and `scratch_space`.
        self.iters += 1;

        Ok(false)
    }

    /// Get the number of iterations.
    ///
    /// This is initially zero, and is incremented every time [`Self::step`] returns `Ok`.
    /// This is always `≤ max_iters`.
    #[must_use]
    pub fn iters(&self) -> u64 {
        self.iters
    }

    /// Get the number of times the function has been called.
    /// This is always `≤ max_calls`.
    #[must_use]
    pub fn calls(&self) -> u64 {
        self.calls
    }

    /// Get the best-known `x`-value.
    #[must_use]
    pub fn best_x(&self) -> &V {
        self.simplex().next().unwrap()
    }

    /// Get the full current simplex.
    #[must_use]
    pub fn simplex(&self) -> impl Clone + ExactSizeIterator<Item = &V> {
        let ([sim_a, sim_b], _) = split_buf::<V>(self.config.n, &self.config.buf);
        let sim = match self.iters % 2 != 0 {
            false => sim_b,
            true => sim_a,
        };
        // SAFETY: Initialized by `begin_minimize`.
        unsafe { crate::util::slice_assume_init(sim) }
            .chunks_exact(self.config.n)
            .map(V::from_slice)
    }

    /// Get the best-known `f(x)`-value.
    #[must_use]
    pub fn best_f(&self) -> f64 {
        self.f_simplex().next().unwrap()
    }

    /// Get the `f`-values of the simplex.
    #[must_use]
    pub fn f_simplex(&self) -> impl Clone + ExactSizeIterator<Item = f64> {
        // SAFETY: Initialized by `begin_minimize`.
        let fsim = unsafe {
            crate::util::slice_assume_init(split_buf::<V>(self.config.n, &self.config.buf).1)
        };
        fsim.iter().map(|pair| pair.f)
    }

    /// Regardless of current convergence state, take the best-known `x`-value and get its output.
    #[must_use]
    pub fn into_output(self) -> Output<V> {
        let f = self.best_f();
        Output {
            x: unsafe { V::buf_into_output(self.config.n, self.config.buf) },
            f,
            iters: self.iters,
            calls: self.calls,
        }
    }
}

// `function` not used
#[expect(clippy::missing_fields_in_debug)]
impl<V: ?Sized + Vector + Debug, F: FnMut(&V) -> Result<f64, E>, E> Debug for Minimizer<V, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut s = f.debug_struct("nelder_mead::Minimizer");
        s.field("iters", &self.iters);
        s.field("calls", &self.calls);
        s.field("simplex", &DebugIter(self.simplex()));
        s.field("f_simplex", &DebugIter(self.f_simplex()));
        config_fields(&mut s, &self.config);
        s.finish()
    }
}

struct DebugIter<I: Clone + Iterator<Item: Debug>>(I);
impl<I: Clone + Iterator<Item: Debug>> Debug for DebugIter<I> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.0.clone()).finish()
    }
}

#[derive(Clone, Copy)]
#[repr(C)]
struct Pair {
    index: u64,
    f: f64,
}

fn split_buf<V: ?Sized + Vector>(
    n: usize,
    buf: &V::Buf,
) -> ([&[MaybeUninit<f64>]; 2], &[MaybeUninit<Pair>]) {
    let buf = buf.borrow();
    assert_eq!(buf.len(), n * (n + 1) * 2 + 2 * (n + 1));
    let (sims, fsim) = buf.split_at(n * (n + 1) * 2);
    let (sim_a, sim_b) = sims.split_at(n * (n + 1));
    let fsim = unsafe { slice::from_raw_parts(fsim.as_ptr().cast(), n + 1) };
    ([sim_a, sim_b], fsim)
}

fn split_buf_mut<V: ?Sized + Vector>(
    n: usize,
    buf: &mut V::Buf,
) -> ([&mut [MaybeUninit<f64>]; 2], &mut [MaybeUninit<Pair>]) {
    let buf = buf.borrow_mut();
    assert_eq!(buf.len(), n * (n + 1) * 2 + 2 * (n + 1));
    let (sims, fsim) = buf.split_at_mut(n * (n + 1) * 2);
    let (sim_a, sim_b) = sims.split_at_mut(n * (n + 1));
    let fsim = unsafe { slice::from_raw_parts_mut(fsim.as_mut_ptr().cast(), n + 1) };
    ([sim_a, sim_b], fsim)
}

/// The output of successful [`NelderMead`] minimization.
#[derive(Debug)]
#[non_exhaustive]
pub struct Output<V: ?Sized + Vector> {
    /// The output `x`-vector.
    pub x: V::Owned,
    /// `f(x)`.
    pub f: f64,
    /// The number of iterations required to achieve the result.
    /// This is always `≤ max_iters`.
    pub iters: u64,
    /// The number of function calls required to achieve the result.
    /// This is always `≤ max_calls`.
    pub calls: u64,
}

/// An abstraction over `[f64; N]` and `[f64]`.
/// The latter uses heap allocation.
///
/// # Safety
///
/// `Buf`’s implementation of `Borrow`/`BorrowMut` must return the same data every time.
pub unsafe trait Vector: BorrowMut<[f64]> + ToOwned {
    /// Reconstruct a vector from a slice of floats.
    ///
    /// # Panics
    ///
    /// Is allowed to panic if the slice has the wrong length.
    fn from_slice(vector: &[f64]) -> &Self;

    /// Reconstruct a vector from a mutable slice of floats.
    ///
    /// # Panics
    ///
    /// Is allowed to panic if the slice has the wrong length.
    fn from_slice_mut(vector: &mut [f64]) -> &mut Self;

    /// Return an owned vector representing an unbounded bound. This may either:
    /// - return an empty vector, or
    /// - if `positive`, return a vector of `+∞`, and otherwise a vector of `-∞`.
    fn unbounded(positive: bool) -> Self::Owned;

    /// A buffer with `2N² + 4N + 2` floats.
    type Buf: BorrowMut<[MaybeUninit<f64>]>;

    /// Allocate an uninitialized buffer.
    fn allocate(n: usize) -> Self::Buf;

    /// Construct the output of the computation from the first `n` `f64`s of the buffer.
    ///
    /// # Safety
    ///
    /// The first `n` `f64`s are initialized.
    unsafe fn buf_into_output(n: usize, buf: Self::Buf) -> Self::Owned;
}

unsafe impl<const N: usize> Vector for [f64; N] {
    fn from_slice(vector: &[f64]) -> &Self {
        vector.try_into().unwrap()
    }
    fn from_slice_mut(vector: &mut [f64]) -> &mut Self {
        vector.try_into().unwrap()
    }
    fn unbounded(positive: bool) -> Self::Owned {
        [if positive {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        }; N]
    }
    type Buf = array_helper::Buf<N>;
    fn allocate(_: usize) -> Self::Buf {
        unsafe { <MaybeUninit<Self::Buf>>::uninit().assume_init() }
    }
    unsafe fn buf_into_output(_: usize, mut buf: Self::Buf) -> Self {
        let borrowed: &mut [MaybeUninit<f64>] = buf.borrow_mut();
        unsafe { crate::util::slice_assume_init_mut(&mut borrowed[..N]) }
            .try_into()
            .unwrap()
    }
}

mod array_helper {
    #[repr(C)]
    #[expect(missing_debug_implementations)]
    pub struct Buf<const N: usize> {
        a: [[[MaybeUninit<f64>; N]; N]; 2],
        b: [[MaybeUninit<f64>; N]; 4],
        c: [MaybeUninit<f64>; 2],
    }
    impl<const N: usize> Borrow<[MaybeUninit<f64>]> for Buf<N> {
        fn borrow(&self) -> &[MaybeUninit<f64>] {
            let this: *const Self = self;
            unsafe { slice::from_raw_parts(this.cast(), size_of::<Self>() / size_of::<f64>()) }
        }
    }
    impl<const N: usize> BorrowMut<[MaybeUninit<f64>]> for Buf<N> {
        fn borrow_mut(&mut self) -> &mut [MaybeUninit<f64>] {
            let this: *mut Self = self;
            unsafe { slice::from_raw_parts_mut(this.cast(), size_of::<Self>() / size_of::<f64>()) }
        }
    }

    use core::slice;
    use std::borrow::Borrow;
    use std::borrow::BorrowMut;
    use std::mem::MaybeUninit;
}

unsafe impl Vector for [f64] {
    fn from_slice(vector: &[f64]) -> &Self {
        vector
    }
    fn from_slice_mut(vector: &mut [f64]) -> &mut Self {
        vector
    }
    fn unbounded(_: bool) -> Self::Owned {
        Vec::new()
    }
    type Buf = Box<[MaybeUninit<f64>]>;
    fn allocate(n: usize) -> Self::Buf {
        Box::new_uninit_slice(2 * n * n + 4 * n + 2)
    }
    unsafe fn buf_into_output(n: usize, buf: Self::Buf) -> Vec<f64> {
        let ptr = Box::into_raw(buf);
        unsafe { Vec::from_raw_parts(ptr.cast::<f64>(), n, ptr.len()) }
    }
}

mod hidden {}

/// An error returned by [`NelderMead`].
#[derive(Debug)]
pub enum Error<E = Infallible> {
    /// A tolerance parameter was invalid.
    InvalidTol(InvalidTol),

    /// A bound was invalid.
    InvalidBound(InvalidBound),

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
            Self::InvalidBound(e) => Some(e),
            Self::Convergence(e) => Some(e),
            Self::Nan(e) => Some(e),
            Self::Function(e) => Some(e),
        }
    }
}

use super::ConvergenceError;
use super::InvalidBound;
use super::InvalidTol;
use super::NanError;
use core::slice;
use std::borrow::Borrow;
use std::borrow::BorrowMut;
use std::convert::Infallible;
use std::fmt;
use std::fmt::Debug;
use std::fmt::Display;
use std::fmt::Formatter;
use std::mem::MaybeUninit;
