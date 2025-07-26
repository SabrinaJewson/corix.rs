Accurate numerics for Rust.
Because like Corixidae, we float well.

In the author’s experience,
existing numerics libraries in Rust
tend to be buggy and low-accuracy
due to a propensity for people to write their own code
instead of copying existing algorithms.
This library takes the opposite approach,
mostly porting or binding to existing code
unless we can show we have better accuracy.

Sources include:
- Many algorithms copied from [Scipy](https://scipy.org/).
- Bindings to [xsf](https://github.com/scipy/xsf), Scipy’s C++ library of special functions,
	which itself is a “best of” of several other libraries, including
	[CEPHES](https://www.netlib.org/cephes/),
	[amos](https://www.netlib.org/amos/),
	[specfun](https://www.netlib.org/specfun/) and
	[Boost](https://www.boost.org/).

We don’t use R because it’s GPL.

## Implemented algorithms

- Optimization
	- [`Brent`](opt::Brent): Brent’s method (univariate optimization)
	- [`NelderMead`](opt::NelderMead): Nelder–Mead (multivariate optimization)
- Root-finding
	- [`Brent`](root::Brent): Brent (univariate root-finding)
- Statistics
	- [`sum`]:
		Float summation that is an order of magnitude faster and more accurate than naïve methods.
	- [`mean`] / [`mean_unbiased_sd`]:
		Mean and unbiased and biased sample standard deviation, variance and covariance,
		based on the above summation algorithms.
- Special functions
	- [`erf`], [`gamma`]: Error function, gamma function
		(at the time of writing, the methods on `f64` itself are unstable).
	- [`airy`](mod@airy): The Airy functions Ai and Bi,
		their first derivatives, zeros and integrals.
	- [`bessel`]: Bessel functions.
	- [`ellip::jacobi`]: Jacobi elliptic functions sn, cn, dn, am.
	- [`ellip`]: Complete and incomplete elliptic integrals of the first and second kind.

## Some other crates

A list of crates that occupy vaguely the same space as this one.
Limited to crates that aren’t obviously terrible
(i.e. have docs, use standard Rust style, don’t have excessive dependencies).
For the algorithms we both implement, this crate will likely have better accuracy.

- Linear algebra
	- [`nalgebra`](https://docs.rs/nalgebra):
		Matrices and vectors. They can have statically known widths/heights.
	- [`ndarray`](https://docs.rs/ndarray):
		Full N-dimensional arrays. They only support dynamic dimensions.
	- [`faer`](https://docs.rs/faer):
		Matrices and vectors, but no static dimensions.
- Linear regression
	- [`linregress`](https://docs.rs/linregress): Depends on `nalgebra`.
	- [`linreg`](https://docs.rs/linreg): Dead simple. No big dependencies.
	- [`ndarray-glm`](https://docs.rs/ndarray-glm): Generalized linear models.
- General (including crates with overlap)
	- [`compute`](https://docs.rs/compute): Has lots of overlap with this library.
	- [`std-dev`](https://docs.rs/std-dev): Lots of kinds of regression.
	- [`scirs2`](https://docs.rs/scirs2):
		Meets my criteria, but allegedly wrote 1.5 million lines of code in 74 days,
		(average of 24 000 lines per day), so probably LLM-generated.

## Stuff that would be nice to have

- Bindings to (or a port of) [odrpack](https://www.netlib.org/odrpack/).
- More root-finding and optimization methods.
- Linear programming.
- Forward and reverse automatic differentiation.
- Various distributions (based on xsf functions).
- Numeric derivatives. An `enum Accuracy` parameter to switch between two implementations:
	- low-accuracy based on Scipy’s [`approx_derivative`](https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_numdiff.py);
	- high-accuracy based on Scipy’s [`derivative`](https://github.com/scipy/scipy/blob/v1.16.0/scipy/differentiate/_differentiate.py#L60-L589).
