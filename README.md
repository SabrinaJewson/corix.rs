Accurate numerics for Rust.
Because like Corixidae, we float well.

In the author’s experience,
existing numerics libraries in Rust
tend to be buggy and low-accuracy
due to a propensity for people to write their own code
instead of copying existing algorithms.
This library takes the opposite approach,
exclusively porting or binding to existing code –
sources include
[CEPHES](https://www.netlib.org/cephes/) and
[Scipy](https://scipy.org/).
We don’t use R because it’s GPL.

## Implemented algorithms

- Optimization
	- [`Brent`](opt::Brent): Brent’s method (univariate optimization)
	- [`NelderMead`](opt::NelderMead): Nelder–Mead (multivariate optimization)
- Root-finding
	- [`Brent`](root::Brent): Brent (univariate root-finding)
- Special functions
	- [`sum`](sum):
		Float summation that is an order of magnitude faster and more accurate than naïve methods.
	- [`mean`](mean) / [`mean_unbiased_sd`](mean_unbiased_sd):
		Mean and unbiased and biased sample standard deviation and variance,
		based on the above summation algorithms.
	- [`erf`](erf), [`gamma`](gamma): Error function, gamma function
		(at the time of writing, the methods on `f64` itself are unstable).

## Stuff that would be nice to have

- Bindings to [CEPHES](https://www.netlib.org/cephes/).
- Bindings to (or a port of) [cdflib](https://people.sc.fsu.edu/~jburkardt/f_src/cdflib/cdflib.html).
- Bindings to (or a port of) [odrpack](https://www.netlib.org/odrpack/).
- Various things from [xsf](https://github.com/scipy/xsf).
- More root-finding and optimization methods.
- Linear programming.
- Forward and reverse automatic differentiation.
- Once we have CEPHES and cdflib, various distributions.
- Numeric derivatives. An `enum Accuracy` parameter to switch between two implementations:
	- low-accuracy based on Scipy’s [`approx_derivative`](https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_numdiff.py);
	- high-accuracy based on Scipy’s [`derivative`](https://github.com/scipy/scipy/blob/v1.16.0/scipy/differentiate/_differentiate.py#L60-L589).
