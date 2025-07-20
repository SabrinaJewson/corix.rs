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
[Scipy](https://scipy.org/),
[CEPHES](https://www.netlib.org/cephes/) and
[cdflib](https://people.sc.fsu.edu/~jburkardt/f_src/cdflib/cdflib.html).
We don’t use R because it’s GPL.

## Implemented algorithms

- Optimization
	- [`Brent`](opt::Brent): Brent’s method (univariate optimization)
	- [`NelderMead`](opt::NelderMead): Nelder–Mead (multivariate optimization)
- Root-finding
	- [`Brent`](root::Brent): Brent (univariate root-finding)

## CEPHES porting notes

```sh
find "$PWD" -name '*.c' -print0 | jq -Rs '[split("\u0000")[] | if . == "" then empty end | {directory:"'"$PWD"'",file:.,arguments:["/usr/bin/clang","-std=c89",.]}]' > compile_commands.json
c2rust transpile compile_commands.json
sed -i 's/std::ffi::c_double/f64/g' *
sed -i 's/std::ffi::c_int/usize/g' *
sed -i 's/unsafe extern "C" //g' *
sed -i 's/#\[no_mangle\]//g' *
```
