# `corix.rs`

Accurate numerics for Rust.
Because like Corixidae, we float well.

In the author’s experience,
existing numerics libraries in Rust
tend to be buggy and low-accuracy
due to a propensity for people writing their own code
instead of copying existing algorithms.
This library takes the opposite approach,
exclusively porting or binding to existing code –
sources include
[Scipy](https://scipy.org/),
[CEPHES](https://www.netlib.org/cephes/) and
[cdflib](https://people.sc.fsu.edu/~jburkardt/f_src/cdflib/cdflib.html).
We don’t use R because it’s GPL.
