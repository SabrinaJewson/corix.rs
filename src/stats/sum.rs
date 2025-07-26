/// Sum a slice of floats at high accuracy and speed.
///
/// On my Zen 4 machine, compared to a naïve `.iter().sum()`, this is:
/// - 5× faster and 2× as accurate for small (≈50) datasets, and
/// - 10× faster and 60x more accurate for large (≈100 000) datasets.
///
/// See also [`sum_with`], a generalized version of this function.
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::sum(&[5.0, 2.0, -6.5]), 0.5_f64);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn sum<F>(arr: &[F]) -> F
where
    F: Default + Clone + for<'f> AddAssign<&'f F> + for<'f> SubAssign<&'f F>,
{
    sum_with(arr, |acc, val| *acc += val)
}

/// Calculate the mean of a list of values.
///
/// Returns [`f64::NAN`] if the list is empty.
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::mean(&[5.0, 2.0, -1.3]), 1.9);
/// assert!(corix::mean(&[]).is_nan());
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn mean(arr: &[f64]) -> f64 {
    sum(arr) / (arr.len() as f64)
}

/// Calculate the mean and unbiased sample variance of a list of values.
///
/// The denominator is `n - 1` – see [`mean_biased_var`] for the `n` version.
///
/// See also [`mean_unbiased_sd`] for the standard deviation.
///
/// Returns [`f64::NAN`] if the list is empty.
///
/// # Examples
///
/// ```
/// let (mean, var) = corix::mean_unbiased_var(&[5.0, 2.0, -1.3]);
/// assert_float_relative_eq!(mean, 1.9);
/// assert_float_relative_eq!(var, 9.93);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn mean_unbiased_var(arr: &[f64]) -> (f64, f64) {
    let mean = mean(arr);
    let var =
        sum_with(arr, |acc: &mut f64, val| *acc += (val - mean).powi(2)) / (arr.len() as f64 - 1.0);
    (mean, var)
}

/// Calculate the mean and biased sample variance of a list of values.
///
/// The denominator is `n` – see [`mean_unbiased_var`] for the `n - 1` version.
///
/// See also [`mean_biased_sd`] for the standard deviation.
///
/// Returns [`f64::NAN`] if the list is empty.
///
/// # Examples
///
/// ```
/// let (mean, var) = corix::mean_biased_var(&[5.0, 2.0, -1.3]);
/// assert_float_relative_eq!(mean, 1.9);
/// assert_float_relative_eq!(var, 6.62);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn mean_biased_var(arr: &[f64]) -> (f64, f64) {
    let mean = mean(arr);
    let var = sum_with(arr, |acc: &mut f64, val| *acc += (val - mean).powi(2)) / (arr.len() as f64);
    (mean, var)
}

/// Calculate the mean and unbiased sample standard deviation of a list of values.
///
/// The denominator is `n - 1` – see [`mean_biased_sd`] for the `n` version.
///
/// See also [`mean_unbiased_var`] for the variance.
///
/// Returns [`f64::NAN`] if the list is empty.
///
/// # Examples
///
/// ```
/// let (mean, sd) = corix::mean_unbiased_sd(&[5.0, 2.0, -1.3]);
/// assert_float_relative_eq!(mean, 1.9);
/// assert_float_relative_eq!(sd, 3.151_190_251_317_746);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn mean_unbiased_sd(arr: &[f64]) -> (f64, f64) {
    let (mean, var) = mean_unbiased_var(arr);
    (mean, var.sqrt())
}

/// Calculate the mean and biased sample standard deviation of a list of values.
///
/// The denominator is `n` – see [`mean_unbiased_sd`] for the `n - 1` version.
///
/// See also [`mean_biased_var`] for the variance.
///
/// Returns [`f64::NAN`] if the list is empty.
///
/// # Examples
///
/// ```
/// let (mean, sd) = corix::mean_biased_sd(&[5.0, 2.0, -1.3]);
/// assert_float_relative_eq!(mean, 1.9);
/// assert_float_relative_eq!(sd, 2.572_936_066_053_721);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn mean_biased_sd(arr: &[f64]) -> (f64, f64) {
    let (mean, var) = mean_biased_var(arr);
    (mean, var.sqrt())
}

/// Calculate the mean and unbiased sample covariance of two lists.
///
/// Returns the mean of `a`, the mean of `b`, and the covariance.
///
/// The denominator is `n - 1` – see [`mean_biased_cov`] for the `n` version.
///
/// Returns [`f64::NAN`] if both lists are empty.
///
/// # Panics
///
/// Panics if `a` and `b` are different lengths.
///
/// # Examples
///
/// ```
/// let (mean_a, mean_b, cov) = corix::mean_unbiased_cov(&[5.0, 2.0, -1.3], &[5.2, 1.8, -1.6]);
/// assert_float_relative_eq!(mean_a, 1.9);
/// assert_float_relative_eq!(mean_b, 1.8);
/// assert_float_relative_eq!(cov, 10.71);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn mean_unbiased_cov(a: &[f64], b: &[f64]) -> (f64, f64, f64) {
    assert_eq!(a.len(), b.len());
    let [a_mean, b_mean] = [mean(a), mean(b)];
    let mut b = b.iter().copied();
    let cov = sum_with(a, |acc: &mut f64, val| {
        *acc += (val - a_mean) * (b.next().unwrap() - b_mean);
    }) / (a.len() as f64 - 1.0);
    (a_mean, b_mean, cov)
}

/// Calculate the mean and biased sample covariance of two lists.
///
/// Returns the mean of `a`, the mean of `b`, and the covariance.
///
/// The denominator is `n` – see [`mean_unbiased_cov`] for the `n - 1` version.
///
/// Returns [`f64::NAN`] if both lists are empty.
///
/// # Panics
///
/// Panics if `a` and `b` are different lengths.
///
/// # Examples
///
/// ```
/// let (mean_a, mean_b, cov) = corix::mean_biased_cov(&[5.0, 2.0, -1.3], &[5.2, 1.8, -1.6]);
/// assert_float_relative_eq!(mean_a, 1.9);
/// assert_float_relative_eq!(mean_b, 1.8);
/// assert_float_relative_eq!(cov, 7.14);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn mean_biased_cov(a: &[f64], b: &[f64]) -> (f64, f64, f64) {
    assert_eq!(a.len(), b.len());
    let [a_mean, b_mean] = [mean(a), mean(b)];
    let mut b = b.iter().copied();
    let cov = sum_with(a, |acc: &mut f64, val| {
        *acc += (val - a_mean) * (b.next().unwrap() - b_mean);
    }) / (a.len() as f64);
    (a_mean, b_mean, cov)
}

/// Sum a slice of data into a single float at high accuracy and speed.
/// Like [`sum`], but allows passing in a custom function.
///
/// The function is expected to add a value to its first argument.
/// It will be called for every element of the slice, in order.
///
/// Works best if the function is vectorizable.
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(
///     corix::sum_with(&[2.0_f64, 3.0, -4.0], |acc: &mut f64, v| *acc += v.powi(2)),
///     29.0
/// );
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn sum_with<T, F>(arr: &[T], mut f: impl FnMut(&mut F, &T)) -> F
where
    F: Default + Clone + for<'f> AddAssign<&'f F> + for<'f> SubAssign<&'f F>,
{
    // Using the `sum_orlp` algorithm from <https://orlp.net/blog/taming-float-sums/>, but adapted
    // to work on stable. Thanks orlp!

    // To compensate for the reduced accuracy in `sum_block` for large `F`, we perform error
    // correction more frequently.
    let mut chunks = arr.chunks_exact(match size_of::<F>() {
        ..=256 => 256,
        257.. => 192,
    });
    let mut sum = F::default();
    let mut c = F::default();
    let mut t = F::default();
    for chunk in &mut chunks {
        // Add the current value onto `c`, giving us `c + current`.
        sum_block(chunk, &mut c, &mut f);

        // Remember the current sum in `t`.
        t.clone_from(&sum);

        // Add the current value and `c` onto the full sum, giving `sum + current + c + ERROR`.
        sum += &c;

        // Subtract the new sum from the old sum, giving us `-(current + c + ERROR)`
        t -= &sum;

        // Add this `t` term onto `c + current`. The `current` and `c` cancel out, giving us just
        // the `ERROR`.
        c += &t;
    }
    sum_block(chunks.remainder(), &mut c, &mut f);
    sum += &c;
    sum
}

#[test]
#[expect(clippy::float_cmp)]
fn sum_many() {
    assert_eq!(sum(&vec![1.0; 1_000]), 1_000.0);
    assert_eq!(sum(&vec![1.0; 1_000_000]), 1_000_000.0);
}

fn sum_block<T, F>(arr: &[T], total: &mut F, f: &mut impl FnMut(&mut F, &T))
where
    F: Default + for<'f> AddAssign<&'f F>,
{
    // We try to avoid overflowing the stack, at the cost of some accuracy, in case `F` is very
    // large.
    let acc: [F; 0] = [];
    if size_of::<F>() <= 256 {
        let (acc, arr) = sum_step::<T, F, 32>(&acc, arr, f);
        let (acc, arr) = sum_step::<T, F, 16>(&acc, arr, f);
        let (acc, arr) = sum_step::<T, F, 8>(&acc, arr, f);
        acc.iter().for_each(|val| *total += val);
        arr.iter().for_each(|val| f(total, val));
    } else if size_of::<F>() <= 512 {
        let (acc, arr) = sum_step::<T, F, 16>(&acc, arr, f);
        let (acc, arr) = sum_step::<T, F, 8>(&acc, arr, f);
        acc.iter().for_each(|val| *total += val);
        arr.iter().for_each(|val| f(total, val));
    } else {
        let (acc, arr) = sum_step::<T, F, 8>(&acc, arr, f);
        acc.iter().for_each(|val| *total += val);
        arr.iter().for_each(|val| f(total, val));
    }
}

fn sum_step<'a, T, F: Default + for<'f> AddAssign<&'f F>, const N: usize>(
    acc: &[F],
    arr: &'a [T],
    f: &mut impl FnMut(&mut F, &T),
) -> ([F; N], &'a [T]) {
    let mut total = array::from_fn(|_| F::default());
    assert_eq!(acc.len() % N, 0);
    for chunk in acc.chunks_exact(N) {
        for (total, val) in total.iter_mut().zip(chunk) {
            *total += val;
        }
    }
    let mut chunks = arr.chunks_exact(N);
    for chunk in &mut chunks {
        for (total, val) in total.iter_mut().zip(chunk) {
            f(total, val);
        }
    }
    (total, chunks.remainder())
}

use std::array;
use std::ops::AddAssign;
use std::ops::SubAssign;
