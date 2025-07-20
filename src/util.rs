/// A 2-tuple, but `repr(C)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(C)]
pub struct ReprCTuple<T, U> {
    /// The first element of the tuple.
    pub a: T,
    /// The second element of the tuple.
    pub b: U,
}

impl<T, const N: usize> ReprCTuple<[T; N], T> {
    pub(crate) fn as_first_array_flattened_mut(&mut self) -> &mut [T] {
        let ptr: *mut Self = self;
        // SAFETY: We are `repr(C)`, and arrays are guaranteed to be layed out as slices, so
        // overall we are layed out as a slice.
        unsafe { slice::from_raw_parts_mut(ptr.cast(), N + 1) }
    }
}

pub(crate) fn array_try_from_fn<T, E, const N: usize>(
    mut f: impl FnMut(usize) -> Result<T, E>,
) -> Result<[T; N], E> {
    let mut array = <MaybeUninit<[T; N]>>::uninit();
    // TODO: stop memory leaks
    for i in 0..N {
        unsafe { array.as_mut_ptr().cast::<T>().add(i).write(f(i)?) };
    }
    Ok(unsafe { array.assume_init() })
}

use core::slice;
use std::mem::MaybeUninit;
