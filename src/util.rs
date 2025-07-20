pub(crate) fn write_copy_of_slice<'dest, T: Copy>(
    dest: &'dest mut [MaybeUninit<T>],
    src: &[T],
) -> &'dest mut [T] {
    let src: &[MaybeUninit<T>] = unsafe { slice::from_raw_parts(src.as_ptr().cast(), src.len()) };
    dest.copy_from_slice(src);
    unsafe { slice_assume_init_mut(dest) }
}

pub(crate) unsafe fn slice_assume_init<T>(slice: &[MaybeUninit<T>]) -> &[T] {
    unsafe { slice::from_raw_parts(slice.as_ptr().cast(), slice.len()) }
}

pub(crate) unsafe fn slice_assume_init_mut<T>(slice: &mut [MaybeUninit<T>]) -> &mut [T] {
    unsafe { slice::from_raw_parts_mut(slice.as_mut_ptr().cast(), slice.len()) }
}

pub(crate) fn slice_try_init_with<T: Copy, E>(
    slice: &mut [MaybeUninit<T>],
    mut f: impl FnMut(usize) -> Result<T, E>,
) -> Result<&mut [T], E> {
    for (i, item) in slice.iter_mut().enumerate() {
        *item = MaybeUninit::new(f(i)?);
    }
    Ok(unsafe { slice_assume_init_mut(slice) })
}

macro_rules! setter {
    ($setter:ident => $field:ident: $fty:ty) => {
        #[doc = concat!("Set [`Self::", stringify!($field), "`].")]
        #[must_use]
        pub fn $setter(mut self, $field: $fty) -> Self {
            self.$field = $field;
            self
        }
    };
}
pub(crate) use setter;

use core::slice;
use std::mem::MaybeUninit;
