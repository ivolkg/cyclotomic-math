// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module includes core functionality to sample according to the
//! uniform random distribution.

use crate::{error::MathError, integer::Z};
use rand::{rngs::ThreadRng, RngCore};

#[allow(dead_code)]
static mut RNG: Option<ThreadRng> = None;

/// Returns a pointer to a static [`ThreadRng`] to sample uniform at random.
/// If the Random Generator was not used/ instantiated before,
/// it initializes the static variable.
///
/// **WARNING:** This function is strictly for single-threaded use!
/// Otherwise access needs to be synchronized.
///
/// # Examples
/// ```compile_fail
/// use qfall_math::utils::sample::get_rng;
///
/// let rng = get_rng();
/// ```
#[allow(dead_code)]
pub(crate) fn get_rng() -> &'static mut ThreadRng {
    if unsafe { RNG.is_none() } {
        // Instantiates a fresh thread-local RNG initialized with
        // a random seed chosen from an OS-dependent RNG.
        // ThreadRng uses ChaCha12, which is considered cryptographically secure.
        unsafe { RNG = Some(ThreadRng::default()) };
    }

    unsafe { RNG.as_mut().unwrap() }
}

/// Computes a uniform at random chosen [`Z`] sample in `[0, interval_size)`.
///
/// Parameters:
/// - `interval_size`: specifies the size of the interval
/// over which the samples are drawn
///
/// Returns a uniform at random chosen [`Z`] instance in `[0, interval_size)`.
///
/// # Examples
/// ```compile_fail
/// use qfall_math::utils::sample::sample_uniform_rejection;
/// use qfall_math::integer::Z;
/// let interval_size = Z::from(20);
///
/// let sample = sample_uniform_rejection(&interval_size).unwrap();
///
/// assert!(Z::ZERO <= sample);
/// assert!(sample < interval_size);
/// ```
///
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`InvalidInterval`](MathError::InvalidInterval)
/// if the interval is chosen smaller than or equal to `1`.
#[allow(dead_code)]
pub(crate) fn sample_uniform_rejection(interval_size: &Z) -> Result<Z, MathError> {
    if interval_size <= &Z::ONE {
        return Err(MathError::InvalidInterval(format!(
            "An invalid interval size {} was provided.",
            interval_size
        )));
    }

    let bit_size = interval_size.bits() as usize;
    let bit_vector = sample_bits_uniform(bit_size);

    let random = Z::from(&bit_vector);
    if &random < interval_size {
        Ok(random)
    } else {
        sample_uniform_rejection(interval_size)
    }
}

/// Computes `nr_bits` many uniform at random chosen bits.
///
/// Parameters:
/// - `nr_bits`: specifies the number of bits to be chosen uniform at random
///
/// Returns a [`Vec<u8>`] including `nr_bits` many uniform at random chosen bits,
/// filled with `0` bits if needed for the last byte resp. [`u8`].
///
/// # Examples
/// ```compile_fail
/// use qfall_math::utils::sample::sample_bits_uniform;
/// let nr_bits = 14;
///
/// let byte_vector = sample_bits_uniform(nr_bits);
///
/// assert_eq!(byte_vector[1] < 64);
/// assert_eq!(2, byte_vector.len());
/// ```
#[allow(dead_code)]
fn sample_bits_uniform(nr_bits: usize) -> Vec<u8> {
    let rng = get_rng();

    // sample ⌈ nr_bits / 8 ⌉ bytes
    let mut byte_vector;
    if nr_bits % 8 == 0 {
        byte_vector = vec![0u8; nr_bits / 8];
    } else {
        byte_vector = vec![0u8; nr_bits / 8 + 1];
    }
    rng.fill_bytes(&mut byte_vector);

    // set superfluous bits at the end to `0`
    if nr_bits % 8 != 0 {
        let last_index = byte_vector.len() - 1;
        byte_vector[last_index] %= 2u8.pow(nr_bits as u32 % 8);
    }

    byte_vector
}

#[cfg(test)]
mod test_get_rng {
    use rand::RngCore;

    use super::{get_rng, RNG};

    /// Checks whether the first initialization of the static RNG variable works
    /// correctly for usage
    #[test]
    fn init() {
        get_rng();
        let rng = unsafe { RNG.as_mut().unwrap() };
        let _ = rng.next_u32();
    }

    // TODO: Check whether the same instance is reused once `get_rng` is used a second time.
    // problem: Eq is not implemented by `ThreadRng`
}

#[cfg(test)]
mod test_sample_uniform_rejection {
    use super::{sample_uniform_rejection, Z};

    /// Ensures that the doc tests works correctly.
    #[test]
    fn doc_test() {
        let interval_size = Z::from(20);

        let sample = sample_uniform_rejection(&interval_size).unwrap();

        assert!(Z::ZERO <= sample);
        assert!(sample < interval_size);
    }

    /// Checks whether sampling works fine for small interval sizes
    #[test]
    fn small_interval() {
        let size_2 = Z::from(2);
        let size_7 = Z::from(7);

        let sample_0 = sample_uniform_rejection(&size_2).unwrap();
        let sample_1 = sample_uniform_rejection(&size_2).unwrap();
        let sample_2 = sample_uniform_rejection(&size_2).unwrap();
        let sample_3 = sample_uniform_rejection(&size_7).unwrap();
        let sample_4 = sample_uniform_rejection(&size_7).unwrap();
        let sample_5 = sample_uniform_rejection(&size_7).unwrap();

        assert!(sample_0 >= Z::ZERO);
        assert!(sample_0 < size_2);
        assert!(sample_1 >= Z::ZERO);
        assert!(sample_1 < size_2);
        assert!(sample_2 >= Z::ZERO);
        assert!(sample_2 < size_2);
        assert!(sample_3 >= Z::ZERO);
        assert!(sample_3 < size_7);
        assert!(sample_4 >= Z::ZERO);
        assert!(sample_4 < size_7);
        assert!(sample_5 >= Z::ZERO);
        assert!(sample_5 < size_7);
    }

    /// Checks whether sampling works fine for large interval sizes
    #[test]
    fn large_interval() {
        let size = Z::from(i64::MAX);

        for _i in 0..u8::MAX {
            let sample = sample_uniform_rejection(&size).unwrap();

            assert!(sample >= Z::ZERO);
            assert!(sample < size);
        }
    }

    /// Checks whether interval sizes smaller than 2 result in an error
    #[test]
    fn invalid_interval() {
        assert!(sample_uniform_rejection(&Z::ONE).is_err());
        assert!(sample_uniform_rejection(&Z::ZERO).is_err());
        assert!(sample_uniform_rejection(&Z::MINUS_ONE).is_err());
    }
}

#[cfg(test)]
mod test_sample_bits_uniform {
    use super::sample_bits_uniform;

    /// Ensures that the doc tests works correctly.
    #[test]
    fn doc_test() {
        let nr_bits = 14;

        let byte_vector = sample_bits_uniform(nr_bits);

        assert!(byte_vector[1] < 64);
        assert_eq!(2, byte_vector.len());
    }

    /// Checks whether random bit sampling works appropriate for full byte orders
    #[test]
    fn full_bytes() {
        let bits_0 = sample_bits_uniform(8);
        let bits_1 = sample_bits_uniform(16);
        let bits_2 = sample_bits_uniform(256);

        assert_eq!(1, bits_0.len());
        assert_eq!(2, bits_1.len());
        assert_eq!(32, bits_2.len());
    }

    /// Checks whether random bit sampling works appropriate for partial bytes
    /// and that superfluous bits are cutoff appropriately
    #[test]
    fn partial_bytes() {
        let bits_0 = sample_bits_uniform(7);
        let bits_1 = sample_bits_uniform(14);
        let bits_2 = sample_bits_uniform(250);

        assert_eq!(1, bits_0.len());
        assert!(128 > bits_0[0]);
        assert_eq!(2, bits_1.len());
        assert!(64 > bits_1[1]);
        assert_eq!(32, bits_2.len());
        assert!(4 > bits_2[31])
    }

    /// Check whether all necessary bits are chosen uniform at random.
    /// This test could possibly fail with probability 2^-64
    #[test]
    fn random_bits() {
        let mut samples: Vec<u8> = vec![];
        for _i in 0..64 {
            samples.push(sample_bits_uniform(7)[0]);
        }

        for position in 0..7 {
            let mut found_1 = false;
            // check for all 64 samples whether the `position`-th bit is `1`
            for sample in &samples {
                if (sample >> position & 1) % 2 == 1 {
                    found_1 = true;
                    break;
                }
            }

            if !found_1 {
                panic!(
                    "None of the inspected 64 random bits at position {} was 1. 
                    This seems suspicious.",
                    position
                );
            }
        }
    }
}