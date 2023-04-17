// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module includes functionality about properties of [`Z`] instances.

use super::Z;
use crate::rational::Q;
use flint_sys::{
    fmpq::{fmpq, fmpq_inv},
    fmpz::{fmpz, fmpz_abs, fmpz_is_prime},
};

impl Z {
    /// Checks if a [`Z`] is prime.
    ///
    /// Returns true if the value is prime.
    ///
    /// ```
    /// use qfall_math::integer::Z;
    ///
    /// let value = Z::from(17);
    /// assert!(value.is_prime())
    /// ```
    pub fn is_prime(&self) -> bool {
        1 == unsafe { fmpz_is_prime(&self.value) }
    }

    /// Returns the given [`Z`] instance with its absolute value.
    ///
    /// # Example
    /// ```
    /// use qfall_math::integer::Z;
    /// let mut value = Z::from(-1);
    ///
    /// let value = value.abs();
    ///
    /// assert_eq!(Z::ONE, value);
    /// ```
    pub fn abs(mut self) -> Self {
        unsafe {
            fmpz_abs(&mut self.value, &self.value);
        }
        self
    }

    /// Returns the inverse of `self` as a fresh [`Q`] instance.
    ///
    /// As the inverse of `0` is undefined, it returns `None` in case `self == 0`.
    ///
    /// # Example
    /// ```
    /// use qfall_math::{integer::Z, rational::Q};
    /// let value = Z::from(4);
    ///
    /// let inverse = value.inv().unwrap();
    ///
    /// assert_eq!(Q::try_from((&1, &4)).unwrap(), inverse);
    /// ```
    pub fn inv(&self) -> Option<Q> {
        if self == &Z::ZERO {
            return None;
        }

        let mut out = Q::ZERO;
        // the manual construction of fmpq removes the need to clone self's value/
        // the numerator. the new fmpz value does not need to be cleared manually
        // as it's small the fmpq instance does neither as the fmpq value is
        // dropped automatically, but the numerator/ self's value is kept alive
        let self_fmpq = fmpq {
            num: self.value,
            den: fmpz(1),
        };
        unsafe { fmpq_inv(&mut out.value, &self_fmpq) };
        Some(out)
    }
}

#[cfg(test)]
mod test_abs {
    use super::Z;

    /// Checks whether `abs` returns the positive value for small values correctly
    #[test]
    fn small_values() {
        let pos = Z::ONE;
        let zero = Z::ZERO;
        let neg = Z::from(-15);

        assert_eq!(Z::ONE, pos.abs());
        assert_eq!(Z::ZERO, zero.abs());
        assert_eq!(Z::from(15), neg.abs());
    }

    /// Checks whether `abs` returns the positive value for large values correctly
    #[test]
    fn large_values() {
        let pos = Z::from(i64::MAX);
        let neg = Z::from(i64::MIN);

        assert_eq!(Z::from(i64::MAX), pos.abs());
        assert_eq!(Z::from(i64::MIN) * Z::from(-1), neg.abs());
    }
}

#[cfg(test)]
mod test_is_prime {
    use super::Z;

    /// ensure that primes are correctly detected
    #[test]
    fn prime_detection() {
        let small = Z::from(2_i32.pow(16) + 1);
        let large = Z::from(u64::MAX - 58);
        assert!(small.is_prime());
        assert!(large.is_prime());
    }

    /// ensure that non-primes are correctly detected
    #[test]
    fn non_prime_detection() {
        let small = Z::from(2_i32.pow(16));
        let large = Z::from(i64::MAX);
        assert!(!small.is_prime());
        assert!(!large.is_prime());
    }
}

#[cfg(test)]
mod test_inv {
    use super::{Q, Z};

    /// Checks whether the inverse is correctly computed for small values
    #[test]
    fn small_values() {
        let val_0 = Z::from(4);
        let val_1 = Z::from(-7);

        let inv_0 = val_0.inv().unwrap();
        let inv_1 = val_1.inv().unwrap();

        assert_eq!(Q::try_from((&1, &4)).unwrap(), inv_0);
        assert_eq!(Q::try_from((&-1, &7)).unwrap(), inv_1);
    }

    /// Checks whether the inverse is correctly computed for large values
    #[test]
    fn large_values() {
        let val_0 = Z::from(i64::MAX);
        let val_1 = Z::from(i64::MIN);

        let inv_0 = val_0.inv().unwrap();
        let inv_1 = val_1.inv().unwrap();

        assert_eq!(Q::try_from((&1, &i64::MAX)).unwrap(), inv_0);
        assert_eq!(Q::try_from((&1, &i64::MIN)).unwrap(), inv_1);
    }

    /// Checks whether the inverse of `0` returns `None`
    #[test]
    fn inv_zero_none() {
        let zero = Z::ZERO;

        let inv_zero = zero.inv();

        assert!(inv_zero.is_none());
    }
}
