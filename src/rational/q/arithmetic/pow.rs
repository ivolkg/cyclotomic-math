// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module provides an implementation of the [`Pow`] trait for [`Q`].

use crate::{
    error::MathError,
    integer::Z,
    macros::for_others::{implement_for_others, implement_for_owned},
    rational::Q,
    traits::Pow,
};
use flint_sys::fmpq::fmpq_pow_fmpz;

impl Pow<&Z> for Q {
    type Output = Q;

    /// Raises the value of `self` to the power of an integer `exp`.
    ///
    /// Parameters:
    /// - `exp`: specifies the exponent to which the value is raised
    ///
    /// Returns the value of `self` powered by `exp` as a new [`Q`] instance.
    ///
    /// # Example
    /// ```
    /// use qfall_math::{rational::Q, integer::Z};
    /// use qfall_math::traits::*;
    ///
    /// let base = Q::try_from((&3, &1)).unwrap();
    /// let exp = Z::from(-2);
    ///
    /// let powered_value = base.pow(&exp).unwrap();
    ///
    /// assert_eq!(Q::try_from((&1, &9)).unwrap(), powered_value);
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidExponent`](MathError::InvalidExponent)
    /// if the provided exponent is negative and the base value of `self` is not invertible.
    fn pow(&self, exp: &Z) -> Result<Self::Output, MathError> {
        if self == &Q::ZERO && exp < &Z::ZERO {
            return Err(MathError::InvalidExponent(format!(
                "A negative exponent {} was used for a zero value. There's no inverse for zero values.",
                exp,
            )));
        }

        let mut out = Q::ZERO;
        unsafe { fmpq_pow_fmpz(&mut out.value, &self.value, &exp.value) };
        Ok(out)
    }
}

implement_for_owned!(Z, Q, Pow);
implement_for_others!(Z, Q, Pow for u8 u16 u32 u64 i8 i16 i32 i64);

#[cfg(test)]
mod test_pow {
    use super::*;

    /// Ensure that `pow` works correctly for zero values
    #[test]
    fn zero() {
        let zero = Q::ZERO;

        assert_eq!(Q::ONE, zero.pow(0).unwrap());
        assert_eq!(Q::ZERO, zero.pow(1).unwrap());
        assert!(zero.pow(-1).is_err());
    }

    /// Ensure that `pow` works correctly for base values `1` and `-1`
    #[test]
    fn one() {
        let base_pos = Q::try_from((&1, &1)).unwrap();
        let base_neg = Q::try_from((&-1, &1)).unwrap();

        assert_eq!(Q::ONE, base_pos.pow(0).unwrap());
        assert_eq!(Q::ONE, base_pos.pow(1).unwrap());
        assert_eq!(Q::ONE, base_pos.pow(-2).unwrap());
        assert_eq!(Q::ONE, base_pos.pow(5).unwrap());
        assert_eq!(Q::ONE, base_neg.pow(0).unwrap());
        assert_eq!(Q::MINUS_ONE, base_neg.pow(1).unwrap());
        assert_eq!(Q::ONE, base_neg.pow(-2).unwrap());
        assert_eq!(Q::MINUS_ONE, base_neg.pow(5).unwrap());
    }

    /// Ensure that `pow` works for [`Q`] properly for small values
    #[test]
    fn small() {
        let base_0 = Q::try_from((&2, &1)).unwrap();
        let base_1 = Q::try_from((&1, &2)).unwrap();
        let exp_pos = Z::from(4);

        let res_0 = base_0.pow(&exp_pos).unwrap();
        let res_1 = base_0.pow(0).unwrap();
        let res_2 = base_1.pow(&exp_pos).unwrap();
        let res_3 = base_1.pow(0).unwrap();

        assert_eq!(Q::try_from((&16, &1)).unwrap(), res_0);
        assert_eq!(Q::try_from((&1, &1)).unwrap(), res_1);
        assert_eq!(Q::try_from((&1, &16)).unwrap(), res_2);
        assert_eq!(Q::try_from((&1, &1)).unwrap(), res_3);
    }

    /// Ensure that `pow` works for [`Q`] properly for large values
    #[test]
    fn large() {
        let base_0 = Q::try_from((&i64::MIN, &1)).unwrap();
        let base_1 = Q::try_from((&1, &i64::MIN)).unwrap();
        let exp_pos = Z::from(3);
        let cmp_0 = &base_0 * &base_0 * &base_0;
        let cmp_1 = &base_1 * &base_1 * &base_1;

        let res_0 = base_0.pow(&exp_pos).unwrap();
        let res_1 = base_0.pow(0).unwrap();
        let res_2 = base_0.pow(-1).unwrap();
        let res_3 = base_1.pow(&exp_pos).unwrap();
        let res_4 = base_1.pow(0).unwrap();
        let res_5 = base_1.pow(-1).unwrap();

        assert_eq!(cmp_0, res_0);
        assert_eq!(Q::ONE, res_1);
        assert_eq!(Q::try_from((&1, &i64::MIN)).unwrap(), res_2);
        assert_eq!(cmp_1, res_3);
        assert_eq!(Q::ONE, res_4);
        assert_eq!(Q::try_from((&i64::MIN, &1)).unwrap(), res_5);
    }

    /// Ensures that the `pow` trait is available for other types
    #[test]
    fn availability() {
        let base = Q::try_from((&i64::MAX, &1)).unwrap();
        let exp = Z::from(4);

        let _ = base.pow(exp);
        let _ = base.pow(2_i8);
        let _ = base.pow(2_i16);
        let _ = base.pow(2_i32);
        let _ = base.pow(2_i64);
        let _ = base.pow(2_u8);
        let _ = base.pow(2_u16);
        let _ = base.pow(2_u32);
        let _ = base.pow(2_u64);
    }

    /// Ensures that `pow` returns an error if a non-invertible basis,
    /// i.e. for [`Q`] only 0, is powered by a negative exponent
    #[test]
    fn non_invertible_detection() {
        let base = Q::try_from((&0, &1)).unwrap();

        assert!(base.pow(-1).is_err());
    }
}
