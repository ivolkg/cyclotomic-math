// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module provides an implementation of the [`Pow`] trait for [`Zq`].

use crate::{
    error::MathError,
    integer::Z,
    integer_mod_q::Zq,
    macros::for_others::{implement_for_others, implement_for_owned},
    traits::Pow,
};
use flint_sys::fmpz_mod::fmpz_mod_pow_fmpz;

impl Pow<&Z> for Zq {
    type Output = Zq;

    /// Raises the value of `self` to the power of an integer `exp`.
    ///
    /// Parameters:
    /// - `exp`: specifies the exponent to which the value is raised
    ///
    /// Returns the value of `self` powered by `exp` as a new `Output` instance.
    ///
    /// # Example
    /// ```
    /// use qfall_math::integer::Z;
    /// use qfall_math::integer_mod_q::Zq;
    /// use qfall_math::traits::*;
    ///
    /// let base = Zq::try_from((2, 9)).unwrap();
    /// let exp = Z::from(4);
    ///
    /// let powered_value = base.pow(&exp).unwrap();
    ///
    /// let cmp = Zq::try_from((7, 9)).unwrap();
    /// assert_eq!(cmp, powered_value);
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidExponent`](MathError::InvalidExponent)
    /// if the provided exponent is negative and the base value of `self` is not invertible.
    fn pow(&self, exp: &Z) -> Result<Self::Output, MathError> {
        let mut out = self.clone();
        if 0 == unsafe {
            fmpz_mod_pow_fmpz(
                &mut out.value.value,
                &self.value.value,
                &exp.value,
                self.modulus.get_fmpz_mod_ctx_struct(),
            )
        } {
            return Err(MathError::InvalidExponent(format!(
                "The negative exponent {} was used for a non-invertible base {}",
                exp, self
            )));
        }
        Ok(out)
    }
}

implement_for_owned!(Z, Zq, Pow);
implement_for_others!(Z, Zq, Pow for u8 u16 u32 u64 i8 i16 i32 i64);

#[cfg(test)]
mod test_pow {
    use super::*;

    /// Ensure that `pow` works for [`Zq`] properly for small values
    #[test]
    fn small() {
        let base = Zq::try_from((2, 9)).unwrap();
        let exp_0 = Z::from(4);
        let exp_1 = Z::ZERO;
        let exp_2 = Z::from(-2);

        let res_0 = base.pow(&exp_0).unwrap();
        let res_1 = base.pow(&exp_1).unwrap();
        let res_2 = base.pow(&exp_2).unwrap();

        assert_eq!(Zq::try_from((7, 9)).unwrap(), res_0);
        assert_eq!(Zq::try_from((1, 9)).unwrap(), res_1);
        assert_eq!(Zq::try_from((7, 9)).unwrap(), res_2);
    }

    /// Ensure that `pow` works for [`Zq`] properly for large values
    #[test]
    fn large() {
        let base = Zq::try_from((i64::MAX, u64::MAX)).unwrap();
        let exp_0 = Z::from(4);
        let exp_1 = Z::ZERO;
        let exp_2 = Z::from(-1);
        let cmp_0 = Zq::try_from((1152921504606846976_i64, u64::MAX)).unwrap();
        let cmp_1 = Zq::try_from((1, u64::MAX)).unwrap();
        let cmp_2 = Zq::try_from((18446744073709551613_u64, u64::MAX)).unwrap();

        let res_0 = base.pow(&exp_0).unwrap();
        let res_1 = base.pow(&exp_1).unwrap();
        let res_2 = base.pow(&exp_2).unwrap();

        assert_eq!(cmp_0, res_0);
        assert_eq!(cmp_1, res_1);
        assert_eq!(cmp_2, res_2);
    }

    /// Ensures that the `pow` trait is available for other types
    #[test]
    fn availability() {
        let base = Zq::try_from((i64::MAX, u64::MAX)).unwrap();
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

    /// Ensures that `pow` returns an error if a non-invertible basis is
    /// powered by a negative exponent
    #[test]
    fn non_invertible_detection() {
        let base = Zq::try_from((2, 4)).unwrap();

        assert!(base.pow(-1).is_err());
    }
}
