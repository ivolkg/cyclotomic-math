// Copyright © 2023 Marvin Beckmann, Sven Moog
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implementations to create a [`Modulus`] value from other types.
//! For each reasonable type, an explicit function with the format
//! `from_<type_name>` and the [`From`] trait should be implemented.
//!
//! The explicit functions contain the documentation.

use super::Modulus;
use crate::{error::MathError, integer::Z};
use flint_sys::fmpz_mod::{fmpz_mod_ctx, fmpz_mod_ctx_init};
use std::{mem::MaybeUninit, rc::Rc, str::FromStr};

impl Modulus {
    /// Create a [`Modulus`] from [`Z`].
    ///
    /// Parameters:
    /// - `value`: the value of the modulus
    ///
    /// Returns a [`Modulus`] or a [`MathError`]
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::Modulus;
    /// use std::str::FromStr;
    ///
    /// let modulus = Modulus::from_str("42").unwrap();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntToModulus`](MathError::InvalidIntToModulus)
    /// if the provided value is not greater than `1`.
    pub fn try_from_z(value: &Z) -> Result<Self, MathError> {
        Modulus::from_fmpz_ref(&value.value)
    }

    /// Initialize a [`Modulus`] from an fmpz reference.
    ///
    /// Parameters:
    /// - `value`: the initial value the modulus should have.
    ///   It must be larger than one.
    ///
    /// Returns the new [`Modulus`] or an error, if the
    /// provided value was not greater than `1`.
    ///
    /// # Safety
    /// Since the parameter is a reference, it still has to be
    /// properly cleared outside this function.
    /// For example, by the drop trait of [`Z`].
    ///
    /// # Examples
    /// ```compile_fail
    /// use qfall_math::integer_mod_q::Modulus;
    /// use qfall_math::integer::Z;
    ///
    /// let z = Z::from(42);
    /// let modulus = Modulus::from_fmpz_ref(&z.value);
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntToModulus`](MathError::InvalidIntToModulus)
    /// if the provided value is not greater than `1`.
    pub(crate) fn from_fmpz_ref(value: &fmpz) -> Result<Self, MathError> {
        if (unsafe { fmpz_cmp_si(value, 1) } <= 0) {
            let z = Z::from(value);
            return Err(MathError::InvalidIntToModulus(z.to_string()));
        }

        let mut ctx = MaybeUninit::uninit();
        unsafe {
            fmpz_mod_ctx_init(ctx.as_mut_ptr(), value);
            Ok(Self {
                modulus: Rc::new(ctx.assume_init()),
            })
        }
    }
}

// TODO: write macro to generate [`TryFrom`] trait.
impl TryFrom<&Z> for Modulus {
    type Error = MathError;
    /// Create [`Modulus`] from [`Z`] reference using [`try_from_z`](Modulus::try_from_z)
    fn try_from(value: &Z) -> Result<Self, Self::Error> {
        Modulus::try_from_z(value)
    }
}

impl FromStr for Modulus {
    type Err = MathError;

    /// Create a [`Modulus`] from a string with a decimal number.
    ///
    /// Parameters:
    /// - `s`: the modulus of form: "[0-9]+" and not all zeros
    ///
    /// Returns a [`Modulus`] or an [`MathError`], if the provided string is not
    /// a valid modulus.
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::Modulus;
    /// use std::str::FromStr;
    ///
    /// let modulus = Modulus::from_str("42").unwrap();
    /// ```
    /// # Errors and Failures
    ///
    /// - Returns a [`MathError`] of type
    /// [`InvalidStringToZInput`](MathError::InvalidStringToZInput) if the
    /// provided string was not formatted correctly, e.g. not a correctly
    /// formatted [`Z`].
    /// - Returns a [`MathError`] of type
    /// [`InvalidIntToModulus`](MathError::InvalidIntToModulus)
    /// if the provided value is not greater than `1`.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let z = Z::from_str(s)?;

        Modulus::from_fmpz_ref(&z.value)
    }
}

#[cfg(test)]
mod test_try_from_z {
    use super::Modulus;
    use crate::integer::Z;

    /// Demonstrate different ways to create a [`Modulus`] from a [`Z`] reference.
    #[test]
    fn z_to_modulus_reference() {
        let value = Z::from(10);

        let _: Modulus = (&value).try_into().unwrap();
        let _: Modulus = <&Z>::try_into(&value).unwrap();
        assert!(Modulus::try_from_z(&value).is_ok());
        assert!(Modulus::try_from(&value).is_ok());
    }

    /// Test [`Modulus`] creation with small positive [`Z`] (valid).
    #[test]
    fn z_valid_small() {
        let z = Z::from(10);

        assert!(Modulus::try_from_z(&z).is_ok());
    }

    /// Test [`Modulus`] creation with large [`Z`] (valid).
    /// (uses FLINT's pointer representation)
    #[test]
    fn z_valid_large() {
        let z = Z::from(u64::MAX);

        assert!(Modulus::try_from_z(&z).is_ok());
    }

    /// Test [`Modulus`] creation with small negative [`Z`] (invalid).
    #[test]
    fn z_negative_small() {
        let z = Z::from(-10);

        assert!(Modulus::try_from_z(&z).is_err());
    }

    /// Test [`Modulus`] creation with large negative [`Z`] (invalid).
    /// (uses FLINT's pointer representation)
    #[test]
    fn z_negative_large() {
        let z = Z::from(i64::MIN);

        assert!(Modulus::try_from_z(&z).is_err());
    }
}

#[cfg(test)]
mod test_from_fmpz_ref {
    use super::*;
    use crate::integer::Z;

    /// Tests whether a small modulus ist instantiated correctly.
    #[test]
    fn working_example() {
        let z = Z::from(100);
        assert!(Modulus::from_fmpz_ref(&z.value).is_ok())
    }

    /// Tests whether a large modulus (> 64 bits) is instantiated correctly.
    #[test]
    fn large_modulus() {
        let z = &Z::from(u64::MAX);
        assert!(Modulus::from_fmpz_ref(&z.value).is_ok());
    }

    /// Tests whether a negative input value returns an error.
    #[test]
    fn negative_modulus() {
        let z = &Z::from(-42);
        assert!(Modulus::from_fmpz_ref(&z.value).is_err())
    }

    /// Tests whether a zero as input value returns an error.
    #[test]
    fn zero_modulus() {
        let z = &Z::ZERO;
        assert!(Modulus::from_fmpz_ref(&z.value).is_err())
    }
}

#[cfg(test)]
mod test_from_str {

    use super::Modulus;
    use std::str::FromStr;

    /// tests whether a correctly formatted string outputs an instantiation of a
    /// Modulus, i.e. does not return an error
    #[test]
    fn working_example() {
        assert!(Modulus::from_str("42").is_ok());
    }

    /// tests whether a large value (> 64 bits) is instantiated correctly
    #[test]
    fn large_value() {
        assert!(Modulus::from_str(&"1".repeat(65)).is_ok())
    }

    /// tests whether a falsely formatted string (wrong whitespaces) returns an
    /// error
    #[test]
    fn false_format_whitespaces() {
        assert!(Modulus::from_str("4 2").is_err());
    }

    /// tests whether a falsely formatted string (wrong symbols) returns an error
    #[test]
    fn false_format_symbols() {
        assert!(Modulus::from_str("b a").is_err());
    }

    /// tests whether a false string (negative) returns an error
    #[test]
    fn false_sign() {
        assert!(Modulus::from_str("-42").is_err());
    }
}
