// Copyright © 2023 Marcel Luca Schmidt, Sven Moog
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implementations to create a [`Q`] value from other types.
//! For each reasonable type, an explicit function with the format
//! `from_<type_name>` and the [`From`] trait should be implemented.
//!
//! The explicit functions contain the documentation.

use super::Q;
use crate::{
    error::MathError,
    integer::Z,
    macros::from::{from_trait, from_type},
};
use flint_sys::{
    fmpq::{fmpq, fmpq_canonicalise, fmpq_clear, fmpq_set_str},
    fmpz::{fmpz_is_zero, fmpz_set, fmpz_swap},
};
use fraction::Fraction;
use std::{ffi::CString, str::FromStr};

impl FromStr for Q {
    type Err = MathError;

    /// Create a [`Q`] rational from a [`String`]
    /// In the string should be two decimal numbers separated by `/`.
    /// Optionally, before one or both of them can be a `-`.
    /// The format of that string looks like this `-12/53`.
    ///
    /// If the number is an integer, the string can be in the format of one.
    /// The format of that string looks like this `-12`.
    /// It is automatically transformed to `-12/1`.
    ///
    /// Parameters:
    /// - `s`: the rational value
    ///
    /// Returns a [`Q`] or an error, if the provided string was not formatted
    /// correctly.
    ///
    /// # Examples
    /// ```
    /// use std::str::FromStr;
    /// use qfall_math::rational::Q;
    ///  
    /// let a: Q = "100/3".parse().unwrap();
    /// let b: Q = Q::from_str("100/3").unwrap();
    /// ```
    ///
    /// ```
    /// use std::str::FromStr;
    /// use qfall_math::rational::Q;
    ///  
    /// let q: Q = Q::from_str("-10/3").unwrap();
    /// let b: Q = Q::from_str("10/-3").unwrap();
    /// ```
    ///
    /// ```
    /// use std::str::FromStr;
    /// use qfall_math::rational::Q;
    ///  
    /// let q: Q = Q::from_str("-10").unwrap();
    /// let b: Q = Q::from_str("10").unwrap();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type
    /// [`InvalidStringToQInput`](MathError::InvalidStringToQInput)
    /// if the provided string was not formatted correctly.
    /// - Returns a [`MathError`] of type
    /// [`DivisionByZeroError`](MathError::DivisionByZeroError)
    /// if the provided string has `0` as the denominator.
    fn from_str(s: &str) -> Result<Self, MathError> {
        if s.contains(char::is_whitespace) {
            return Err(MathError::InvalidStringToQInput(s.to_owned()));
        }

        // `fmpq::default()` returns the value '0/0'
        let mut value = fmpq::default();

        let c_string = CString::new(s)?;

        // -1 is returned if the string is an invalid input.
        //
        // Given the documentation `c_string.as_ptr()` is freed once c_string is deallocated
        // 'The pointer will be valid for as long as `self` is'
        // For reading more look at the documentation of `.as_ptr()`.
        //
        // since value is set to `0`, if an error occurs, we do not need to free
        // the allocated space manually
        if -1 == unsafe { fmpq_set_str(&mut value, c_string.as_ptr(), 10) } {
            return Err(MathError::InvalidStringToQInput(s.to_owned()));
        };

        // canonical form is expected by other functions
        unsafe { fmpq_canonicalise(&mut value) };

        // if `value.den` is set to `0`, `value.num` is not necessarily `0` as well.
        // hence we do need to free the allocated space of the numerator
        // manually by using `fmpq_clear`
        match unsafe { fmpz_is_zero(&value.den) } {
            0 => Ok(Q { value }),
            _ => {
                unsafe {
                    fmpq_clear(&mut value);
                }
                Err(MathError::DivisionByZeroError(s.to_owned()))
            }
        }
    }
}

impl Q {
    /// Create a [`Q`] from two references that can be converted to [`Z`].
    /// For example, [`&Z`].
    ///
    /// Warning: The interface of this function will likely change in the future
    /// or it will entirely be removed. Therefore use [`Q::try_from`] instead.
    ///
    /// Parameters:
    /// - `numerator` of the new [`Q`].
    /// - `denominator` of the new [`Q`].
    ///
    /// Returns a [`Q`] or a [`MathError`]
    ///
    /// # Example
    /// ```ignore (private function)
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer::Z;
    ///
    /// let num = Z::from(100);
    ///
    /// let a = Q::try_from_int_int(&num, &i64::MAX).unwrap();
    /// let b = Q::try_from_int_int(&num, &i64::MAX).unwrap();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`DivisionByZeroError`](MathError::DivisionByZeroError)
    /// if the denominator is zero.
    fn try_from_int_int(
        numerator: &(impl Into<Z> + Clone),
        denominator: &(impl Into<Z> + Clone),
    ) -> Result<Self, MathError> {
        let mut numerator: Z = numerator.to_owned().into();
        let mut denominator: Z = denominator.to_owned().into();

        // TODO: this is not as efficient as possible when passing values that
        // internally include a [`fmpz`], e.g. [`Z`] or [`Zq`].
        // In those cases it would be faster to use `fmpq_set_fmpz_frac`.
        // The best way would probably be to use `Into<fmpz>` and do the
        // performance improvements there. However, this takes longer to implement and this
        // functionality is now required for others to make progress.

        if denominator == Z::ZERO {
            return Err(MathError::DivisionByZeroError(format!(
                "{}/{}",
                numerator, denominator
            )));
        }

        let mut res = Q::default();

        unsafe {
            fmpz_swap(&mut res.value.num, &mut numerator.value);
            fmpz_swap(&mut res.value.den, &mut denominator.value);
            fmpq_canonicalise(&mut res.value);
        }
        Ok(res)
    }

    /// Create a new Integer that can grow arbitrary large.
    ///
    /// Parameters:
    /// - `value`: the initial value the integer should have
    ///
    /// Returns the new integer.
    ///
    /// # Example
    /// ```
    /// use qfall_math::integer::Z;
    /// use qfall_math::rational::Q;
    ///
    /// let m = Z::from(17);
    ///
    /// let a: Q = Q::from_int(m);
    /// let b: Q = Q::from_int(17);
    /// ```
    pub fn from_int(value: impl Into<Z>) -> Self {
        let value = value.into();
        // this efficient implementation depends on Q::default instantiating 1 as denominator
        let mut out = Q::default();
        unsafe { fmpz_set(&mut out.value.num, &value.value) }
        out
    }

    /// Create a new rational number of type [`Q`] from a [`f64`].
    ///
    /// Input parameters:
    /// - `value` : The value the rational number will have, provided as a [`f64`]
    ///
    /// Returns a [`Q`].
    ///
    /// # Example
    /// ```rust
    /// use qfall_math::rational::Q;
    ///
    /// let a: Q = Q::from_f64(0.3);
    /// let a: Q = Q::from_f64(-123.4567);
    /// ```
    pub fn from_f64(value: f64) -> Self {
        let f = Fraction::from(value);
        let sign = f
            .sign()
            .expect("Got None element instead of a fraction, may be overflow error (NaN)")
            .is_positive();
        match sign {
            true => Q::try_from((f.numer().unwrap(), f.denom().unwrap())).unwrap(),
            false => Q::try_from((f.numer().unwrap(), f.denom().unwrap())).unwrap() * Q::MINUS_ONE,
        }
    }

    from_type!(f32, f64, Q, Q::from_f64);
}

impl<T1: Into<Z> + Clone, T2: Into<Z> + Clone> TryFrom<(&T1, &T2)> for Q {
    type Error = MathError;

    /// Create a [`Q`] from two values that can be converted to [`Z`].
    /// For example, [`Z`] and [`u32`].
    ///
    /// Parameters:
    /// - `num_den_tuple`
    ///     - first value: numerator of the new [`Q`].
    ///     - second value: denominator of the new [`Q`].
    ///
    /// Returns a [`Q`] or a [`MathError`]
    ///
    /// # Example
    /// ```rust
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer::Z;
    ///
    /// let a = Q::try_from((&42, &2)).unwrap();
    /// let b = Q::try_from((&Z::from(21), &Z::from(1))).unwrap();
    /// assert_eq!(a,b);
    ///
    /// let c = Q::try_from((&10,&0)); // Division by zero
    /// assert!(c.is_err());
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`DivisionByZeroError`](MathError::DivisionByZeroError)
    /// if the denominator is zero.
    fn try_from(num_den_tuple: (&T1, &T2)) -> Result<Self, Self::Error> {
        Q::try_from_int_int(num_den_tuple.0, num_den_tuple.1)
    }
}

impl<T: Into<Z>> From<T> for Q {
    /// Create a new Integer that can grow arbitrary large.
    ///
    /// Parameters:
    /// - `value`: the initial value the integer should have
    ///
    /// Returns the new integer.
    ///
    /// # Example
    /// ```
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer::Z;
    ///
    /// let a: Q = Q::from(17);
    /// let b: Q = Q::from(Z::from(17));
    /// ```
    fn from(value: T) -> Self {
        Q::from_int(value)
    }
}

impl From<f64> for Q {
    /// Create a new rational number of type [`Q`] from a [`f64`].
    ///
    /// Input parameters:
    /// - `value` : The value the rational number will have, provided as a [`f64`]
    ///
    /// Returns a [`Q`].
    ///
    /// # Example
    /// ```rust
    /// use qfall_math::rational::Q;
    ///
    /// let a: Q = Q::from(0.3);
    /// let a: Q = Q::from(-123.4567);
    /// ```
    fn from(value: f64) -> Self {
        Q::from_f64(value)
    }
}

from_trait!(f32, Q, Q::from_f32);

#[cfg(test)]
mod tests_from_str {

    use crate::rational::Q;
    use std::str::FromStr;

    /// Ensure that initialization with large numerators and denominators works.
    #[test]
    fn max_int_positive() {
        let mut s1 = (i64::MAX).to_string();
        s1.push('/');
        s1.push_str(&(i64::MAX).to_string());

        let mut s2 = ("1/").to_string();
        s2.push_str(&(i64::MAX).to_string());

        assert!(Q::from_str(&(i64::MAX).to_string()).is_ok());
        assert!(Q::from_str(&s1).is_ok());
        assert!(Q::from_str(&s2).is_ok());
    }

    /// Ensure that initialization with large numerators and denominators
    /// (larger than i64) works.
    #[test]
    fn big_positive() {
        let mut s1 = "1".repeat(65);
        s1.push('/');
        s1.push_str(&"1".repeat(65));

        let mut s2 = ("1/").to_string();
        s2.push_str(&"1".repeat(65));

        assert!(Q::from_str(&"1".repeat(65)).is_ok());
        assert!(Q::from_str(&s1).is_ok());
        assert!(Q::from_str(&s2).is_ok());
    }

    /// Ensure that initialization with large negative numerators and
    /// denominators works.
    #[test]
    fn max_int_negative() {
        let mut s1 = (i64::MIN).to_string();
        s1.push('/');
        s1.push_str(&(i64::MIN).to_string());

        let mut s2 = ("1/").to_string();
        s2.push_str(&(i64::MIN).to_string());

        assert!(Q::from_str(&(i64::MIN).to_string()).is_ok());
        assert!(Q::from_str(&s1).is_ok());
        assert!(Q::from_str(&s2).is_ok());
    }

    /// Ensure that initialization with large negative numerators and
    /// denominators (larger than [`i64`]) works.
    #[test]
    fn big_negative() {
        let mut s1 = "-".to_string();
        s1.push_str(&"1".repeat(65));
        s1.push('/');
        s1.push_str(&"1".repeat(65));

        let mut s2 = ("-1/").to_string();
        s2.push_str(&"1".repeat(65));

        assert!(Q::from_str(&"1".repeat(65)).is_ok());
        assert!(Q::from_str(&s1).is_ok());
        assert!(Q::from_str(&s2).is_ok());
    }

    /// Ensure that an initialization with two minus works.
    #[test]
    fn no_error_both_minus() {
        assert!(Q::from_str("-3/-2").is_ok());
    }

    /// Ensure that wrong initialization yields an Error.
    #[test]
    fn error_wrong_letters() {
        assert!(Q::from_str("hbrkt35itu3gg").is_err());
    }

    /// Ensure that wrong initialization yields an Error.
    #[test]
    fn error_wrong_order() {
        assert!(Q::from_str("3/2-").is_err());
    }

    /// Ensure that wrong initialization yields an Error.
    #[test]
    fn error_two_divisions() {
        assert!(Q::from_str("3/2/4").is_err());
    }

    /// Ensure that wrong initialization yields an Error.
    #[test]
    fn error_wrong_minus() {
        assert!(Q::from_str("-3-4/2").is_err());
    }

    /// Ensure that wrong initialization yields an Error.
    #[test]
    fn error_whitespace_mid() {
        assert!(Q::from_str("876/ 543").is_err());
    }

    /// Ensure that wrong initialization yields an Error.
    #[test]
    fn error_whitespace_start() {
        assert!(Q::from_str(" 876543").is_err());
    }

    /// Ensure that wrong initialization yields an Error.
    #[test]
    fn error_whitespace_end() {
        assert!(Q::from_str("876543 ").is_err());
    }

    /// Ensure that wrong initialization yields an Error.
    #[test]
    fn error_whitespace_minus() {
        assert!(Q::from_str("- 876543").is_err());
    }

    /// Ensure that values returned by [`Q::from_str()`] are canonical.
    #[test]
    fn canonical_result() {
        let one_1 = Q::from_str("1/1").unwrap();
        let one_2 = Q::from_str("2/2").unwrap();
        let one_3 = Q::from_str("-42/-42").unwrap();

        let zero_1 = Q::from_str("0/1").unwrap();
        let zero_2 = Q::from_str("0/42").unwrap();

        assert_eq!(one_1, one_2);
        assert_eq!(one_1, one_3);
        assert_eq!(zero_1, zero_2);
    }
}

#[cfg(test)]
mod test_from_int_int {

    use crate::integer::Z;
    use crate::rational::Q;

    /// Test the different borrowed parameter types with small numerator and denominator.
    #[test]
    fn test_types_borrowed_small() {
        let numerator = 10;
        let denominator = 15;

        let q_1 = Q::try_from((&(numerator as u8), &(denominator as i8))).unwrap();
        let q_2 = Q::try_from((&(numerator as u16), &(denominator as i16))).unwrap();
        let q_3 = Q::try_from((&(numerator as u32), &(denominator))).unwrap();
        let q_4 = Q::try_from((&(numerator as u64), &(denominator as i64))).unwrap();
        let q_5 = Q::try_from((&Z::from(numerator), &Z::from(denominator))).unwrap();

        let q_6 = Q::try_from((&(numerator as i16), &(denominator as u16))).unwrap();
        let q_7 = Q::try_from((&(numerator as i16), &(denominator as i16))).unwrap();

        assert_eq!(q_1, q_2);
        assert_eq!(q_1, q_3);
        assert_eq!(q_1, q_4);
        assert_eq!(q_1, q_5);
        assert_eq!(q_1, q_6);
        assert_eq!(q_1, q_7);
    }

    /// Ensure that large parameters work (FLINT uses pointer representation).
    #[test]
    fn working_large() {
        let numerator = u64::MAX;
        let denominator = u64::MAX - 1;
        let numerator_z = Z::from(numerator);
        let denominator_z = Z::from(denominator);

        let q_1 = Q::try_from((&numerator, &denominator)).unwrap();
        let q_2 = Q::try_from((&numerator_z, &denominator_z)).unwrap();

        assert_eq!(q_1, q_2);
    }

    /// Test with zero denominator (not valid -> should lead to an error)
    #[test]
    fn divide_by_zero() {
        let numerator = 10;
        let denominator = 0;

        let new_q = Q::try_from((&numerator, &denominator));

        assert!(new_q.is_err());
    }

    /// Test with either negative denominator or numerator
    #[test]
    fn negative_small() {
        let numerator = 10;
        let denominator = -1;

        let q_1 = Q::try_from((&numerator, &denominator)).unwrap();
        let q_2 = Q::try_from((&-numerator, &-denominator)).unwrap();

        assert_eq!(q_1, q_2);
    }

    /// Ensure that the result is canonical for small parameters.
    #[test]
    fn canonical_small() {
        let numerator = 10;
        let denominator = 1;

        let q_1 = Q::try_from((&numerator, &denominator)).unwrap();
        let q_2 = Q::try_from((&-numerator, &-denominator)).unwrap();
        let q_3 = Q::try_from((&(numerator * 2), &(denominator * 2))).unwrap();

        let q_4_negative = Q::try_from((&-numerator, &denominator)).unwrap();
        let q_5_negative = Q::try_from((&numerator, &-denominator)).unwrap();

        assert_eq!(q_1, q_2);
        assert_eq!(q_1, q_3);

        assert_eq!(q_4_negative, q_5_negative);
    }

    /// Ensure that the result is canonical for large parameters.
    #[test]
    fn canonical_large() {
        let numerator = i64::MAX;
        let denominator = i64::MAX - 1;

        let numerator_z = Z::from(numerator);
        let denominator_z = Z::from(denominator);

        let q_1 = Q::try_from((&numerator, &denominator)).unwrap();
        let q_2 = Q::try_from((&-numerator, &-denominator)).unwrap();
        let q_3 = Q::try_from((&numerator_z, &denominator_z)).unwrap();
        let q_4 =
            Q::try_from((&(&numerator_z * Z::from(2)), &(&denominator_z * Z::from(2)))).unwrap();
        let q_5_negative = Q::try_from((&(&numerator_z * Z::from(-1)), &denominator_z)).unwrap();
        let q_6_negative = Q::try_from((&numerator_z, &(&denominator_z * Z::from(-1)))).unwrap();

        assert_eq!(q_1, q_2);
        assert_eq!(q_1, q_3);
        assert_eq!(q_1, q_4);
        assert_eq!(q_5_negative, q_6_negative);
    }
}

#[cfg(test)]
mod test_from_z {
    use super::Q;
    use crate::integer::Z;
    use std::str::FromStr;

    /// Ensure that the `from_int` function is available and works correctly for
    /// small and large instances of [`Z`] and structs implementing [`Into<Z>`].
    #[test]
    fn large_small_numbers() {
        let z_1 = Z::from(u64::MAX);
        let z_2 = Z::from(17);

        assert_eq!(
            Q::from_str(&u64::MAX.to_string()).unwrap(),
            Q::from_int(z_1)
        );
        assert_eq!(Q::from_str("17").unwrap(), Q::from_int(z_2));
    }

    /// Ensure that the [`From`] trait is available and works correctly for
    /// small and large instances of [`Z`].
    #[test]
    fn from_trait() {
        let z_1 = Z::from(u64::MAX);
        let z_2 = Z::from(17);

        assert_eq!(Q::from_str(&u64::MAX.to_string()).unwrap(), Q::from(z_1));
        assert_eq!(Q::from_str("17").unwrap(), Q::from(z_2));
    }

    /// Ensure that all types that can be turned into an [`Z`]
    /// can be used to instantiate a [`Q`]
    #[test]
    fn from_into_z() {
        let _ = Q::from(u8::MAX);
        let _ = Q::from(u16::MAX);
        let _ = Q::from(u32::MAX);
        let _ = Q::from(u64::MAX);

        let _ = Q::from(i8::MIN);
        let _ = Q::from(i16::MIN);
        let _ = Q::from(i32::MIN);
        let _ = Q::from(i64::MIN);
    }
}

#[cfg(test)]
mod test_from_float {
    use super::Q;
    use std::f64::consts::{E, LN_10, LN_2};

    /// Enure that the from works correctly for positive values
    #[test]
    fn positive() {
        let numerator = 150001;
        let denominator = 16;

        let value = Q::from(numerator as f64 / denominator as f64);

        let cmp = Q::try_from((&numerator, &denominator)).unwrap();
        assert_eq!(cmp, value)
    }

    /// Enure that the from works correctly for positive values
    #[test]
    fn negative() {
        let numerator = 150001;
        let denominator = -8;

        let value = Q::from(numerator as f64 / denominator as f64);

        let cmp = Q::try_from((&numerator, &denominator)).unwrap();
        assert_eq!(cmp, value)
    }

    /// Ensure that the [`From`] trait is available for [`f64`] constants
    #[test]
    fn from_trait() {
        let _ = Q::from(E);
        let _ = Q::from(LN_10);
        let _ = Q::from(LN_2);
    }

    /// test availability for [`f32`]
    #[test]
    fn from_f32_available() {
        let f: f32 = 42.17;

        let _ = Q::from(f);
        let _ = Q::from_f32(f);
    }
}
