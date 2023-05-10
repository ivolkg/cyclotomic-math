// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implementations to create a [`PolyOverZ`] value from other types.
//! For each reasonable type, an explicit function with the format
//! `from_<type_name>` and the [`From`] trait should be implemented.
//!
//! The explicit functions contain the documentation.

use super::PolyOverZ;
use crate::{
    error::MathError,
    integer_mod_q::{PolyOverZq, PolynomialRingZq},
};
use flint_sys::{fmpz_mod_poly::fmpz_mod_poly_get_fmpz_poly, fmpz_poly::fmpz_poly_set_str};
use std::{ffi::CString, str::FromStr};

impl FromStr for PolyOverZ {
    type Err = MathError;

    // TODO: the second whitespace is not shown in the Rust-documentation
    /// Create a new polynomial with arbitrarily many coefficients of type
    /// [`Z`](crate::integer::Z).
    ///
    /// Parameters:
    /// - `s`: the polynomial of form: `"[#number of coefficients]⌴⌴[0th coefficient]⌴[1st coefficient]⌴..."`.
    ///
    /// Note that the `[#number of coefficients]` and `[0th coefficient]`
    /// are divided by two spaces and the input string is trimmed, i.e. all whitespaces
    /// before and after are removed.
    ///
    /// Returns a [`PolyOverZ`] or an error, if the provided string was not formatted
    /// correctly.
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer::PolyOverZ;
    /// use std::str::FromStr;
    ///
    /// let poly = PolyOverZ::from_str("4  0 1 2 3").unwrap();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`MathError::InvalidStringToPolyInput`]
    /// if the provided string was not formatted correctly or the number of
    /// coefficients was smaller than the number provided at the start of the
    /// provided string.
    /// - Returns a [`MathError`] of type
    /// [`InvalidStringToPolyMissingWhitespace`](MathError::InvalidStringToPolyMissingWhitespace)
    /// if the provided value did not contain two whitespaces.
    /// - Returns a [`MathError`] of type
    /// [`InvalidStringToCStringInput`](MathError::InvalidStringToCStringInput)
    /// if the provided string contains a Null Byte.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // remove whitespaces at the start and at the end
        let s_trimmed = s.trim();

        if s_trimmed == "0" {
            return Ok(Self::default());
        }

        // fmpz_poly_set_str just skips the two symbols after the first space
        // behind the number of coefficients (even if not a space), hence
        // it has to be checked here to ensure that no number is lost.
        // We only have to check it once, because for every other position it checks
        // whether there is only one space.
        if !s_trimmed.contains("  ") {
            return Err(MathError::InvalidStringToPolyMissingWhitespace(
                s.to_owned(),
            ));
        };

        let mut res = Self::default();

        let c_string = CString::new(s_trimmed)?;

        match unsafe { fmpz_poly_set_str(&mut res.poly, c_string.as_ptr()) } {
            0 => Ok(res),
            _ => Err(MathError::InvalidStringToPolyInput(s.to_owned())),
        }
    }
}

impl PolyOverZ {
    /// Create a [`PolyOverZ`] from a [`PolyOverZq`].
    ///
    /// Parameters:
    /// - `poly`: the polynomial from which the coefficients are copied
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer::PolyOverZ;
    /// use qfall_math::integer_mod_q::PolyOverZq;
    /// use std::str::FromStr;
    ///
    /// let poly = PolyOverZq::from_str("4  0 1 102 3 mod 101").unwrap();
    ///
    /// let poly_z = PolyOverZ::from_poly_over_zq(&poly);
    ///
    /// # let cmp_poly = PolyOverZ::from_str("4  0 1 1 3").unwrap();
    /// # assert_eq!(cmp_poly, poly_z);
    /// ```
    pub fn from_poly_over_zq(poly: &PolyOverZq) -> Self {
        let mut out = Self::default();
        unsafe {
            fmpz_mod_poly_get_fmpz_poly(
                &mut out.poly,
                &poly.poly,
                poly.modulus.get_fmpz_mod_ctx_struct(),
            )
        };
        out
    }

    /// Create a [`PolyOverZ`] from a [`PolynomialRingZq`].
    ///
    /// Parameters:
    /// - `poly_ring`: the polynomial from which the coefficients are copied
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer::PolyOverZ;
    /// use qfall_math::integer_mod_q::{PolynomialRingZq, ModulusPolynomialRingZq};
    /// use std::str::FromStr;
    ///
    /// let modulus = ModulusPolynomialRingZq::from_str("4  1 0 0 1 mod 17").unwrap();
    /// let poly = PolyOverZ::from_str("4  -1 0 1 1").unwrap();
    /// let poly_ring = PolynomialRingZq::from((&poly, &modulus));
    ///
    /// let poly_z = PolyOverZ::from_polynomial_ring_zq(&poly_ring);
    ///
    /// # let cmp_poly = PolyOverZ::from_str("3  15 0 1").unwrap();
    /// # assert_eq!(cmp_poly, poly_z);
    /// ```
    pub fn from_polynomial_ring_zq(poly_ring: &PolynomialRingZq) -> Self {
        poly_ring.poly.clone()
    }
}

impl From<&PolyOverZq> for PolyOverZ {
    /// Converts a polynomial of type [`PolyOverZq`] to a [`PolyOverZ`] using
    /// [`PolyOverZ::from_poly_over_zq`].
    fn from(poly: &PolyOverZq) -> Self {
        Self::from_poly_over_zq(poly)
    }
}

impl From<&PolynomialRingZq> for PolyOverZ {
    /// Converts a polynomial ring element
    ///  of type [`PolynomialRingZq`] to a [`PolyOverZ`] using
    /// [`PolyOverZ::from_polynomial_ring_zq`].
    fn from(poly_ring: &PolynomialRingZq) -> Self {
        Self::from_polynomial_ring_zq(poly_ring)
    }
}

#[cfg(test)]
mod test_from_str {
    use super::PolyOverZ;
    use std::str::FromStr;

    /// Ensure that zero-coefficients are reduced
    #[test]
    fn reduce_zero_coeff() {
        let one_1 = PolyOverZ::from_str("2  24 1").unwrap();
        let one_2 = PolyOverZ::from_str("3  24 1 0").unwrap();

        assert_eq!(one_1, one_2)
    }

    /// tests whether the same string yields the same polynomial
    #[test]
    fn same_string() {
        let string1 = format!("3  1 {} {}", u64::MAX, i64::MIN);

        let poly_1 = PolyOverZ::from_str(&string1).unwrap();
        let poly_2 = PolyOverZ::from_str(&string1).unwrap();

        assert_eq!(poly_1, poly_2)
    }

    /// tests whether a correctly formatted string outputs an instantiation of a
    /// polynomial, i.e. does not return an error
    #[test]
    fn working_example() {
        assert!(PolyOverZ::from_str("3  1 2 -3").is_ok());
    }

    /// tests whether a falsely formatted string (missing double-space) returns
    /// an error
    #[test]
    fn missing_whitespace() {
        assert!(PolyOverZ::from_str("3 12 2 -3").is_err());
        assert!(PolyOverZ::from_str("2 17 42").is_err());
        assert!(PolyOverZ::from_str("2 17  42").is_err());
        assert!(PolyOverZ::from_str("2 17 42  ").is_err());
        assert!(PolyOverZ::from_str("  2 17 42").is_err());
    }

    /// tests whether a falsely formatted string (too many whitespaces) returns
    /// an error
    #[test]
    fn too_many_whitespaces() {
        assert!(PolyOverZ::from_str("3  1  2  -3").is_err());
    }

    /// tests whether a falsely formatted string (wrong number of total
    /// coefficients) returns an error
    #[test]
    fn false_number_of_coefficient() {
        assert!(PolyOverZ::from_str("4  1 2 -3").is_err());
    }

    /// ensure that the input works with strings that have to be trimmed
    #[test]
    fn trim_input() {
        let poly = PolyOverZ::from_str("                   4  1 2 3 -4                  ");
        assert!(poly.is_ok());
        assert_eq!(PolyOverZ::from_str("4  1 2 3 -4").unwrap(), poly.unwrap());
    }
}

#[cfg(test)]
mod test_from_poly_over_zq {
    use crate::{integer::PolyOverZ, integer_mod_q::PolyOverZq};
    use std::str::FromStr;

    /// ensure that the conversion works with positive large entries
    #[test]
    fn large_positive() {
        let poly = PolyOverZq::from_str(&format!("4  0 1 102 {} mod {}", u64::MAX - 58, u64::MAX))
            .unwrap();

        let poly_q = PolyOverZ::from(&poly);

        let cmp_poly =
            PolyOverZ::from_str(&format!("4  0 1 102 {} mod {}", u64::MAX - 58, u64::MAX)).unwrap();
        assert_eq!(cmp_poly, poly_q);
    }
}

#[cfg(test)]
mod from_polynomial_ring_zq {
    use crate::{
        integer::PolyOverZ,
        integer_mod_q::{ModulusPolynomialRingZq, PolynomialRingZq},
    };
    use std::str::FromStr;

    /// ensure that the conversion works with positive large entries
    #[test]
    fn large_positive() {
        let modulus =
            ModulusPolynomialRingZq::from_str(&format!("4  1 0 0 1 mod {}", u64::MAX - 58))
                .unwrap();
        let poly = PolyOverZ::from_str("4  -1 0 1 1").unwrap();
        let poly_ring = PolynomialRingZq::from((&poly, &modulus));

        let poly_z = PolyOverZ::from(&poly_ring);

        let cmp_poly = PolyOverZ::from_str(&format!("3  {} 0 1", u64::MAX - 60)).unwrap();
        assert_eq!(cmp_poly, poly_z);
    }
}
