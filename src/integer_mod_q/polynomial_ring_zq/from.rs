// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implementations to create a [`PolynomialRingZq`] value from other types.
//! For each reasonable type, an explicit function with the format
//! `from_<type_name>` and the [`From`] trait should be implemented.
//!
//! The explicit functions contain the documentation.

use super::PolynomialRingZq;
use crate::{integer::PolyOverZ, integer_mod_q::ModulusPolynomialRingZq};

impl From<(&PolyOverZ, &ModulusPolynomialRingZq)> for PolynomialRingZq {
    /// Create a new polynomial ring element of type [`PolynomialRingZq`].
    ///
    /// Parameters:
    /// - `value` is a tuple `(polynomial, modulus)`
    ///   The modulus defines the ring.
    ///
    /// Returns a new element inside the polynomial ring.
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::PolynomialRingZq;
    /// use qfall_math::integer_mod_q::ModulusPolynomialRingZq;
    /// use qfall_math::integer::PolyOverZ;
    /// use std::str::FromStr;
    ///
    /// let modulus = ModulusPolynomialRingZq::from_str("4  1 0 0 1 mod 17").unwrap();
    /// let poly = PolyOverZ::from_str("4  -1 0 1 1").unwrap();
    /// let poly_ring = PolynomialRingZq::from((&poly, &modulus));
    /// ```
    fn from(value: (&PolyOverZ, &ModulusPolynomialRingZq)) -> Self {
        let mut out = Self {
            poly: value.0.clone(),
            modulus: value.1.clone(),
        };
        out.reduce();
        out
    }
}

#[cfg(test)]
mod test_from_poly_over_z_modulus_polynomial_ring_zq {
    use crate::{
        integer::PolyOverZ,
        integer_mod_q::{ModulusPolynomialRingZq, PolynomialRingZq},
    };
    use std::str::FromStr;

    const BITPRIME64: u64 = u64::MAX - 58;

    /// Ensure that the modulus is applied with a large prime and large coefficients
    #[test]
    fn is_reduced_large() {
        let modulus =
            ModulusPolynomialRingZq::from_str(&format!("4  1 0 0 1 mod {}", BITPRIME64)).unwrap();

        let poly = PolyOverZ::from_str(&format!("4  {} {} 1 1", BITPRIME64 + 2, u64::MAX)).unwrap();
        let poly_ring = PolynomialRingZq::from((&poly, &modulus));

        let cmp_poly = PolyOverZ::from_str("3  1 58 1").unwrap();
        let cmp_poly_ring = PolynomialRingZq::from((&cmp_poly, &modulus));

        assert_eq!(poly_ring, cmp_poly_ring);
    }

    /// Ensure that two ring elements that are instantiated the same way are equal
    #[test]
    fn same_instantiation() {
        let modulus =
            ModulusPolynomialRingZq::from_str(&format!("4  1 0 0 1 mod {}", BITPRIME64)).unwrap();
        let poly = PolyOverZ::from_str(&format!("4  {} {} 1 1", BITPRIME64 + 2, u64::MAX)).unwrap();

        let poly_ring_1 = PolynomialRingZq::from((&poly, &modulus));
        let poly_ring_2 = PolynomialRingZq::from((&poly, &modulus));

        assert_eq!(poly_ring_1, poly_ring_2);
    }
}
