//! Implementations to evaluate a [`PolyOverZ`].
//! For each reasonable type, an implementation
//! of the [`Evaluate`] trait should be implemented.

use super::PolyOverZ;
use crate::{
    integer::Z, macros::for_others::implement_for_others_evaluate, rational::Q, traits::Evaluate,
};
use flint_sys::fmpz_poly::{fmpz_poly_evaluate_fmpq, fmpz_poly_evaluate_fmpz};

impl Evaluate<Z, Z> for PolyOverZ {
    /// Evaluates a [`PolyOverZ`] on a given input.
    ///
    /// Parameters:
    /// - `value`: the value with which to evaluate the polynomial.
    ///
    /// Returns the evaluation of the polynomial as a [`Z`].
    ///
    /// # Examples
    /// ```rust
    /// use math::traits::Evaluate;
    /// use math::integer::Z;
    /// use math::integer::PolyOverZ;
    /// use std::str::FromStr;
    ///
    /// let poly = PolyOverZ::from_str("5  0 1 2 -3 1").unwrap();
    /// let value = Z::from(3);
    /// let res = poly.evaluate(value);
    /// ```
    ///
    /// ```rust
    /// use math::traits::Evaluate;
    /// use math::integer::PolyOverZ;
    /// use std::str::FromStr;
    ///
    /// let poly = PolyOverZ::from_str("5  0 1 2 -3 1").unwrap();
    /// let value = 3;
    /// let res = poly.evaluate(value);
    /// ```
    ///
    /// ```rust
    /// use math::traits::Evaluate;
    /// use math::integer::PolyOverZ;
    /// use std::str::FromStr;
    ///
    /// let poly = PolyOverZ::from_str("5  0 1 2 -3 1").unwrap();
    /// let value = i64::MAX;
    /// let res = poly.evaluate(value);
    /// ```
    fn evaluate(&self, value: Z) -> Z {
        self.evaluate(&value)
    }
}

impl Evaluate<&Z, Z> for PolyOverZ {
    /// Evaluates a [`PolyOverZ`] on a given input of [`Z`]. Note that the
    /// [`Z`] in this case is only a reference.
    ///
    /// Parameters:
    /// - `value`: the value with which to evaluate the polynomial.
    ///
    /// Returns the evaluation of the polynomial as a [`Z`].
    ///
    /// # Example
    /// ```rust
    /// use math::traits::Evaluate;
    /// use math::integer::Z;
    /// use math::integer::PolyOverZ;
    /// use std::str::FromStr;
    ///
    /// let poly = PolyOverZ::from_str("5  0 1 2 -3 1").unwrap();
    /// let value = Z::from(3);
    /// let res = poly.evaluate(&value);
    /// ```
    fn evaluate(&self, value: &Z) -> Z {
        let mut res = Z::default();

        unsafe { fmpz_poly_evaluate_fmpz(&mut res.value, &self.poly, &value.value) };

        res
    }
}

impl Evaluate<Q, Q> for PolyOverZ {
    /// Evaluates a [`PolyOverZ`] on a given input of [`Q`].
    ///
    /// Parameters:
    /// - `value`: the value with which to evaluate the polynomial.
    ///
    /// Returns the evaluation of the polynomial as a [`Q`].
    ///
    /// # Example
    /// ```rust
    /// use math::traits::Evaluate;
    /// use math::rational::Q;
    /// use math::integer::PolyOverZ;
    /// use std::str::FromStr;
    ///
    /// let poly = PolyOverZ::from_str("5  0 1 2 -3 1").unwrap();
    /// let value = Q::from_str("3/2").unwrap();
    /// let res = poly.evaluate(value);
    /// ```
    fn evaluate(&self, value: Q) -> Q {
        self.evaluate(&value)
    }
}

impl Evaluate<&Q, Q> for PolyOverZ {
    /// Evaluates a [`PolyOverZ`] on a given input of [`Q`]. Note that the
    /// [`Q`] in this case is only a reference.
    ///
    /// Parameters:
    /// - `value`: the value with which to evaluate the polynomial.
    ///
    /// Returns the evaluation of the polynomial as a [`Q`].
    ///
    /// # Example
    /// ```
    /// use math::traits::Evaluate;
    /// use math::rational::Q;
    /// use math::integer::PolyOverZ;
    /// use std::str::FromStr;
    ///
    /// let poly = PolyOverZ::from_str("5  0 1 2 -3 1").unwrap();
    /// let value = Q::from_str("3/2").unwrap();
    /// let res = poly.evaluate(&value);
    /// ```
    fn evaluate(&self, value: &Q) -> Q {
        let mut res = Q::default();

        unsafe { fmpz_poly_evaluate_fmpq(&mut res.value, &self.poly, &value.value) };

        res
    }
}

implement_for_others_evaluate!(i64, Z, Z, PolyOverZ);
implement_for_others_evaluate!(i32, Z, Z, PolyOverZ);
implement_for_others_evaluate!(i16, Z, Z, PolyOverZ);
implement_for_others_evaluate!(i8, Z, Z, PolyOverZ);
implement_for_others_evaluate!(u64, Z, Z, PolyOverZ);
implement_for_others_evaluate!(u32, Z, Z, PolyOverZ);
implement_for_others_evaluate!(u16, Z, Z, PolyOverZ);
implement_for_others_evaluate!(u8, Z, Z, PolyOverZ);

#[cfg(test)]
mod test_evaluate_z {

    use crate::integer::{PolyOverZ, Z};
    use crate::traits::Evaluate;
    use std::str::FromStr;

    /// tests if evaluate works for [`Z`] as input
    #[test]
    fn eval_z() {
        let poly = PolyOverZ::from_str("2  1 2").unwrap();

        let res = poly.evaluate(Z::from(3));

        assert_eq!(Z::from(7), res)
    }

    /// tests if evaluate with a reference works
    #[test]
    fn eval_z_ref() {
        let poly = PolyOverZ::from_str("2  1 2").unwrap();

        let res = poly.evaluate(&Z::from(3));

        assert_eq!(Z::from(7), res)
    }

    /// tests if evaluate works with negative values
    #[test]
    fn eval_z_negative() {
        let poly = PolyOverZ::from_str("2  1 2").unwrap();

        let res = poly.evaluate(&Z::from(-5));

        assert_eq!(Z::from(-9), res)
    }

    /// tests if evaluate works with large integers
    #[test]
    fn eval_z_large() {
        let poly = PolyOverZ::from_str("2  1 2").unwrap();

        let res = poly.evaluate(&Z::from_str(&"1".repeat(65)).unwrap());

        let mut res_cmp_str = "2".repeat(64);
        res_cmp_str.push('3');
        assert_eq!(Z::from_str(&res_cmp_str).unwrap(), res)
    }

    /// test if evaluate works with max of [`i64`], [`i32`], ...
    #[test]
    fn eval_max() {
        let poly = PolyOverZ::from_str("2  1 2").unwrap();

        // signed
        let _ = poly.evaluate(i64::MAX);
        let _ = poly.evaluate(i32::MAX);
        let _ = poly.evaluate(i16::MAX);
        let _ = poly.evaluate(i8::MAX);

        //unsigned
        let _ = poly.evaluate(u64::MAX);
        let _ = poly.evaluate(u32::MAX);
        let _ = poly.evaluate(u16::MAX);
        let _ = poly.evaluate(u8::MAX);
    }

    /// test if evaluate works with min of [`i64`], [`i32`], ...
    #[test]
    fn eval_min() {
        let poly = PolyOverZ::from_str("2  1 2").unwrap();

        // signed
        let _ = poly.evaluate(i64::MIN);
        let _ = poly.evaluate(i32::MIN);
        let _ = poly.evaluate(i16::MIN);
        let _ = poly.evaluate(i8::MIN);

        // unsigned
        let _ = poly.evaluate(u64::MIN);
        let _ = poly.evaluate(u32::MIN);
        let _ = poly.evaluate(u16::MIN);
        let _ = poly.evaluate(u8::MIN);
    }
}

#[cfg(test)]
mod test_evaluate_q {
    use crate::{integer::PolyOverZ, rational::Q, traits::Evaluate};
    use std::str::FromStr;

    /// ensures that positive values return expected evaluation
    #[test]
    fn evaluate_positive() {
        let poly = PolyOverZ::from_str("2  1 3").unwrap();
        let value = Q::from_str("3/2").unwrap();

        let res_ref = poly.evaluate(&value);
        let res = poly.evaluate(value);

        assert_eq!(Q::from_str("11/2").unwrap(), res);
        assert_eq!(res_ref, res);
    }

    /// ensures that negative values return expected evaluation
    #[test]
    fn evaluate_negative() {
        let poly = PolyOverZ::from_str("2  1 3").unwrap();
        let value = Q::from_str("-3/2").unwrap();

        let res_ref = poly.evaluate(&value);
        let res = poly.evaluate(value);

        assert_eq!(Q::from_str("-7/2").unwrap(), res);
        assert_eq!(res_ref, res);
    }

    /// ensures that positive large values return expected evaluation
    #[test]
    fn evaluate_large_positive() {
        let poly = PolyOverZ::from_str(&format!("2  {} 1", (u64::MAX - 1) / 2)).unwrap();
        let value = Q::from_str(&((u64::MAX - 1) / 2).to_string()).unwrap();

        let res_ref = poly.evaluate(&value);
        let res = poly.evaluate(value);

        assert_eq!(Q::from_str(&(u64::MAX - 1).to_string()).unwrap(), res);
        assert_eq!(res_ref, res);
    }

    /// ensures that negative large values return expected evaluation
    #[test]
    fn evaluate_large_negative() {
        let poly = PolyOverZ::from_str(&format!("2  {} 2", u64::MAX)).unwrap();
        let value = Q::from_str(&format!("-{}", (u64::MAX - 1) / 2)).unwrap();

        let res_ref = poly.evaluate(&value);
        let res = poly.evaluate(value);

        assert_eq!(Q::from_str("1").unwrap(), res);
        assert_eq!(res_ref, res);
    }
}