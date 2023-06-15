// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains algorithms for sampling according to the discrete Gaussian distribution.

use crate::{
    error::MathError,
    integer::Z,
    integer_mod_q::MatZq,
    rational::Q,
    traits::{GetNumColumns, GetNumRows, SetEntry},
    utils::sample::discrete_gauss::sample_z,
};
use std::fmt::Display;

impl MatZq {
    /// Initializes a new matrix with dimensions `num_rows` x `num_columns` and with each entry
    /// sampled independently according to the discrete Gaussian distribution,
    /// using [`Z::sample_discrete_gauss`].
    ///
    /// Parameters:
    /// - `num_rows`: specifies the number of rows the new matrix should have
    /// - `num_cols`: specifies the number of columns the new matrix should have
    /// - `n`: specifies the range from which [`Z::sample_discrete_gauss`] samples
    /// - `center`: specifies the positions of the center with peak probability
    /// - `s`: specifies the Gaussian parameter, which is proportional
    /// to the standard deviation `sigma * sqrt(2 * pi) = s`
    ///
    /// Returns a matrix with each entry sampled independently from the
    /// specified discrete Gaussian distribution.
    ///
    /// # Example
    /// ```
    /// use qfall_math::integer_mod_q::MatZq;
    ///
    /// let sample = MatZq::sample_discrete_gauss(3, 1, 83, 1024, 0, 1.25f32).unwrap();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type
    /// [`InvalidMatrix`](MathError::InvalidMatrix)
    /// if the number of rows or columns is `0`.
    /// - Returns a [`MathError`] of type [`OutOfBounds`](MathError::OutOfBounds)
    /// if the number of rows or columns is negative or it does not fit into an [`i64`].
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if the `n <= 1` or `s <= 0`.
    pub fn sample_discrete_gauss<T, T1, T2, T3>(
        num_rows: impl TryInto<i64> + Display,
        num_cols: impl TryInto<i64> + Display,
        modulus: T,
        n: T1,
        center: T2,
        s: T3,
    ) -> Result<MatZq, MathError>
    where
        T: Into<Z>,
        T1: Into<Z>,
        T2: Into<Q>,
        T3: Into<Q>,
    {
        let modulus: Z = modulus.into();
        let n: Z = n.into();
        let center: Q = center.into();
        let s: Q = s.into();
        let mut out = Self::new(num_rows, num_cols, modulus)?;

        for row in 0..out.get_num_rows() {
            for col in 0..out.get_num_columns() {
                let sample = sample_z(&n, &center, &s)?;
                out.set_entry(row, col, sample).unwrap();
            }
        }

        Ok(out)
    }
}

#[cfg(test)]
mod test_sample_discrete_gauss {
    use crate::{
        integer::Z,
        integer_mod_q::{MatZq, Modulus},
        rational::Q,
    };

    // This function only allows for a broader availability, which is tested here.

    /// Checks whether `sample_discrete_gauss` is available for all types
    /// implementing Into<Z>, i.e. u8, u16, u32, u64, i8, ...
    /// or Into<Q>, i.e. u8, i16, f32, Z, Q, ...
    #[test]
    fn availability() {
        let n = Z::from(1024);
        let center = Q::from(0);
        let s = Q::ONE;
        let modulus = Modulus::try_from(&Z::from(83)).unwrap();

        let _ = MatZq::sample_discrete_gauss(2u64, 3i8, &modulus, &16u16, 0f32, &1u16);
        let _ = MatZq::sample_discrete_gauss(3u8, 2i16, &83u8, &2u32, &center, &1u8);
        let _ = MatZq::sample_discrete_gauss(1, 1, &n, &2u64, &center, &1u32);
        let _ = MatZq::sample_discrete_gauss(1, 1, &83i8, &2i8, &center, &1u64);
        let _ = MatZq::sample_discrete_gauss(1, 1, &83, &2i16, &center, &1i64);
        let _ = MatZq::sample_discrete_gauss(1, 1, &83, &2i32, &center, &1i32);
        let _ = MatZq::sample_discrete_gauss(1, 1, &83, &2i64, &center, &1i16);
        let _ = MatZq::sample_discrete_gauss(1, 1, &83, &n, &center, &1i8);
        let _ = MatZq::sample_discrete_gauss(1, 1, &83, &2u8, &center, &1i64);
        let _ = MatZq::sample_discrete_gauss(1, 1, &83, &2, &center, &n);
        let _ = MatZq::sample_discrete_gauss(1, 1, &83, &2, &center, &s);
        let _ = MatZq::sample_discrete_gauss(1, 1, &83, 2, &center, 1.25f64);
        let _ = MatZq::sample_discrete_gauss(1, 1, &83, 2, 0, 1.25f64);
        let _ = MatZq::sample_discrete_gauss(1, 1, 83, &2, center, 15.75f32);
    }
}