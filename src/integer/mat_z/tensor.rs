// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains the implementation of the `tensor` product.

use super::MatZ;
use crate::traits::{GetNumColumns, GetNumRows, Tensor};
use flint_sys::{
    fmpz::{fmpz_is_zero, fmpz_mul},
    fmpz_mat::fmpz_mat_entry,
};

impl Tensor for MatZ {
    /// Computes the tensor product of `self` with `other`
    ///
    /// Parameters:
    /// - `other`: the value with which the tensor product is computed.
    ///
    /// Returns the tensor product of `self` with `other`.
    ///
    /// # Example
    /// ```
    /// use qfall_math::integer::MatZ;
    /// use qfall_math::traits::Tensor;
    /// use std::str::FromStr;
    ///
    /// let mat_1 = MatZ::from_str("[[1, 0, 0],[0, 1, 0],[0, 0, 1]]").unwrap();
    /// let mat_2 = MatZ::from_str("[[1, 2, 1],[3, 4, 1]]").unwrap();
    ///
    /// let mat_3 = mat_1.tensor(&mat_2);
    /// ```
    fn tensor(&self, other: &Self) -> Self {
        let columns_other = other.get_num_columns();
        let rows_other = other.get_num_rows();

        // we can unwrap since we know that the dimensions are positive
        let mut out = MatZ::new(
            self.get_num_rows() * rows_other,
            self.get_num_columns() * columns_other,
        )
        .unwrap();

        for i in 0..self.get_num_rows() {
            for j in 0..self.get_num_columns() {
                let entry = unsafe { fmpz_mat_entry(&self.matrix, i, j) };

                if unsafe { 1 != fmpz_is_zero(entry) } {
                    for i_other in 0..rows_other {
                        for j_other in 0..columns_other {
                            unsafe {
                                fmpz_mul(
                                    fmpz_mat_entry(
                                        &out.matrix,
                                        i * rows_other + i_other,
                                        j * columns_other + j_other,
                                    ),
                                    entry,
                                    fmpz_mat_entry(&other.matrix, i_other, j_other),
                                )
                            }
                        }
                    }
                }
            }
        }

        out
    }
}

#[cfg(test)]
mod test_tensor {
    use crate::{
        integer::MatZ,
        traits::{GetNumColumns, GetNumRows, Tensor},
    };
    use std::str::FromStr;

    /// ensure that the dimensions of the tensor product are taken over correctly.
    #[test]
    fn dimensions_fit() {
        let mat_1 = MatZ::new(17, 13).unwrap();
        let mat_2 = MatZ::new(3, 4).unwrap();

        let mat_3 = mat_1.tensor(&mat_2);

        assert_eq!(51, mat_3.get_num_rows());
        assert_eq!(52, mat_3.get_num_columns());
    }

    /// ensure that the tensor works correctly with identity
    #[test]
    fn identity() {
        let identity = MatZ::from_str("[[1, 0],[0, 1]]").unwrap();
        let mat_1 =
            MatZ::from_str(&format!("[[1, {}, 1],[0, {}, -1]]", u64::MAX, i64::MIN)).unwrap();

        let mat_2 = identity.tensor(&mat_1);
        let mat_3 = mat_1.tensor(&identity);

        let cmp_mat_2 = MatZ::from_str(&format!(
            "[[1, {}, 1, 0, 0, 0],[0, {}, -1, 0, 0, 0],[0, 0, 0, 1, {}, 1],[0, 0, 0, 0, {}, -1]]",
            u64::MAX,
            i64::MIN,
            u64::MAX,
            i64::MIN
        ))
        .unwrap();
        let cmp_mat_3 = MatZ::from_str(&format!(
            "[[1, 0, {}, 0, 1, 0],[0, 1, 0, {}, 0, 1],[0, 0, {}, 0, -1, 0],[0, 0, 0, {}, 0, -1]]",
            u64::MAX,
            u64::MAX,
            i64::MIN,
            i64::MIN
        ))
        .unwrap();

        assert_eq!(cmp_mat_2, mat_2);
        assert_eq!(cmp_mat_3, mat_3);
    }

    /// ensure the tensor product works where one is a vector and the other is a matrix
    #[test]
    fn vector_matrix() {
        let vector = MatZ::from_str("[[1],[-1]]").unwrap();
        let mat_1 =
            MatZ::from_str(&format!("[[1, {}, 1],[0, {}, -1]]", u64::MAX, i64::MAX)).unwrap();

        let mat_2 = vector.tensor(&mat_1);
        let mat_3 = mat_1.tensor(&vector);

        let cmp_mat_2 = MatZ::from_str(&format!(
            "[[1, {}, 1],[0, {}, -1],[-1, -{}, -1],[0, -{}, 1]]",
            u64::MAX,
            i64::MAX,
            u64::MAX,
            i64::MAX
        ))
        .unwrap();
        let cmp_mat_3 = MatZ::from_str(&format!(
            "[[1, {}, 1],[-1, -{}, -1],[0, {}, -1],[0, -{}, 1]]",
            u64::MAX,
            u64::MAX,
            i64::MAX,
            i64::MAX
        ))
        .unwrap();

        assert_eq!(cmp_mat_2, mat_2);
        assert_eq!(cmp_mat_3, mat_3);
    }

    /// ensure that the tensor product works correctly with two vectors
    #[test]
    fn vector_vector() {
        let vec_1 = MatZ::from_str("[[2],[1]]").unwrap();
        let vec_2 =
            MatZ::from_str(&format!("[[{}],[{}]]", (u64::MAX - 1) / 2, i64::MIN / 2)).unwrap();

        let vec_3 = vec_1.tensor(&vec_2);
        let vec_4 = vec_2.tensor(&vec_1);

        let cmp_vec_3 = MatZ::from_str(&format!(
            "[[{}],[{}],[{}],[{}]]",
            u64::MAX - 1,
            i64::MIN,
            (u64::MAX - 1) / 2,
            i64::MIN / 2
        ))
        .unwrap();
        let cmp_vec_4 = MatZ::from_str(&format!(
            "[[{}],[{}],[{}],[{}]]",
            u64::MAX - 1,
            (u64::MAX - 1) / 2,
            i64::MIN,
            i64::MIN / 2
        ))
        .unwrap();

        assert_eq!(cmp_vec_3, vec_3);
        assert_eq!(cmp_vec_4, vec_4);
    }
}
