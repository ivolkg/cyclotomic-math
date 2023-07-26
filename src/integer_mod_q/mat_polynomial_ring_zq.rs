// Copyright © 2023 Marcel Luca Schmidt
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! [`MatPolynomialRingZq`] is a type of matrix with entries of [`PolynomialRingZq`](crate::integer_mod_q::PolynomialRingZq).
//! This implementation uses the [FLINT](https://flintlib.org/) library.

use super::ModulusPolynomialRingZq;
use crate::integer::MatPolyOverZ;
use derive_more::Display;
use serde::{Deserialize, Serialize};

mod arithmetic;
mod concat;
mod default;
mod from;
mod get;
mod reduce;
mod sample;
mod set;
mod transpose;
mod vector;

/// [`MatPolynomialRingZq`] is a matrix with entries of type [`PolynomialRingZq`](crate::integer_mod_q::PolynomialRingZq).
///
/// Attributes:
/// - `matrix`: holds the [`MatPolyOverZ`](crate::integer::MatPolyOverZ) matrix
/// - `modulus` : holds the [`ModulusPolynomialRingZq`](crate::integer_mod_q::ModulusPolynomialRingZq)
/// modulus of the matrix
///
/// TODO: Add Examples
#[derive(PartialEq, Eq, Debug, Serialize, Deserialize, Display)]
#[display(fmt = "{matrix} / {modulus}")]
pub struct MatPolynomialRingZq {
    pub(crate) matrix: MatPolyOverZ,
    pub(crate) modulus: ModulusPolynomialRingZq,
}
