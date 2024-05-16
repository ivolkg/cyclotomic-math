// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-math.
//
// qFALL-math is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! [`PolynomialRingZq`] is a type of ring over PolyOverZq/f(X).
//! Where f(X) is a [`PolyOverZq`](crate::integer_mod_q::PolyOverZq).
//! This implementation uses the [FLINT](https://flintlib.org/) library.
//!
//! For **DEVELOPERS**: Many functions assume that the [`PolynomialRingZq`] instances are reduced.
//! To avoid unnecessary checks and reductions, always return canonical/reduced
//! values. The end-user should be unable to obtain a non-reduced value.
//! Therefore, the DEVELOPER has to call the [`PolynomialRingZq::reduce`], whenever
//! a computation may exceed the modulus, because it is not reduced automatically

use super::{ModulusPolynomialRingZq, Zq};
use crate::{
    integer::{PolyOverZ, Z},
    traits::{GetCoefficient, Pow},
};
use derive_more::Display;
use serde::{Deserialize, Serialize};

mod arithmetic;
mod from;
mod get;
mod reduce;
mod sample;
mod set;

/// [`PolynomialRingZq`] represents polynomials over the finite field
/// [`PolyOverZq`](crate::integer_mod_q::PolyOverZq)/f(X) where f(X) is a polynomial over [`Zq`](super::Zq).
///
/// Attributes
/// - `poly`: holds the value
/// - `modulus`: holds the modulus q and f(X)
///
/// # Examples
/// ```
/// # use qfall_math::error::MathError;
/// use qfall_math::integer::PolyOverZ;
/// use qfall_math::integer_mod_q::ModulusPolynomialRingZq;
/// use qfall_math::integer_mod_q::PolyOverZq;
/// use qfall_math::integer_mod_q::PolynomialRingZq;
/// use std::str::FromStr;
///
/// let poly_mod = PolyOverZq::from_str("3  1 0 1 mod 17").unwrap();
/// let modulus = ModulusPolynomialRingZq::from(poly_mod);
///
/// // instantiation
/// let a = PolynomialRingZq::from((PolyOverZ::from(5), &modulus));
/// let b = PolynomialRingZq::from((PolyOverZ::from_str("2  1 5").unwrap(), &modulus));
/// let _ = a.clone();
///
/// // arithmetics
/// let _ = &a + &b;
/// let q_ = &a * &b;
///
/// // to_string incl. (de-)serialization
/// assert_eq!("1  5 / 3  1 0 1 mod 17", &a.to_string());
/// let _ = serde_json::to_string(&a).unwrap();
///
/// # Ok::<(), MathError>(())
/// ```
#[derive(PartialEq, Eq, Debug, Serialize, Deserialize, Display, Clone)]
#[display(fmt = "{poly} / {modulus}")]
pub struct PolynomialRingZq {
    pub(crate) poly: PolyOverZ,
    pub(crate) modulus: ModulusPolynomialRingZq,
}

impl PolynomialRingZq {
    fn evaluate(&self, prime: usize, prime_power: usize, rou: Zq) -> PolynomialRingZqNTTBasis {
        let evals = self.ntt(prime, prime_power, rou);
        PolynomialRingZqNTTBasis { evals }
    }
    fn ntt(&self, prime: usize, prime_power: usize, rou: Zq) -> Vec<Zq> {
        // Currently only using prime power decomposition, as we might not need decomposition for all factors in NTT
        // get rou from ModulusPoly?
        let q = rou.get_mod().to_string().parse::<u32>().unwrap();
        let domain_size = prime.pow(prime_power as u32) as u32;
        let resized_power = q / domain_size;
        let omega = rou.pow(resized_power).unwrap();
        let mut current_omega_power = Zq::from_z_modulus(&Z::from(1), q);
        let mut omega_powers = Vec::new();
        for i in 0..prime.pow(prime_power as u32) {
            omega_powers.push(current_omega_power.clone());
            current_omega_power = &current_omega_power * &omega;
        }

        let poly = &self.poly;
        let mut coeffs = (0..domain_size) // Unclear if the coeficients are order from the least power to the greater, probably is
            .map(|i| {
                let coeff_z = poly.get_coeff(i).unwrap();
                Zq::from_z_modulus(&coeff_z, q)
            })
            .collect::<Vec<_>>();
        radixp_ntt(prime, prime_power, &omega_powers, &mut coeffs);

        coeffs
    }
}

fn radixp_ntt(prime: usize, prime_power: usize, omega_powers: &[Zq], coeffs: &mut [Zq]) {
    let n = coeffs.len();
    assert_eq!(prime.pow(prime_power as u32), n);
    let q = omega_powers.first().unwrap().get_mod();
    if n == 1 {
        return; // This might need to change for p instead of 2
    }
    let zero = Zq::from_z_modulus(&Z::from(0), q);
    let mut decomposed_coeffs = (0..prime)
        .map(|_| vec![zero.clone(); n / prime])
        .collect::<Vec<_>>();

    for i in 0..n / prime {
        for j in 0..decomposed_coeffs.len() {
            decomposed_coeffs[j][i] = coeffs[prime * i + j].clone();
        }
    }

    let primed_omegas = omega_powers
        .iter()
        .step_by(prime)
        .map(|w| w.clone())
        .collect::<Vec<_>>();

    for i in 0..prime {
        radixp_ntt(
            prime,
            prime_power - 1,
            &primed_omegas,
            &mut decomposed_coeffs[i],
        );
    }

    for q in 0..n / prime {
        // review and optimize
        for s in 0..prime {
            let mut sum = zero.clone();
            for l in 0..prime {
                sum =
                    sum + &omega_powers[(l * (q + s * (n / prime))) % n] * &decomposed_coeffs[l][q];
            }
            coeffs[q + s * (n / prime)] = sum;
        }
    }
}

pub struct PolynomialRingZqNTTBasis {
    pub(crate) evals: Vec<Zq>,
}

impl PolynomialRingZqNTTBasis {
    fn intt(&self, prime: usize, prime_power: usize, rou: Zq) -> Vec<Zq> {
        let q = rou.get_mod().to_string().parse::<u32>().unwrap(); //compute domain outside
        let domain_size = prime.pow(prime_power as u32) as u32;
        let resized_power = q / domain_size;
        let omega = rou.pow(resized_power).unwrap();
        let mut current_omega_power = Zq::from_z_modulus(&Z::from(1), q);
        let mut omega_powers = Vec::new();
        for i in 0..prime.pow(prime_power as u32) {
            omega_powers.push(current_omega_power.clone());
            current_omega_power = &current_omega_power * &omega;
        }

        // inverse domain
        omega_powers.reverse();
        omega_powers.rotate_right(1);

        let mut evals = self.evals.clone();
        radixp_ntt(prime, prime_power, &omega_powers, &mut evals);

        let nq_inv = Zq::from_z_modulus(&Z::from(domain_size), q)
            .inverse()
            .unwrap();

        let coeffs = evals.iter_mut().map(|eval| &*eval * &nq_inv).collect();

        coeffs
    }
}

pub struct PolynomialRingZqCRTBasis {
    pub(crate) coeffs: Vec<Zq>,
}
