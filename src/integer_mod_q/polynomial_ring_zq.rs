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

use super::{MatZq, ModulusPolynomialRingZq, Zq};
use crate::{
    integer::{PolyOverZ, Z},
    traits::{GetCoefficient, GetEntry, Pow, SetEntry},
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

    fn to_crt_basis(&self, prime: usize, prime_power: usize, rou: Zq) -> PolynomialRingZqCRTBasis {
        let poly = &self.poly;
        let m = self.get_mod().get_degree() as usize;
        let q = rou.get_mod().to_string().parse::<u32>().unwrap();
        let domain_size = prime.pow(prime_power as u32) as u32;

        let mut coeffs = (0..domain_size) // Unclear if the coeficients are order from the least power to the greater, probably is
            .map(|i| {
                let coeff_z = poly.get_coeff(i).unwrap();
                Zq::from_z_modulus(&coeff_z, q)
            })
            .collect::<Vec<_>>();
        crt(&mut coeffs, prime, prime_power, rou);
        PolynomialRingZqCRTBasis { coeffs }
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

fn crt(coeffs: &mut [Zq], prime: usize, prime_power: usize, rou: Zq) {
    // Currently only using prime power decomposition, as we might noot need decomposition for all factors in NTT
    // get rou from ModulusPoly?
    let q = rou.get_mod().to_string().parse::<u32>().unwrap();
    let m = coeffs.len(); // Double check
    let m_prime = m / prime;
    let domain_size = prime.pow(prime_power as u32) as u32;
    let resized_power = q / domain_size;
    let omega = rou.pow(resized_power).unwrap();
    let mut current_omega_power = Zq::from_z_modulus(&Z::from(1), q);
    let mut omega_powers = Vec::new();
    for i in 0..prime.pow(prime_power as u32) {
        omega_powers.push(current_omega_power.clone());
        current_omega_power = &current_omega_power * &omega;
    }

    let prime_omegas = omega_powers
        .iter()
        .step_by(m_prime)
        .map(|w| w.clone())
        .collect::<Vec<_>>();

    let crt_prime = crt_prime(&prime_omegas);
    if prime_power == 1 {
        let coeffs_mat = MatZq::new(prime - 1, 1, q); // Theres got to be  abetter way or implment Mat * Vec multiplication
        let res = crt_prime * coeffs_mat;
        for i in 0..prime - 1 {
            coeffs[i] = res.get_entry(i, 0).unwrap();
        }
        return;
    }

    stride_permutation(prime, coeffs);

    // _ variables because it should modify coeffs
    let _ = coeffs // TODO: Parallelize
        .chunks_exact_mut(prime - 1) // phi(p) = p - 1
        .map(|coeffs_chunk| {
            let mut coeffs_mat = MatZq::new(prime, 1, q);
            for (i, coeff) in coeffs_chunk.iter().enumerate() {
                coeffs_mat.set_entry(i, 0, coeff).unwrap();
            }
            let res = crt_prime.clone() * coeffs_mat;
            for i in 0..coeffs_chunk.len() {
                coeffs_chunk[i] = res.get_entry(i, 0).unwrap();
            }
        });

    inverse_stride_permutation(prime, coeffs);

    let t_hat = twiddle_hat_factors(prime, &omega_powers);
    for (twiddle_factor, coeff) in t_hat.iter().zip(coeffs.iter_mut()) {
        *coeff = &*coeff * twiddle_factor;
    }

    let m_prime_omega_powers = omega_powers
        .iter()
        .step_by(prime)
        .map(|w| w.clone())
        .collect::<Vec<_>>();

    let _ = coeffs.chunks_exact_mut(m_prime).map(|coeffs_chunk| {
        // TODO: Parallelize
        radixp_ntt(prime, prime_power - 1, &m_prime_omega_powers, coeffs_chunk);
    });

    // Note that the last stride permutation is not done
    // This might not be needed because we will perform the inverse as well
}

// Double check
fn twiddle_hat_factors(prime: usize, omega_powers: &[Zq]) -> Vec<Zq> {
    let euler_totient_p = prime - 1;
    let (euler_totient_m, m_relative_primes) = euler_totient(omega_powers.len());
    let m_prime = euler_totient_m / euler_totient_p;
    let q = omega_powers.first().unwrap().get_mod();
    let zero = Zq::from_z_modulus(&Z::from(0), q);
    let mut t_hat = vec![zero; euler_totient_m];
    for r in m_relative_primes {
        for j in 0..euler_totient_m {
            let index = (m_prime * r) + j;
            t_hat[index] = omega_powers[(r * j) % (omega_powers.len())].clone();
        }
    }
    t_hat
}

fn crt_prime(prime_omegas: &[Zq]) -> MatZq {
    let (euler_totient_m, prime_relatives) = euler_totient(prime_omegas.len());
    assert_eq!(
        prime_omegas.len() - 1,
        euler_totient_m,
        "lenght is not a prime"
    );
    let q = prime_omegas.first().unwrap().get_mod();
    let mut crt = MatZq::new(euler_totient_m, euler_totient_m, q);
    for i in 1..prime_omegas.len() {
        for j in 0..euler_totient_m {
            crt.set_entry(i, j, prime_omegas[(i * j) % prime_omegas.len()].clone())
                .unwrap();
        }
    }
    crt
}

fn inverse_crt_prime(prime_omegas: &[Zq]) -> MatZq {
    crt_prime(prime_omegas).inverse().unwrap()
}

fn euler_totient(m: usize) -> (usize, Vec<usize>) {
    // (euler_totient, Z_m* i.e. prime relatives to m)
    todo!()
}

fn stride_permutation(prime: usize, input: &mut [Zq]) {
    let m = input.len();
    assert_eq!(m % prime, 0);
    let d = m / prime;

    let q = input.first().unwrap().get_mod();
    let zero = Zq::from_z_modulus(&Z::from(0), q);
    let mut temp = vec![zero; m];

    for i in 0..m - 1 {
        let new_index = (i * d) % (m - 1);
        temp[new_index] = input[i].clone(); // TODO: Remove clone
    }
    temp[m - 1] = input[m - 1].clone(); // TODO: Remove clone

    input.clone_from_slice(&temp);
}

fn inverse_stride_permutation(prime: usize, input: &mut [Zq]) {
    let m = input.len();
    assert_eq!(m % prime, 0);

    let q = input.first().unwrap().get_mod();
    let zero = Zq::from_z_modulus(&Z::from(0), q);
    let mut temp = vec![zero; m];

    for i in 0..m - 1 {
        let new_index = (i * prime) % (m - 1);
        temp[new_index] = input[i].clone(); // TODO: Remove clone
    }
    temp[m - 1] = input[m - 1].clone(); // TODO: Remove clone

    input.clone_from_slice(&temp);
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

fn icrt(evals: &mut [Zq], prime: usize, prime_power: usize, rou: Zq) {
    // Currently only using prime power decomposition, as we might noot need decomposition for all factors in NTT
    // get rou from ModulusPoly?
    let q = rou.get_mod().to_string().parse::<u32>().unwrap();
    let m = evals.len(); // Double check
    let m_prime = m / prime;
    let domain_size = prime.pow(prime_power as u32) as u32;
    let resized_power = q / domain_size;
    let omega = rou.pow(resized_power).unwrap();
    let mut current_omega_power = Zq::from_z_modulus(&Z::from(1), q);
    let mut omega_powers = Vec::new();
    for i in 0..prime.pow(prime_power as u32) {
        omega_powers.push(current_omega_power.clone());
        current_omega_power = &current_omega_power * &omega;
    }
    // Note that the first stride permutation is not done
    // This might not be needed because we will perform the inverse as well

    // Define crt_prime here to handle the prime case fast
    let prime_omegas = omega_powers
        .iter()
        .step_by(m_prime)
        .map(|w| w.clone())
        .collect::<Vec<_>>();

    let inv_crt_prime = crt_prime(&prime_omegas).inverse().unwrap();
    if prime_power == 1 {
        let evals_mat = MatZq::new(prime - 1, 1, q); // Theres got to be  abetter way or implment Mat * Vec multiplication
        let res = inv_crt_prime * evals_mat;
        for i in 0..prime - 1 {
            evals[i] = res.get_entry(i, 0).unwrap();
        }
        return;
    }

    // Perform inverse NTT_m_prime

    let mut m_prime_omega_powers = omega_powers
        .iter()
        .step_by(prime)
        .map(|w| w.clone())
        .collect::<Vec<_>>();
    m_prime_omega_powers.reverse();
    m_prime_omega_powers.rotate_right(1);

    let _ = evals.chunks_exact_mut(m_prime).map(|coeffs_chunk| {
        // TODO: Parallelize
        radixp_ntt(prime, prime_power - 1, &m_prime_omega_powers, coeffs_chunk);
        let nq_inv = Zq::from_z_modulus(&Z::from(m_prime_omega_powers.len() as u32), q)
            .inverse()
            .unwrap();
        for coeff in coeffs_chunk {
            *coeff = &*coeff * nq_inv.clone();
        }
    });

    let mut t_hat = twiddle_hat_factors(prime, &omega_powers);
    for twiddle in t_hat.iter_mut() {
        *twiddle = twiddle.inverse().unwrap();
    }
    for (twiddle_factor, coeff) in t_hat.iter().zip(evals.iter_mut()) {
        *coeff = &*coeff * twiddle_factor;
    }

    stride_permutation(prime, evals);

    // _ variables because it should modify coeffs
    let _ = evals // TODO: Parallelize
        .chunks_exact_mut(prime - 1) // phi(p) = p - 1
        .map(|coeffs_chunk| {
            let mut coeffs_mat = MatZq::new(prime, 1, q);
            for (i, coeff) in coeffs_chunk.iter().enumerate() {
                coeffs_mat.set_entry(i, 0, coeff).unwrap();
            }
            let res = inv_crt_prime.clone() * coeffs_mat;
            for i in 0..coeffs_chunk.len() {
                coeffs_chunk[i] = res.get_entry(i, 0).unwrap();
            }
        });

    inverse_stride_permutation(prime, evals);
}
