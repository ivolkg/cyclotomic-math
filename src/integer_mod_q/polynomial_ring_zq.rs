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

use std::{ops::Mul, str::FromStr};

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
    pub fn evaluate(&self, prime: usize, prime_power: usize, rou: Zq) -> PolynomialRingZqNTTBasis {
        let evals = self.ntt(prime, prime_power, rou);
        PolynomialRingZqNTTBasis { evals }
    }

    pub fn to_crt_basis(
        &self,
        prime: usize,
        prime_power: usize,
        rou: Zq,
    ) -> PolynomialRingZqCRTBasis {
        // Pass rou for m and not the general rou
        let poly = &self.poly;
        let q = rou.get_mod().to_string().parse::<u32>().unwrap();
        let euler_totient_m = euler_totient(prime, prime_power).count();

        let mut coeffs = (0..euler_totient_m) // Unclear if the coeficients are order from the least power to the greater, probably is
            .map(|i| {
                let coeff_z = poly.get_coeff(i).unwrap();
                Zq::from_z_modulus(&coeff_z, q)
            })
            .collect::<Vec<_>>();
        crt(&mut coeffs, prime, prime_power, rou);
        PolynomialRingZqCRTBasis { coeffs }
    }
    pub fn ntt(&self, prime: usize, prime_power: usize, rou: Zq) -> Vec<Zq> {
        // Currently only using prime power decomposition, as we might not need decomposition for all factors in NTT
        // get rou from ModulusPoly?
        let q = rou.get_mod().to_string().parse::<u32>().unwrap();
        let domain_size = prime.pow(prime_power as u32) as u32;
        let resized_power = q / domain_size;
        let omega = rou.pow(resized_power).unwrap();
        let mut current_omega_power = Zq::from_z_modulus(&Z::from(1), q);
        let mut omega_powers = Vec::new();
        for _ in 0..prime.pow(prime_power as u32) {
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

    fn from_vec(coeffs: Vec<Zq>, modulus: &ModulusPolynomialRingZq) -> Self {
        let mut poly_desc = format!("{} ", coeffs.len());
        for coeff in coeffs {
            poly_desc = poly_desc + &format!(" {}", coeff.value);
        }
        let poly = PolyOverZ::from_str(&poly_desc).unwrap();
        PolynomialRingZq::from((&poly, modulus))
    }
}

fn crt(coeffs: &mut [Zq], prime: usize, prime_power: usize, rou: Zq) {
    // Currently only using prime power decomposition, as we might noot need decomposition for all factors in NTT
    let q = rou.get_mod().to_string().parse::<u64>().unwrap();
    let varphi_m = coeffs.len(); // Double check
    let m_prime = varphi_m / (prime - 1);
    let domain_size = prime.pow(prime_power as u32) as u32;
    assert_eq!((q - 1) % domain_size as u64, 0);
    let resize_power = (q - 1) / domain_size as u64;
    let omega = rou.pow(resize_power).unwrap();
    let one = Zq::from_z_modulus(&Z::from(1), q);
    assert_eq!(
        omega.pow(domain_size).unwrap(),
        one,
        "{}-th root of unity not correct",
        domain_size
    );
    let mut current_omega_power = Zq::from_z_modulus(&Z::from(1), q);
    let mut omega_powers = Vec::new();
    for _ in 0..prime.pow(prime_power as u32) {
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

    if prime != 2 {
        stride_permutation(prime - 1, coeffs);

        for coeffs_chunk in coeffs.chunks_exact_mut(prime - 1) {
            let mut coeffs_mat = MatZq::new(prime - 1, 1, q);
            for (i, coeff) in coeffs_chunk.iter().enumerate() {
                coeffs_mat.set_entry(i, 0, coeff).unwrap();
            }
            let res = crt_prime.clone() * coeffs_mat;
            for i in 0..coeffs_chunk.len() {
                coeffs_chunk[i] = res.get_entry(i, 0).unwrap();
            }
        }

        inverse_stride_permutation(prime - 1, coeffs);
    }

    let t_hat = twiddle_hat_factors(prime, &omega_powers); // This could be omitted for powers of
                                                           // two

    for (twiddle_factor, coeff) in t_hat.iter().zip(coeffs.iter_mut()) {
        *coeff = &*coeff * twiddle_factor;
    }

    let m_prime_omega_powers = omega_powers
        .iter()
        .step_by(prime)
        .map(|w| w.clone())
        .collect::<Vec<_>>();

    for coeffs_chunk in coeffs.chunks_exact_mut(m_prime) {
        radixp_ntt(prime, prime_power - 1, &m_prime_omega_powers, coeffs_chunk);
    }

    assert_eq!(coeffs.len(), varphi_m);
    stride_permutation(prime - 1, coeffs);
}

// Double check
fn twiddle_hat_factors(prime: usize, omega_powers: &[Zq]) -> Vec<Zq> {
    let m_prime = omega_powers.len() / prime;
    let mut t_hat = Vec::new();
    for r in 1..prime {
        for j in 0..m_prime {
            t_hat.push(omega_powers[r * j].clone());
        }
    }
    t_hat
}

fn crt_prime(prime_omegas: &[Zq]) -> MatZq {
    let prime = prime_omegas.len();
    let euler_totient_p = prime - 1;
    let prime_relatives = 1..prime;
    assert_eq!(
        prime_omegas.len() - 1,
        euler_totient_p,
        "lenght is not a prime"
    );
    let q = prime_omegas.first().unwrap().get_mod();
    let mut crt = MatZq::new(euler_totient_p, euler_totient_p, q);
    for (i, rel) in prime_relatives.enumerate() {
        for j in 0..euler_totient_p {
            crt.set_entry(i, j, prime_omegas[(rel * j) % prime_omegas.len()].clone())
                .unwrap();
        }
    }
    crt
}

fn inverse_crt_prime(prime_omegas: &[Zq]) -> MatZq {
    crt_prime(prime_omegas).inverse().unwrap()
}

fn euler_totient(prime: usize, prime_power: usize) -> impl Iterator<Item = usize> {
    let m_prime = prime.pow(prime_power as u32-1);
    let relative_set = (1..prime)
        .chain(
            (1..(m_prime)).flat_map(move |i| i * prime + 1..prime*(i + 1))
        );
    relative_set
}

// Not necesary for m a prime power
// fn gcd(a: usize, b: usize) -> usize {
//     if b == 0 {
//         a
//     } else {
//         gcd(b, a % b)
//     }
// }

/// Note that this stride permutation will only be used for CRT and not DFT
/// so we use varphi instead or prime
fn inverse_stride_permutation(varphi_p: usize, input: &mut [Zq]) {
    let varphi_m = input.len();
    assert_eq!(varphi_m % varphi_p, 0);
    let d = varphi_m / varphi_p;

    let q = input.first().unwrap().get_mod();
    let zero = Zq::from_z_modulus(&Z::from(0), q);
    let mut temp = vec![zero; varphi_m];

    for i in 0..varphi_m - 1 {
        let new_index = (i * d) % (varphi_m - 1);
        temp[new_index] = input[i].clone(); // TODO: Remove clone
    }
    temp[varphi_m - 1] = input[varphi_m - 1].clone(); // TODO: Remove clone

    input.clone_from_slice(&temp);
}

fn stride_permutation(varphi_p: usize, input: &mut [Zq]) {
    let varphi_m = input.len();
    assert_eq!(varphi_m % varphi_p, 0);

    let q = input.first().unwrap().get_mod();
    let zero = Zq::from_z_modulus(&Z::from(0), q);
    let mut temp = vec![zero; varphi_m];

    for i in 0..varphi_m - 1 {
        let new_index = (i * varphi_p) % (varphi_m - 1);
        temp[new_index] = input[i].clone(); // TODO: Remove clone
    }
    temp[varphi_m - 1] = input[varphi_m - 1].clone(); // TODO: Remove clone

    input.clone_from_slice(&temp);
}

fn radixp_ntt(prime: usize, prime_power: usize, omega_powers: &[Zq], coeffs: &mut [Zq]) {
    let n = coeffs.len();
    if n == 1 {
        return;
    }
    assert_eq!(prime.pow(prime_power as u32), n);
    let q = omega_powers.first().unwrap().get_mod();
    let zero = Zq::from_z_modulus(&Z::from(0), q);
    let mut decomposed_coeffs = vec![vec![zero.clone(); n/prime]; prime];

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
    pub fn intt(&self, prime: usize, prime_power: usize, rou: Zq) -> Vec<Zq> {
        let q = rou.get_mod().to_string().parse::<u32>().unwrap(); //compute domain outside
        let domain_size = prime.pow(prime_power as u32) as u32;
        let resized_power = q / domain_size;
        let omega = rou.pow(resized_power).unwrap();
        let mut current_omega_power = Zq::from_z_modulus(&Z::from(1), q);
        let mut omega_powers = Vec::new();
        for _ in 0..prime.pow(prime_power as u32) {
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

impl PolynomialRingZqCRTBasis {}

impl PolynomialRingZqCRTBasis {
    pub fn to_powerful_basis(
        &self,
        prime: usize,
        prime_power: usize,
        rou: Zq,
        modulus: &ModulusPolynomialRingZq,
    ) -> PolynomialRingZq {
        let mut coeffs = self.coeffs.clone();
        icrt(&mut coeffs, prime, prime_power, rou);
        PolynomialRingZq::from_vec(coeffs, modulus)
    }
}

impl Mul for PolynomialRingZqCRTBasis {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let coeffs_1 = self.coeffs;
        let coeffs_2 = rhs.coeffs;
        let res_coeffs = coeffs_1
            .iter()
            .zip(coeffs_2.iter())
            .map(|(c1, c2)| c1 * c2)
            .collect::<Vec<_>>();
        PolynomialRingZqCRTBasis { coeffs: res_coeffs }
    }
}

fn icrt(crt_coeffs: &mut [Zq], prime: usize, prime_power: usize, rou: Zq) {
    // Currently only using prime power decomposition, as we might noot need decomposition for all factors in NTT
    let q = rou.get_mod().to_string().parse::<u64>().unwrap();
    let varphi_m = crt_coeffs.len(); // Double check
    let m_prime = varphi_m / (prime - 1);
    let domain_size = prime.pow(prime_power as u32) as u32;
    assert_eq!((q - 1) % domain_size as u64, 0);
    let resize_power = (q - 1) / domain_size as u64;
    let omega = rou.pow(resize_power).unwrap();
    let one = Zq::from_z_modulus(&Z::from(1), q);
    assert_eq!(
        omega.pow(domain_size).unwrap(),
        one,
        "{}-th root of unity not correct",
        domain_size
    );
    let mut current_omega_power = Zq::from_z_modulus(&Z::from(1), q);
    let mut omega_powers = Vec::new();
    for _ in 0..prime.pow(prime_power as u32) {
        omega_powers.push(current_omega_power.clone());
        current_omega_power = &current_omega_power * &omega;
    }
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
            crt_coeffs[i] = res.get_entry(i, 0).unwrap();
        }
        return;
    }

    inverse_stride_permutation(prime - 1, crt_coeffs);

    // Perform inverse NTT_m_prime

    let mut m_prime_omega_powers = omega_powers
        .iter()
        .step_by(prime)
        .map(|w| w.clone())
        .collect::<Vec<_>>();
    m_prime_omega_powers.reverse();
    m_prime_omega_powers.rotate_right(1);

    let nq_inv = Zq::from_z_modulus(&Z::from(m_prime_omega_powers.len() as u32), q)
        .inverse()
        .unwrap();
    for crt_coeffs_chunk in crt_coeffs.chunks_exact_mut(m_prime) {
        radixp_ntt(
            prime,
            prime_power - 1,
            &m_prime_omega_powers,
            crt_coeffs_chunk,
        );
        for coeff in crt_coeffs_chunk {
            *coeff = &*coeff * nq_inv.clone();
        }
    }

    let t_hat = twiddle_hat_factors(prime, &omega_powers);
    for (twiddle_factor, coeff) in t_hat.iter().zip(crt_coeffs.iter_mut()) {
        *coeff = &*coeff * twiddle_factor.inverse().unwrap();
    }

    if prime != 2 {
        stride_permutation(prime - 1, crt_coeffs);

        for crt_coeffs_chunk in crt_coeffs.chunks_exact_mut(prime - 1) {
            let mut coeffs_mat = MatZq::new(prime - 1, 1, q);
            for (i, coeff) in crt_coeffs_chunk.iter().enumerate() {
                coeffs_mat.set_entry(i, 0, coeff).unwrap();
            }
            let res = inv_crt_prime.clone() * coeffs_mat;
            for i in 0..crt_coeffs_chunk.len() {
                crt_coeffs_chunk[i] = res.get_entry(i, 0).unwrap();
            }
        }

        inverse_stride_permutation(prime - 1, crt_coeffs);
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    // use rand_distr::num_traits::Pow;

    use crate::traits::Pow;

    use super::*;

    fn modulus_poly_description(prime: usize, prime_power: usize, q: usize) -> String {
        let varphi_m = prime.pow(prime_power as u32 - 1) * (prime - 1);
        let mut coeffs = vec![0; varphi_m + 1];
        let m_prime = prime.pow(prime_power as u32 - 1);
        for i in 0..prime {
            coeffs[i * m_prime] = 1;
        }
        let mut desc = format!("{} ", varphi_m + 1);
        let coeffs_string = coeffs
            .iter()
            .map(|&x| x.to_string())
            .collect::<Vec<String>>()
            .join(" ");
        desc = format!("{} {} mod {}", desc, coeffs_string, q);
        desc
    }

    #[test]
    fn crt_po2_test() {
        let q = 15 * (1 << 27) + 1;
        let rou = Zq::from((76160998, q as u64));
        let prime = 2;
        let prime_power = 10;
        let mut m = 1;
        for _ in 0..prime_power {
            m = m * prime;
        }
        let mod_poly_desc = modulus_poly_description(prime, prime_power, q);
        let modulus_poly = ModulusPolynomialRingZq::from_str(&mod_poly_desc).unwrap();
        let poly = PolynomialRingZq::sample_uniform(modulus_poly.clone());
        let coeffs = (0..poly.get_degree() + 1)
            .map(|i| {
                let coeff = poly.get_coeff(i).unwrap();
                Zq::from_z_modulus(&coeff, q as u32)
            })
            .collect::<Vec<_>>();

        let new_rou = rou.pow(((q - 1) / m) as u32).unwrap();
        let relative_primes = euler_totient(prime, prime_power);
        let crt_coeffs_naive = relative_primes
            .map(|r| {
                let crt_coeff: Zq = coeffs.iter().enumerate().skip(1).fold(
                    coeffs.first().unwrap().clone(),
                    |acc, (i, c)| {
                        let x = c * &new_rou.pow((r * i) as u32).unwrap();
                        let acc = &x + acc;
                        acc
                    },
                );
                crt_coeff
            })
            .collect::<Vec<_>>();

        let poly_crt = poly.to_crt_basis(prime, prime_power, rou.clone());
        let crt_coeffs = poly_crt.coeffs;

        assert_eq!(crt_coeffs.len(), crt_coeffs_naive.len());
        for (i, (coeff, naive_coeff)) in crt_coeffs.iter().zip(crt_coeffs_naive.iter()).enumerate()
        {
            assert_eq!(coeff, naive_coeff, "coeff-{} not equal", i);
        }
    }

    #[test]
    fn crt_po3_test() {
        let q = (1 << 31) - 1; // Change field for a high 3-adicity field
        let rou = Zq::from((7, q as u64)); // Change Root of unity as well
        let prime = 3;
        let prime_power = 2;
        let mut m = 1;
        for _ in 0..prime_power {
            m = m * prime;
        }
        let mod_poly_desc = modulus_poly_description(prime, prime_power, q);
        let modulus_poly = ModulusPolynomialRingZq::from_str(&mod_poly_desc).unwrap();
        // let poly_z = PolyOverZ::from_str("6  1 2 3 4 5 6 7 8 9").unwrap();
        // let poly = PolynomialRingZq::from((&poly_z, &modulus_poly));
        let poly = PolynomialRingZq::sample_uniform(modulus_poly.clone());
        let coeffs = (0..poly.get_degree() + 1)
            .map(|i| {
                let coeff = poly.get_coeff(i).unwrap();
                Zq::from_z_modulus(&coeff, q as u32)
            })
            .collect::<Vec<_>>();

        let new_rou = rou.pow(((q - 1) / m) as u32).unwrap();
        let relative_primes = euler_totient(prime, prime_power);
        let crt_coeffs_naive = relative_primes
            .map(|r| {
                let crt_coeff: Zq = coeffs.iter().enumerate().skip(1).fold(
                    coeffs.first().unwrap().clone(),
                    |acc, (i, c)| {
                        let x = c * &new_rou.pow((r * i) as u32).unwrap();
                        let acc = &x + acc;
                        acc
                    },
                );
                crt_coeff
            })
            .collect::<Vec<_>>();

        let poly_crt = poly.to_crt_basis(prime, prime_power, rou.clone());
        let crt_coeffs = poly_crt.coeffs;

        assert_eq!(crt_coeffs.len(), crt_coeffs_naive.len());
        for (i, (coeff, naive_coeff)) in crt_coeffs.iter().zip(crt_coeffs_naive.iter()).enumerate()
        {
            assert_eq!(coeff, naive_coeff, "coeff-{} not equal", i);
        }
    }

    #[test]
    fn crt_mul_po2_test() {
        // Case for m = 2 not handled
        for prime_power in 4..5 {
            let q = 15 * (1 << 27) + 1;
            let rou = Zq::from((76160998, q as u64)); // Change the root of unity
            let prime = 2;
            let mod_poly_desc = modulus_poly_description(prime, prime_power, q);
            let modulus_poly = ModulusPolynomialRingZq::from_str(&mod_poly_desc).unwrap();
            let poly_1 = PolynomialRingZq::sample_uniform(modulus_poly.clone());
            let poly_2 = PolynomialRingZq::sample_uniform(modulus_poly.clone());

            let poly_1_crt = poly_1.to_crt_basis(prime, prime_power, rou.clone());

            let poly_2_crt = poly_2.to_crt_basis(prime, prime_power, rou.clone());

            let result_crt_poly = poly_1_crt * poly_2_crt;
            let result_poly =
                result_crt_poly.to_powerful_basis(prime, prime_power, rou, &modulus_poly);

            assert_eq!(
                result_poly,
                poly_1 * poly_2,
                "CRT multiplication failed at m = 2^{}",
                prime_power
            );
        }
    }
    #[test]
    fn crt_mul_po3_test() {
        let q = (1 << 31) - 1;
        let rou = Zq::from((7, q as u64)); // Change the root of unity
        let prime = 3;
        let prime_power = 2;
        let mod_poly_desc = modulus_poly_description(prime, prime_power, q);
        let modulus_poly = ModulusPolynomialRingZq::from_str(&mod_poly_desc).unwrap();
        let poly_1 = PolynomialRingZq::sample_uniform(modulus_poly.clone());
        let poly_2 = PolynomialRingZq::sample_uniform(modulus_poly.clone());

        let poly_1_crt = poly_1.to_crt_basis(prime, prime_power, rou.clone());

        let poly_2_crt = poly_2.to_crt_basis(prime, prime_power, rou.clone());

        let result_crt_poly = poly_1_crt * poly_2_crt;
        let result_poly = result_crt_poly.to_powerful_basis(prime, prime_power, rou, &modulus_poly);

        assert_eq!(result_poly, poly_1 * poly_2, "CRT multiplication failed");
    }
}
