use criterion::*;
use qfall_math::integer_mod_q::ModulusPolynomialRingZq;
use qfall_math::integer_mod_q::PolynomialRingZq;
use qfall_math::integer_mod_q::Zq;
use std::str::FromStr;

// // Mersenne31
// const Q: usize = (1 << 31) - 1;
// const ROU: usize = 7;
// const PRIME: usize = 3;
// const PPOWER: usize = 2;

// BabyBear
const Q: usize = 15 * (1 << 27) + 1;
const ROU: usize = 76160998;
const PRIME: usize = 2;
const PPOWER: usize = 8;
pub fn bench_crt_mul(c: &mut Criterion) {
    let rou = Zq::from((ROU as u32, Q as u64)); // Change the root of unity
    let mod_poly_desc = modulus_poly_description(PRIME, PPOWER, Q);
    let modulus_poly = ModulusPolynomialRingZq::from_str(&mod_poly_desc).unwrap();
    let poly_1 = PolynomialRingZq::sample_uniform(modulus_poly.clone());
    let poly_2 = PolynomialRingZq::sample_uniform(modulus_poly.clone());
    let desc = format!("CRT multiplication of {}^{}-th cyclotomics", PRIME, PPOWER);
    c.bench_function(desc.as_str(), |b| {
        b.iter(|| {
            let poly_1_crt = poly_1.to_crt_basis(PRIME, PPOWER, rou.clone());

            let poly_2_crt = poly_2.to_crt_basis(PRIME, PPOWER, rou.clone());

            let result_crt_poly = &poly_1_crt * &poly_2_crt;
                let result_crt_poly = &poly_1_crt * &poly_2_crt;
            result_crt_poly.to_powerful_basis(PRIME, PPOWER, rou.clone(), &modulus_poly);
        })
    });
}
pub fn bench_naive_mul(c: &mut Criterion) {
    let mod_poly_desc = modulus_poly_description(PRIME, PPOWER, Q);
    let modulus_poly = ModulusPolynomialRingZq::from_str(&mod_poly_desc).unwrap();
    let poly_1 = PolynomialRingZq::sample_uniform(modulus_poly.clone());
    let poly_2 = PolynomialRingZq::sample_uniform(modulus_poly.clone());
    let desc = format!("KS multiplication of {}^{}-th cyclotomics", PRIME, PPOWER);
    c.bench_function(desc.as_str(), |b| {
        b.iter(|| {
                &poly_1 * &poly_2;
        })
    });
}

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

criterion_group!(benches, bench_crt_mul, bench_naive_mul);
