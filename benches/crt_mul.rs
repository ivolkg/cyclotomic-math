use criterion::*;
use qfall_math::integer_mod_q::ModulusPolynomialRingZq;
use qfall_math::integer_mod_q::PolynomialRingZq;
use qfall_math::integer_mod_q::Zq;
use std::str::FromStr;

pub fn bench_crt_mul(c: &mut Criterion) {
    let q = 15 * (1 << 27) + 1;
    let rou = Zq::from((761609910, q as u64)); // Change the root of unity
    let prime = 2;
    let prime_power = 10;
    let mod_poly_desc = modulus_poly_description(prime, prime_power, q);
    let modulus_poly = ModulusPolynomialRingZq::from_str(&mod_poly_desc).unwrap();
    let poly_1 = PolynomialRingZq::sample_uniform(modulus_poly.clone());
    let poly_2 = PolynomialRingZq::sample_uniform(modulus_poly.clone());
    c.bench_function("CRT multiplication of 2^10-th cyclotomics", |b|{
        b.iter(||{
            let poly_1_crt = poly_1.to_crt_basis(prime, prime_power, rou.clone());

            let poly_2_crt = poly_2.to_crt_basis(prime, prime_power, rou.clone());

            let result_crt_poly = poly_1_crt * poly_2_crt;
            result_crt_poly.to_powerful_basis(prime, prime_power, rou.clone(), &modulus_poly);
        })
    });

}
pub fn bench_naive_mul(c: &mut Criterion) {
    let q = 15 * (1 << 27) + 1;
    let prime = 2;
    let prime_power = 10;
    let mod_poly_desc = modulus_poly_description(prime, prime_power, q);
    let modulus_poly = ModulusPolynomialRingZq::from_str(&mod_poly_desc).unwrap();
    let poly_1 = PolynomialRingZq::sample_uniform(modulus_poly.clone());
    let poly_2 = PolynomialRingZq::sample_uniform(modulus_poly.clone());
    c.bench_function("KS multiplication of 2^10-th cyclotomics", |b|{
        b.iter(||{
            &poly_1 * &poly_2
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
